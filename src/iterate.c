#include "include/iterate.h"
#include "Include/arm_math_types.h"
#include "Include/dsp/matrix_functions.h"
#include "include/laextension.h"
#include "include/quat.h"
#include "arm_math.h"
#include "math.h"
#include <string.h>

// Line-by-line port of the UKF in the Simulink MATLAB Function block (asbCubeSatACS.slx):
// iterate, ensure_positive_definite, calculate_lambda, get_sigma_points, get_weights,
// error_sigmas_to_quat_sigmas, propagate_quat_sigmas, get_measurements, get_error_vectors,
// gradient_descent. Each function below names the MATLAB function it mirrors.
//
// Conventions (same as the MATLAB):
//  - attitude errors are GLOBAL (ECI frame): q_sigma = dq(err) * q_est, err = rv(q_i * x^-1)
//  - propagation is in the body frame: q <- q * dq((gyro - bias) dt)
//  - Q is added after the unscented transform (P_hat = sum(...) + Q)
//
// The only intentional numerical deviation is quat2rotationvec (quat.c), which uses
// 2*atan2(|v|, w) instead of MATLAB's 2*acos(w). The two are identical in exact arithmetic,
// but acos(w) near w = 1 loses everything under ~1e-3 rad in float32.

#define QUAT_SIGMA_SIZE (STATE_SIZE + 1) // [q (4), bias (3)]
#define GRAD_DESCENT_MAX_ITERS 100       // gradient_descent: for iter = 1:100
#define GRAD_DESCENT_TOL 1e-9f           // gradient_descent: if norm(ave) < 1e-9
#define PD_EPSILON 1e-9f                 // ensure_positive_definite: epsilon = 1e-9
#define SIGMA_JITTER 1e-10f              // get_sigma_points: P + eye(n) * 1e-10
#define JACOBI_MAX_SWEEPS 50

/**
 * \fn cholesky_lower
 *
 * \brief L * L^T = A for a symmetric positive definite n x n A (L^T is MATLAB's chol(A)).
 * Zeroes the strict upper triangle of L, which CMSIS leaves untouched.
 */
static bool cholesky_lower(const float32_t* A, float32_t* L, int n) {
    memset(L, 0, sizeof(float32_t) * n * n);
    arm_matrix_instance_f32 A_mat = {(uint16_t)n, (uint16_t)n, (float32_t*)A};
    arm_matrix_instance_f32 L_mat = {(uint16_t)n, (uint16_t)n, L};
    if (arm_mat_cholesky_f32(&A_mat, &L_mat) != ARM_MATH_SUCCESS) {
        return false;
    }
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            L[i * n + j] = 0.0f;
        }
    }
    return all_finite(L, n * n);
}

/**
 * \fn sym_min_eig
 *
 * \brief Smallest eigenvalue of a symmetric n x n matrix (n <= STATE_SIZE), cyclic Jacobi.
 * Stands in for min(real(eig(P))); CMSIS has no eigen-solver.
 */
static float32_t sym_min_eig(const float32_t* P, int n) {
    float32_t A[STATE_SIZE * STATE_SIZE];
    memcpy(A, P, sizeof(float32_t) * n * n);
    for (int sweep = 0; sweep < JACOBI_MAX_SWEEPS; sweep++) {
        float32_t off = 0.0f, diag_sq = 0.0f;
        for (int i = 0; i < n; i++) {
            diag_sq += A[i * n + i] * A[i * n + i];
            for (int j = i + 1; j < n; j++) {
                off += A[i * n + j] * A[i * n + j];
            }
        }
        if (off <= 1e-14f * diag_sq || off < 1e-37f) {
            break;
        }
        for (int p = 0; p < n; p++) {
            for (int q = p + 1; q < n; q++) {
                float32_t apq = A[p * n + q];
                if (fabsf(apq) < 1e-38f) {
                    continue;
                }
                float32_t theta = (A[q * n + q] - A[p * n + p]) / (2.0f * apq);
                float32_t t = (theta >= 0.0f ? 1.0f : -1.0f) / (fabsf(theta) + sqrtf(theta * theta + 1.0f));
                float32_t c = 1.0f / sqrtf(t * t + 1.0f);
                float32_t s = t * c;
                for (int k = 0; k < n; k++) { // A <- A J (columns p, q)
                    float32_t akp = A[k * n + p], akq = A[k * n + q];
                    A[k * n + p] = c * akp - s * akq;
                    A[k * n + q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; k++) { // A <- J^T A (rows p, q)
                    float32_t apk = A[p * n + k], aqk = A[q * n + k];
                    A[p * n + k] = c * apk - s * aqk;
                    A[q * n + k] = s * apk + c * aqk;
                }
            }
        }
    }
    float32_t m = A[0];
    for (int i = 1; i < n; i++) {
        m = fminf(m, A[i * n + i]);
    }
    return m;
}

/**
 * \fn ensure_psd
 *
 * \brief Port of ensure_positive_definite(P) (epsilon = 1e-9), in place, n <= STATE_SIZE:
 *   P = (P + P') / 2;
 *   if min(eig(P)) < epsilon, P = P + eye(n) * (epsilon - min(eig(P)) + epsilon); end
 *
 * \return false if the matrix contains NaN/inf (MATLAB would propagate them)
 */
bool ensure_psd(float32_t* P, int n) {
    if (n <= 0 || n > STATE_SIZE || !all_finite(P, n * n)) {
        return false;
    }
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            float32_t avg = 0.5f * (P[i * n + j] + P[j * n + i]);
            P[i * n + j] = avg;
            P[j * n + i] = avg;
        }
    }
    float32_t min_eig = sym_min_eig(P, n);
    if (min_eig < PD_EPSILON) {
        float32_t shift = PD_EPSILON - min_eig + PD_EPSILON;
        for (int i = 0; i < n; i++) {
            P[i * n + i] += shift;
        }
    }
    return all_finite(P, n * n);
}

/**
 * \fn calculate_lambda
 *
 * \brief lambda = alpha^2 (n + kappa) - n
 */
float32_t calculate_lambda(int n, float32_t alpha, float32_t kappa) {
    return alpha * alpha * ((float32_t)n + kappa) - (float32_t)n;
}

/**
 * \fn get_sigma_points
 *
 * \brief Port of get_sigma_points(lam, x, P):
 *   P = (P + P')/2 + eye(n)*1e-10;  U = chol((lam + n) * P);
 *   sigmas = [x; x + U(k,:) for k = 1..n; x - U(k,:) for k = 1..n]
 * Rows of MATLAB's upper U are columns of the lower CMSIS factor L = U'.
 *
 * \return false where MATLAB's chol would error (matrix not positive definite)
 */
bool get_sigma_points(float32_t lam, const float32_t* state, const float32_t* P, float32_t* sigmas) {
    float32_t P1[STATE_SIZE * STATE_SIZE];
    for (int i = 0; i < STATE_SIZE; i++) {
        for (int j = 0; j < STATE_SIZE; j++) {
            float32_t sym = 0.5f * (P[i * STATE_SIZE + j] + P[j * STATE_SIZE + i]);
            if (i == j) {
                sym += SIGMA_JITTER;
            }
            P1[i * STATE_SIZE + j] = (lam + STATE_SIZE) * sym;
        }
    }

    float32_t L[STATE_SIZE * STATE_SIZE];
    if (!cholesky_lower(P1, L, STATE_SIZE)) {
        return false;
    }

    memcpy(sigmas, state, sizeof(float32_t) * STATE_SIZE);
    for (int k = 0; k < STATE_SIZE; k++) {
        float32_t* plus = &sigmas[(1 + k) * STATE_SIZE];
        float32_t* minus = &sigmas[(1 + STATE_SIZE + k) * STATE_SIZE];
        for (int i = 0; i < STATE_SIZE; i++) {
            float32_t offset = L[i * STATE_SIZE + k]; // column k of L == row k of U
            plus[i] = state[i] + offset;
            minus[i] = state[i] - offset;
        }
    }
    return true;
}

/**
 * \fn get_weights
 *
 * \brief Port of get_weights(lam, x, alpha, beta).
 */
void get_weights(float32_t lambda, int n, float32_t alpha, float32_t beta, float32_t* cov_weights,
                 float32_t* mean_weights) {
    float32_t c = 0.5f / ((float32_t)n + lambda);

    for (int i = 0; i < 2 * n + 1; i++) {
        cov_weights[i] = c;
        mean_weights[i] = c;
    }

    cov_weights[0] = lambda / ((float32_t)n + lambda) + (1 - alpha * alpha + beta);
    mean_weights[0] = lambda / ((float32_t)n + lambda);
}

/**
 * \fn error_sigmas_to_quat_sigmas
 *
 * \brief Port of error_sigmas_to_quat_sigmas: q_i = normalize(rv2q(err_i) * rotation).
 * Output rows are [q (4), bias (3)].
 */
static void error_sigmas_to_quat_sigmas(const float32_t* error_sigmas, const float32_t* rotation,
                                        float32_t* quat_sigmas) {
    for (int i = 0; i < NUM_SIGMAS; i++) {
        const float32_t* err = &error_sigmas[i * STATE_SIZE];
        float32_t* out = &quat_sigmas[i * QUAT_SIGMA_SIZE];
        float32_t quat_rot[4];
        rotationvec2quat(err, quat_rot);
        quat_multiply(quat_rot, rotation, out);
        quat_norm(out, out);
        out[4] = err[3];
        out[5] = err[4];
        out[6] = err[5];
    }
}

/**
 * \fn propagate_quat_sigmas
 *
 * \brief Port of propagate_quat_sigmas: q <- normalize(q * rv2q((gyro - bias) * dt)).
 */
static void propagate_quat_sigmas(float32_t* quat_sigmas, const float32_t* gyro, float32_t dt) {
    for (int i = 0; i < NUM_SIGMAS; i++) {
        float32_t* s = &quat_sigmas[i * QUAT_SIGMA_SIZE];
        float32_t w_dt[3] = {(gyro[0] - s[4]) * dt, (gyro[1] - s[5]) * dt, (gyro[2] - s[6]) * dt};
        float32_t delta_q[4];
        rotationvec2quat(w_dt, delta_q);
        quat_multiply(s, delta_q, s);
        quat_norm(s, s);
    }
}

/**
 * \fn get_measurements
 *
 * \brief Port of get_measurements: each ref vector rotated by quaternion_inverse(q_i).
 */
static void get_measurements(const float32_t* quat_sigmas, const float32_t* ref, int num_vecs, float32_t* msmts) {
    int m = 3 * num_vecs;
    for (int i = 0; i < NUM_SIGMAS; i++) {
        float32_t eci_to_body[4];
        quat_inv(&quat_sigmas[i * QUAT_SIGMA_SIZE], eci_to_body);
        for (int v = 0; v < num_vecs; v++) {
            quat_apply(eci_to_body, &ref[3 * v], &msmts[i * m + 3 * v]);
        }
    }
}

/**
 * \fn get_error_vectors
 *
 * \brief Port of get_error_vectors(Y, x): rv(Y_i * x^-1) (global error).
 */
static void get_error_vectors(const float32_t* quat_sigmas, const float32_t* x, float32_t* error_vecs) {
    float32_t x_inv[4];
    quat_inv(x, x_inv);
    for (int i = 0; i < NUM_SIGMAS; i++) {
        float32_t err_quat[4];
        quat_multiply(&quat_sigmas[i * QUAT_SIGMA_SIZE], x_inv, err_quat);
        quat2rotationvec(err_quat, &error_vecs[i * 3]);
    }
}

/**
 * \fn gradient_descent
 *
 * \brief Port of gradient_descent(Y, x, weights). Like the MATLAB, x is not renormalized
 * and, if all 100 iterations run, error_vecs belong to the second-to-last x.
 */
static void gradient_descent(const float32_t* quat_sigmas, const float32_t* x0, const float32_t* weights,
                             float32_t* x, float32_t* error_vecs) {
    float32_t w_sum = 0.0f;
    for (int i = 0; i < NUM_SIGMAS; i++) {
        w_sum += weights[i];
    }
    float32_t w[NUM_SIGMAS];
    for (int i = 0; i < NUM_SIGMAS; i++) {
        w[i] = weights[i] / w_sum;
    }

    memcpy(x, x0, sizeof(float32_t) * 4);
    for (int iter = 0; iter < GRAD_DESCENT_MAX_ITERS; iter++) {
        get_error_vectors(quat_sigmas, x, error_vecs);
        float32_t ave[3] = {0, 0, 0};
        for (int i = 0; i < NUM_SIGMAS; i++) {
            for (int j = 0; j < 3; j++) {
                ave[j] += error_vecs[i * 3 + j] * w[i];
            }
        }
        if (l2_norm(ave, 3) < GRAD_DESCENT_TOL) {
            break;
        }
        float32_t ave_quat[4];
        rotationvec2quat(ave, ave_quat);
        quat_multiply(ave_quat, x, x);
    }
}

ukf_status_t iterate(const float* error_state, const float* quat_state, const float* cov,
                     const float* body, const float* ref, int num_vecs, const float* gyro,
                     const float* Q, const float* R, float dt, float* new_err_state,
                     float* new_quat_state, float* new_P) {
    // Copy inputs first so outputs may alias them, and so failures can return the inputs
    float x_in[STATE_SIZE];
    float q_in[4];
    float P_in[STATE_SIZE * STATE_SIZE];
    memcpy(x_in, error_state, sizeof(x_in));
    memcpy(q_in, quat_state, sizeof(q_in));
    memcpy(P_in, cov, sizeof(P_in));

#define UKF_FAIL(code)                                                                             \
    do {                                                                                           \
        memcpy(new_err_state, x_in, sizeof(x_in));                                                 \
        memcpy(new_quat_state, q_in, sizeof(q_in));                                                \
        memcpy(new_P, P_in, sizeof(P_in));                                                         \
        return (code);                                                                             \
    } while (0)

    if (num_vecs < 0 || num_vecs > MAX_MSMT_VECS || !(dt >= 0.0f) || !all_finite(x_in, STATE_SIZE)
        || !all_finite(q_in, 4) || !all_finite(gyro, 3) || !all_finite(Q, STATE_SIZE * STATE_SIZE)) {
        UKF_FAIL(UKF_ERR_ARGS);
    }
    const int m = 3 * num_vecs;
    if (num_vecs > 0 && (!all_finite(body, m) || !all_finite(ref, m) || !all_finite(R, m * m))) {
        UKF_FAIL(UKF_ERR_ARGS);
    }

    // current_cov = ensure_positive_definite(current_cov);
    float current_cov[STATE_SIZE * STATE_SIZE];
    memcpy(current_cov, P_in, sizeof(current_cov));
    if (!ensure_psd(current_cov, STATE_SIZE)) {
        UKF_FAIL(UKF_ERR_CHOLESKY);
    }

    float lam = calculate_lambda(STATE_SIZE, UKF_ALPHA, UKF_KAPPA);
    float sigmas[NUM_SIGMAS * STATE_SIZE];
    if (!get_sigma_points(lam, x_in, current_cov, sigmas)) {
        UKF_FAIL(UKF_ERR_CHOLESKY);
    }
    float quat_sigmas[NUM_SIGMAS * QUAT_SIGMA_SIZE];
    error_sigmas_to_quat_sigmas(sigmas, q_in, quat_sigmas);
    propagate_quat_sigmas(quat_sigmas, gyro, dt); // always propagate

    float covariance_weights[NUM_SIGMAS];
    float mean_weights[NUM_SIGMAS];
    get_weights(lam, STATE_SIZE, UKF_ALPHA, UKF_BETA, covariance_weights, mean_weights);

    // [average_quaternion, propagated_error_vectors] = gradient_descent(..., current_guess, mean_weights)
    float average_quaternion[4];
    float err_vecs[NUM_SIGMAS * 3];
    gradient_descent(quat_sigmas, q_in, mean_weights, average_quaternion, err_vecs);

    float propagated_errors[NUM_SIGMAS * STATE_SIZE];
    for (int i = 0; i < NUM_SIGMAS; i++) {
        for (int j = 0; j < 3; j++) {
            propagated_errors[i * STATE_SIZE + j] = err_vecs[i * 3 + j];
            propagated_errors[i * STATE_SIZE + 3 + j] = quat_sigmas[i * QUAT_SIGMA_SIZE + 4 + j];
        }
    }

    float mean_error[STATE_SIZE] = {0};
    for (int col = 0; col < STATE_SIZE; col++) {
        for (int row = 0; row < NUM_SIGMAS; row++) {
            mean_error[col] += propagated_errors[row * STATE_SIZE + col] * mean_weights[row];
        }
    }

    float dx[NUM_SIGMAS * STATE_SIZE];
    for (int i = 0; i < NUM_SIGMAS; i++) {
        for (int j = 0; j < STATE_SIZE; j++) {
            dx[i * STATE_SIZE + j] = propagated_errors[i * STATE_SIZE + j] - mean_error[j];
        }
    }

    float P_hat[STATE_SIZE * STATE_SIZE] = {0};
    for (int i = 0; i < NUM_SIGMAS; i++) {
        for (int r = 0; r < STATE_SIZE; r++) {
            for (int c = 0; c < STATE_SIZE; c++) {
                P_hat[r * STATE_SIZE + c] += covariance_weights[i] * dx[i * STATE_SIZE + r] * dx[i * STATE_SIZE + c];
            }
        }
    }
    for (int i = 0; i < STATE_SIZE * STATE_SIZE; i++) {
        P_hat[i] += Q[i]; // process noise always accrues, every propagate step
    }

    float x_hat[STATE_SIZE];
    float P[STATE_SIZE * STATE_SIZE];

    if (num_vecs > 0) {
        // --- MEASUREMENT CORRECTION ---
        float measurements[NUM_SIGMAS * MAX_MSMT_SIZE];
        get_measurements(quat_sigmas, ref, num_vecs, measurements);

        float mean_measurement[MAX_MSMT_SIZE] = {0};
        for (int col = 0; col < m; col++) {
            for (int row = 0; row < NUM_SIGMAS; row++) {
                mean_measurement[col] += measurements[row * m + col] * mean_weights[row];
            }
        }

        float dz[NUM_SIGMAS * MAX_MSMT_SIZE];
        for (int i = 0; i < NUM_SIGMAS; i++) {
            for (int j = 0; j < m; j++) {
                dz[i * m + j] = measurements[i * m + j] - mean_measurement[j];
            }
        }

        float P_xz[STATE_SIZE * MAX_MSMT_SIZE] = {0};
        float P_zz[MAX_MSMT_SIZE * MAX_MSMT_SIZE] = {0};
        for (int i = 0; i < NUM_SIGMAS; i++) {
            for (int r = 0; r < STATE_SIZE; r++) {
                for (int c = 0; c < m; c++) {
                    P_xz[r * m + c] += covariance_weights[i] * dx[i * STATE_SIZE + r] * dz[i * m + c];
                }
            }
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < m; c++) {
                    P_zz[r * m + c] += covariance_weights[i] * dz[i * m + r] * dz[i * m + c];
                }
            }
        }

        // P_vv = P_zz + R
        float P_vv[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
        for (int i = 0; i < m * m; i++) {
            P_vv[i] = P_zz[i] + R[i];
        }

        // k = P_xz * inv(P_vv). arm_mat_inverse_f32 destroys its input, so invert a copy.
        float P_vv_scratch[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
        float P_vv_inv[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
        memcpy(P_vv_scratch, P_vv, sizeof(float) * m * m);
        arm_matrix_instance_f32 P_vv_scratch_mat = {(uint16_t)m, (uint16_t)m, P_vv_scratch};
        arm_matrix_instance_f32 P_vv_inv_mat = {(uint16_t)m, (uint16_t)m, P_vv_inv};
        if (arm_mat_inverse_f32(&P_vv_scratch_mat, &P_vv_inv_mat) != ARM_MATH_SUCCESS
            || !all_finite(P_vv_inv, m * m)) {
            UKF_FAIL(UKF_ERR_SINGULAR);
        }
        float K[STATE_SIZE * MAX_MSMT_SIZE];
        arm_matrix_instance_f32 P_xz_mat = {STATE_SIZE, (uint16_t)m, P_xz};
        arm_matrix_instance_f32 K_mat = {STATE_SIZE, (uint16_t)m, K};
        arm_mat_mult_f32(&P_xz_mat, &P_vv_inv_mat, &K_mat);

        // x_hat = mean_error + (k * (body_msmts - mean_measurement)')'
        for (int r = 0; r < STATE_SIZE; r++) {
            float corr = 0.0f;
            for (int c = 0; c < m; c++) {
                corr += K[r * m + c] * (body[c] - mean_measurement[c]);
            }
            x_hat[r] = mean_error[r] + corr;
        }

        // P = P_hat - k * P_vv * k'
        float KPvv[STATE_SIZE * MAX_MSMT_SIZE];
        for (int r = 0; r < STATE_SIZE; r++) {
            for (int c = 0; c < m; c++) {
                float acc = 0.0f;
                for (int j = 0; j < m; j++) {
                    acc += K[r * m + j] * P_vv[j * m + c];
                }
                KPvv[r * m + c] = acc;
            }
        }
        for (int r = 0; r < STATE_SIZE; r++) {
            for (int c = 0; c < STATE_SIZE; c++) {
                float acc = 0.0f;
                for (int j = 0; j < m; j++) {
                    acc += KPvv[r * m + j] * K[c * m + j];
                }
                P[r * STATE_SIZE + c] = P_hat[r * STATE_SIZE + c] - acc;
            }
        }
    } else {
        // --- PREDICT-ONLY ---
        memcpy(x_hat, mean_error, sizeof(x_hat));
        memcpy(P, P_hat, sizeof(P));
    }

    // P = ensure_positive_definite(P);
    if (!ensure_psd(P, STATE_SIZE)) {
        UKF_FAIL(UKF_ERR_CHOLESKY);
    }

    // new_guess = rotation_vec2quaternion(x_hat(1:3)) * average_quaternion (not renormalized, as in MATLAB)
    float x_hat_rot[4];
    rotationvec2quat(x_hat, x_hat_rot);
    float new_guess[4];
    quat_multiply(x_hat_rot, average_quaternion, new_guess);

    // The Simulink block resets error_state(1:3) = 0 right after iterate(); done here instead
    x_hat[0] = 0.0f;
    x_hat[1] = 0.0f;
    x_hat[2] = 0.0f;

    if (!all_finite(x_hat, STATE_SIZE) || !all_finite(new_guess, 4) || !all_finite(P, STATE_SIZE * STATE_SIZE)) {
        UKF_FAIL(UKF_ERR_NONFINITE);
    }

    memcpy(new_err_state, x_hat, sizeof(x_hat));
    memcpy(new_quat_state, new_guess, sizeof(new_guess));
    memcpy(new_P, P, sizeof(P));
    return UKF_OK;
#undef UKF_FAIL
}
