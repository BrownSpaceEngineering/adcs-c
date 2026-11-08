#ifndef ITER
#define ITER
#include "arm_math.h"
#include <stdbool.h>

// Multiplicative UKF on [attitude error rotation vector (3), gyro bias (3)].
// Exact port of the Simulink block UKF. Attitude errors are GLOBAL (ECI frame):
// q_sigma = dq(err) * q_est; propagation is q * dq((gyro - bias) dt).
#define STATE_SIZE 6
#define MAX_MSMT_VECS 2
#define MAX_MSMT_SIZE (3 * MAX_MSMT_VECS)
#define MSMT_SIZE MAX_MSMT_SIZE
#define NUM_SIGMAS (STATE_SIZE * 2 + 1)

// Unscented transform parameters. These match the latest MATLAB (simulink.m): alpha = 1,
// kappa = 0. The older alpha = 0.01, kappa = 3 - n setting (from the double-precision Python
// prototype) gives a centre weight of about -2e4, which cancels catastrophically in float32.
// Overridable at compile time (e.g. -DUKF_ALPHA=0.5f) for tuning experiments.
#ifndef UKF_ALPHA
#define UKF_ALPHA 1.0f
#endif
#ifndef UKF_BETA
#define UKF_BETA 2.0f
#endif
#ifndef UKF_KAPPA
#define UKF_KAPPA 0.0f
#endif

typedef enum {
    UKF_OK = 0,
    UKF_ERR_ARGS,       // bad num_vecs / dt / non-finite inputs
    UKF_ERR_CHOLESKY,   // covariance could not be made positive definite
    UKF_ERR_SINGULAR,   // innovation covariance was singular
    UKF_ERR_NONFINITE,  // NaN/inf in the result
} ukf_status_t;

/**
 * One predict (+ optional update) step. Estimates body->ref.
 *
 * error_state  [6]      attitude part must be 0 (it is folded into quat_state every step);
 *                       bias part is the current gyro bias estimate (rad/s)
 * quat_state   [4]      current body->ref estimate
 * cov          [6x6]    error covariance
 * body, ref    [3*num_vecs] measured body-frame vectors and the matching ref-frame vectors
 * num_vecs     0 = predict only, 1 or 2 vectors otherwise
 * gyro         [3]      raw gyro measurement (rad/s, body frame)
 * Q            [6x6]    process noise added to the propagated covariance (P_hat + Q)
 * R            [(3*num_vecs)^2] measurement noise (unused when num_vecs == 0)
 *
 * Outputs may alias inputs. On failure outputs are copies of the inputs.
 * new_err_state has its attitude part zeroed (already applied to new_quat_state).
 */
ukf_status_t iterate(const float* error_state, const float* quat_state, const float* cov,
                     const float* body, const float* ref, int num_vecs, const float* gyro,
                     const float* Q, const float* R, float dt, float* new_err_state,
                     float* new_quat_state, float* new_P);

// Internal steps, exposed for unit testing
float32_t calculate_lambda(int n, float32_t alpha, float32_t kappa);
void get_weights(float32_t lambda, int n, float32_t alpha, float32_t beta, float32_t* cov_weights,
                 float32_t* mean_weights);
// sigmas is NUM_SIGMAS rows x STATE_SIZE. Returns false if P can't be factored.
bool get_sigma_points(float32_t lam, const float32_t* state, const float32_t* P, float32_t* sigmas);
// Port of ensure_positive_definite: symmetrize P (n x n, n <= STATE_SIZE); if min eig < 1e-9 shift the diagonal.
bool ensure_psd(float32_t* P, int n);
#endif
