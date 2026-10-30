#include "include/test.h"
#include "Include/dsp/basic_math_functions.h"
#include "Include/dsp/matrix_functions.h"
#include "arm_math.h"
#include "include/bdot.h"
#include "include/body.h"
#include "include/down_quat.h"
#include "include/ecef2eci.h"
#include "include/filters.h"
#include "include/iterate.h"
#include "include/kepler.h"
#include "include/laextension.h"
#include "include/magnetosphere.h"
#include "include/moments2currents.h"
#include "include/orbital2eci.h"
#include "include/pd.h"
#include "include/photodiode_determination.h"
#include "include/quat.h"
#include "include/quest.h"
#include "include/sgp4_step.h"
#include "include/sunvec.h"
#include "include/torque2moments.h"
#include "include/matlab_reference.h"
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEG2RAD(x) ((float)(x) * (float)M_PI / 180.0f)
#define RAD2DEG(x) ((float)(x) * 180.0f / (float)M_PI)

/*
 * Tiny test harness. Every check prints a PASS/FAIL line and is counted; nothing ever aborts,
 * so one failure never hides the rest of the results. Expected values marked "python ref" come
 * from tools/reference_values.py (independent double-precision implementations).
 */

static int checks_run = 0;
static int checks_failed = 0;
static int notes = 0;
static const char* current_suite = "";

static void begin_suite(const char* name) {
    current_suite = name;
    printf("\n===== %s =====\n", name);
}

// Prints "[PASS]/[FAIL] suite: name -- detail". Details are always shown for failures and,
// when show_on_pass is set (numeric checks), for passes too so margins are visible.
static bool vreport(bool ok, bool show_on_pass, const char* name, const char* detail_fmt, va_list args) {
    checks_run++;
    if (!ok) {
        checks_failed++;
    }
    printf("[%s] %s: %s", ok ? "PASS" : "FAIL", current_suite, name);
    if (detail_fmt != NULL && (!ok || show_on_pass)) {
        printf(" -- ");
        vprintf(detail_fmt, args);
    }
    printf("\n");
    return ok;
}

// Informational line, not a pass/fail (used for known bugs on the MATLAB side)
static void note(const char* fmt, ...) {
    notes++;
    printf("[NOTE] %s: ", current_suite);
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    printf("\n");
}

static bool report(bool ok, const char* name, const char* detail_fmt, ...) {
    va_list args;
    va_start(args, detail_fmt);
    bool r = vreport(ok, false, name, detail_fmt, args);
    va_end(args);
    return r;
}

static bool report_value(bool ok, const char* name, const char* detail_fmt, ...) {
    va_list args;
    va_start(args, detail_fmt);
    bool r = vreport(ok, true, name, detail_fmt, args);
    va_end(args);
    return r;
}

static bool check_true(bool cond, const char* name) { return report(cond, name, NULL); }

// Exact float equality without tripping -Wfloat-equal (intentional in a few tests)
static bool exactly(float a, float b) { return !(a < b) && !(a > b); }

static bool check_close(const char* name, double got, double expected, double tol) {
    bool ok = isfinite(got) && fabs(got - expected) <= tol;
    return report_value(ok, name, "got %.9g, expected %.9g (tol %.3g)", got, expected, tol);
}

static bool check_less(const char* name, double got, double limit) {
    bool ok = isfinite(got) && got < limit;
    return report_value(ok, name, "got %.6g, limit %.6g", got, limit);
}

static bool check_vec_close(const char* name, const float* got, const float* expected, int n, float tol) {
    bool ok = true;
    for (int i = 0; i < n; i++) {
        if (!isfinite(got[i]) || fabsf(got[i] - expected[i]) > tol) {
            ok = false;
        }
    }
    if (ok) {
        return report(true, name, NULL);
    }
    char buf[512];
    int off = snprintf(buf, sizeof(buf), "got [");
    for (int i = 0; i < n && off < (int)sizeof(buf); i++) {
        off += snprintf(buf + off, sizeof(buf) - off, "%s%.7g", i ? ", " : "", got[i]);
    }
    off += snprintf(buf + off, sizeof(buf) - off, "], expected [");
    for (int i = 0; i < n && off < (int)sizeof(buf); i++) {
        off += snprintf(buf + off, sizeof(buf) - off, "%s%.7g", i ? ", " : "", expected[i]);
    }
    snprintf(buf + off, sizeof(buf) - off, "] (tol %.3g)", tol);
    return report(false, name, "%s", buf);
}

// ---------------------------------------------------------------------------------------------
// Deterministic RNG (xorshift32) so results are identical on every platform, unlike rand()
// ---------------------------------------------------------------------------------------------
static uint32_t rng_state = 0x12345678u;

static void rng_seed(uint32_t seed) { rng_state = seed ? seed : 0x12345678u; }

static uint32_t rng_u32(void) {
    uint32_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_state = x;
    return x;
}

// Uniform in (0, 1)
static float uniform_01(void) { return ((float)(rng_u32() >> 8) + 0.5f) / 16777216.0f; }

static float uniform_pm1(void) { return 2.0f * uniform_01() - 1.0f; }

static float normal_sample(float mean, float stddev) {
    float u1 = uniform_01();
    float u2 = uniform_01();
    float mag = sqrtf(-2.0f * logf(u1));
    float z0 = mag * cosf(2.0f * (float)M_PI * u2);
    return mean + stddev * z0;
}

static void random_unit_vec(float* v) {
    do {
        v[0] = normal_sample(0, 1);
        v[1] = normal_sample(0, 1);
        v[2] = normal_sample(0, 1);
    } while (!normalize_vec(v, 3));
}

static void random_quat(float* q) {
    for (int i = 0; i < 4; i++) {
        q[i] = normal_sample(0, 1);
    }
    quat_norm(q, q);
}

// ---------------------------------------------------------------------------------------------
// Independent (double precision, explicit formula) helpers used to check the library
// ---------------------------------------------------------------------------------------------

// R(q) v for a Hamilton scalar-first quaternion (active rotation)
static void ref_rotate(const float* qf, const float* vf, float* out) {
    double n = sqrt((double)qf[0] * qf[0] + (double)qf[1] * qf[1] + (double)qf[2] * qf[2]
                    + (double)qf[3] * qf[3]);
    double w = qf[0] / n, x = qf[1] / n, y = qf[2] / n, z = qf[3] / n;
    double R[9] = {1 - 2 * (y * y + z * z), 2 * (x * y - w * z),     2 * (x * z + w * y),
                   2 * (x * y + w * z),     1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
                   2 * (x * z - w * y),     2 * (y * z + w * x),     1 - 2 * (x * x + y * y)};
    for (int i = 0; i < 3; i++) {
        out[i] = (float)(R[i * 3] * vf[0] + R[i * 3 + 1] * vf[1] + R[i * 3 + 2] * vf[2]);
    }
}

// R(q)^T v
static void ref_rotate_inv(const float* q, const float* v, float* out) {
    float qc[4] = {q[0], -q[1], -q[2], -q[3]};
    ref_rotate(qc, v, out);
}

// Rotation angle (rad) between two attitudes, insensitive to the q / -q sign ambiguity
static float quat_angle_between(const float* a, const float* b) {
    double na = sqrt((double)a[0] * a[0] + (double)a[1] * a[1] + (double)a[2] * a[2] + (double)a[3] * a[3]);
    double nb = sqrt((double)b[0] * b[0] + (double)b[1] * b[1] + (double)b[2] * b[2] + (double)b[3] * b[3]);
    double d = fabs(((double)a[0] * b[0] + (double)a[1] * b[1] + (double)a[2] * b[2] + (double)a[3] * b[3]) / (na * nb));
    if (d > 1.0) {
        d = 1.0;
    }
    // 2*acos(d), computed as 2*atan2 for accuracy near 0
    return (float)(2.0 * atan2(sqrt(1.0 - d * d), d));
}

static float vec_angle(const float* a, const float* b) {
    float c[3];
    cross(a, b, c);
    return atan2f(l2_norm(c, 3), dot3(a, b));
}

static bool eps_close_matrix(const float32_t* A, const float32_t* B, int rows, int cols, float32_t eps) {
    for (int i = 0; i < rows * cols; i++) {
        if (!(fabsf(A[i] - B[i]) <= eps)) {
            return false;
        }
    }
    return true;
}

static bool matrix_is_symmetric(const float* P, int n, float tol) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (!(fabsf(P[i * n + j] - P[j * n + i]) <= tol)) {
                return false;
            }
        }
    }
    return true;
}

// Positive definite iff a Cholesky factorization succeeds
static bool matrix_is_pd(const float* P, int n) {
    float A[36], L[36];
    memcpy(A, P, sizeof(float) * n * n);
    memset(L, 0, sizeof(L));
    arm_matrix_instance_f32 A_mat = {(uint16_t)n, (uint16_t)n, A};
    arm_matrix_instance_f32 L_mat = {(uint16_t)n, (uint16_t)n, L};
    return arm_mat_cholesky_f32(&A_mat, &L_mat) == ARM_MATH_SUCCESS;
}

static void diag(float* A, int n, const float* values) {
    memset(A, 0, sizeof(float) * n * n);
    for (int i = 0; i < n; i++) {
        A[i * n + i] = values[i];
    }
}

// =============================================================================================

int test_run_all(void) {
    checks_run = 0;
    checks_failed = 0;
    notes = 0;

    test_linalg();
    test_matrix_product();
    test_quaternion();
    test_quest();
    test_ukf_internals();
    test_ukf_behaviour();
    test_iteration_1vec();
    test_iteration_2vec();
    test_bdot();
    test_pd();
    test_torque_and_currents();
    test_kepler();
    test_orbital_to_eci();
    test_ecef_to_eci();
    test_sun_vec();
    test_magnetosphere();
    test_photodiodes();
    test_pointing();
    test_filters();
    test_body();
    test_matlab_parity();

    printf("\n===== SUMMARY =====\n");
    printf("%d checks, %d passed, %d failed (%d notes)\n", checks_run, checks_run - checks_failed, checks_failed,
           notes);
    return checks_failed;
}

// =============================================================================================
// Linear algebra helpers
// =============================================================================================
void test_linalg(void) {
    begin_suite("linalg");

    float ex[3] = {1, 0, 0}, ey[3] = {0, 1, 0}, ez[3] = {0, 0, 1};
    float c[3];
    cross(ex, ey, c);
    check_vec_close("x cross y = z", c, ez, 3, 0);
    cross(ey, ez, c);
    check_vec_close("y cross z = x", c, ex, 3, 0);
    cross(ez, ex, c);
    check_vec_close("z cross x = y", c, ey, 3, 0);

    float a[3] = {1.5f, -2.0f, 0.25f}, b[3] = {-3.0f, 0.5f, 4.0f};
    float ab[3], ba[3];
    cross(a, b, ab);
    cross(b, a, ba);
    float neg_ba[3] = {-ba[0], -ba[1], -ba[2]};
    check_vec_close("cross is anti-commutative", ab, neg_ba, 3, 1e-6f);
    check_close("cross result is orthogonal to a", dot3(ab, a), 0, 1e-5);
    check_close("cross result is orthogonal to b", dot3(ab, b), 0, 1e-5);
    float expected_ab[3] = {-2.0f * 4.0f - 0.25f * 0.5f, 0.25f * -3.0f - 1.5f * 4.0f, 1.5f * 0.5f - (-2.0f * -3.0f)};
    check_vec_close("cross product value", ab, expected_ab, 3, 1e-6f);
    float aliased[3] = {a[0], a[1], a[2]};
    cross(aliased, b, aliased);
    check_vec_close("cross output may alias an input", aliased, ab, 3, 1e-6f);

    check_close("dot3", dot3(a, b), 1.5 * -3.0 + -2.0 * 0.5 + 0.25 * 4.0, 1e-6);

    float arr[5] = {3.0f, -7.5f, 12.25f, 0.0f, 4.0f};
    check_close("max_arr", max_arr(arr, 5), 12.25, 0);
    check_close("min_arr", min_arr(arr, 5), -7.5, 0);
    check_close("max_arr of empty array is 0", max_arr(arr, 0), 0, 0);
    check_close("min_arr of empty array is 0", min_arr(arr, 0), 0, 0);

    float v345[3] = {3, 4, 12};
    check_close("l2_norm(3,4,12) = 13", l2_norm(v345, 3), 13, 1e-6);
    float zero[3] = {0, 0, 0};
    check_close("l2_norm(0) = 0", l2_norm(zero, 3), 0, 0);

    float to_norm[3] = {3, 4, 12};
    bool ok = normalize_vec(to_norm, 3);
    float expected_norm[3] = {3.0f / 13, 4.0f / 13, 12.0f / 13};
    check_true(ok, "normalize_vec succeeds on non-zero vector");
    check_vec_close("normalize_vec value", to_norm, expected_norm, 3, 1e-7f);
    float z2[3] = {0, 0, 0};
    check_true(!normalize_vec(z2, 3) && exactly(z2[0], 0) && exactly(z2[1], 0) && exactly(z2[2], 0),
               "normalize_vec refuses the zero vector and leaves it untouched");
    float nanv[3] = {NAN, 1, 1};
    check_true(!normalize_vec(nanv, 3), "normalize_vec refuses NaN");

    float finite[3] = {1, 2, 3}, inf[3] = {1, INFINITY, 3};
    check_true(all_finite(finite, 3), "all_finite on finite data");
    check_true(!all_finite(inf, 3), "all_finite detects inf");
    check_true(!all_finite(nanv, 3), "all_finite detects NaN");

    float I4[16];
    eye_f32(I4, 4);
    bool eye_ok = true;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (!exactly(I4[i * 4 + j], i == j ? 1.0f : 0.0f)) {
                eye_ok = false;
            }
        }
    }
    check_true(eye_ok, "eye(4) is the identity");

    float M[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    arm_matrix_instance_f32 M_mat = {3, 3, M};
    check_close("trace of 3x3", trace(&M_mat), 15, 0);
    arm_matrix_instance_f32 M_rect = {3, 2, M};
    check_close("trace of non-square matrix is 0", trace(&M_rect), 0, 0);
}

void test_matrix_product(void) {
    begin_suite("matrix product (CMSIS)");

    float A[4] = {1., 2., 3., 4.};
    float B[4] = {5., 6., 7., 8.};
    float C[4] = {0.};
    float C_expected[4] = {19., 22., 43., 50.};

    arm_matrix_instance_f32 A_mat = {2, 2, A};
    arm_matrix_instance_f32 B_mat = {2, 2, B};
    arm_matrix_instance_f32 C_mat = {2, 2, C};
    arm_mat_mult_f32(&A_mat, &B_mat, &C_mat);
    check_true(eps_close_matrix(C, C_expected, 2, 2, FLT_EPSILON), "2x2 product");

    float identity[4 * 4] = {0.};
    float large_result[4 * 4] = {0.};
    eye_f32(identity, 4);
    arm_matrix_instance_f32 identity_mat = {4, 4, identity};
    arm_matrix_instance_f32 large_result_mat = {4, 4, large_result};
    arm_mat_mult_f32(&identity_mat, &identity_mat, &large_result_mat);
    check_true(eps_close_matrix(identity, large_result, 4, 4, FLT_EPSILON), "identity squared");

    float A_large[4 * 4] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.};
    arm_matrix_instance_f32 A_large_mat = {4, 4, A_large};
    arm_mat_mult_f32(&A_large_mat, &identity_mat, &large_result_mat);
    check_true(eps_close_matrix(A_large, large_result, 4, 4, FLT_EPSILON), "identity post-multiplication");
    arm_mat_mult_f32(&identity_mat, &A_large_mat, &large_result_mat);
    check_true(eps_close_matrix(A_large, large_result, 4, 4, FLT_EPSILON), "identity pre-multiplication");

    float B_large[4 * 4] = {5.24829, 6.21496, 3.27374, 3.49223, 1.52040, 3.70849, 7.21884, 0.41667,
                            7.77438, 8.24807, 8.63347, 2.01096, 8.29170, 1.46735, 8.53606, 5.14221};
    float C_large[4 * 4] = {1.92304, 0.14043, 1.64762, 4.97396, 0.68077, 4.99275, 7.04041, 2.44857,
                            8.22049, 1.66745, 7.94150, 0.56302, 2.68638, 7.59450, 1.43236, 3.59834};
    float large_multiplication_expected[4 * 4]
        = {50.6168, 63.7473, 83.4036, 55.732,  65.9102, 33.9305, 86.5396, 22.2066,
           96.939,  71.9404, 142.322, 70.9624, 100.929, 61.7765, 99.1469, 68.1449};
    arm_matrix_instance_f32 B_large_mat = {4, 4, B_large};
    arm_matrix_instance_f32 C_large_mat = {4, 4, C_large};
    arm_mat_mult_f32(&B_large_mat, &C_large_mat, &large_result_mat);
    // Expected values are only given to ~6 significant figures, so FLT_EPSILON (the old
    // tolerance) could never pass; 1e-3 absolute on values ~100 is the precision of the data.
    check_true(eps_close_matrix(large_result, large_multiplication_expected, 4, 4, 1e-3f),
               "4x4 product (6 significant figures)");

    // Documents a CMSIS gotcha that bit iterate.c: arm_mat_inverse_f32 destroys its input
    float S[4] = {4, 7, 2, 6};
    float S_copy[4] = {4, 7, 2, 6};
    float S_inv[4];
    arm_matrix_instance_f32 S_mat = {2, 2, S};
    arm_matrix_instance_f32 S_inv_mat = {2, 2, S_inv};
    arm_status st = arm_mat_inverse_f32(&S_mat, &S_inv_mat);
    float S_inv_expected[4] = {0.6f, -0.7f, -0.2f, 0.4f};
    check_true(st == ARM_MATH_SUCCESS, "2x2 inverse succeeds");
    check_vec_close("2x2 inverse value", S_inv, S_inv_expected, 4, 1e-6f);
    check_true(!eps_close_matrix(S, S_copy, 2, 2, 1e-6f),
               "arm_mat_inverse_f32 overwrites its source (so callers must pass a copy)");
    float sing[4] = {1, 2, 2, 4};
    arm_matrix_instance_f32 sing_mat = {2, 2, sing};
    check_true(arm_mat_inverse_f32(&sing_mat, &S_inv_mat) == ARM_MATH_SINGULAR, "singular matrix is reported");
}

// =============================================================================================
// Quaternions
// =============================================================================================
void test_quaternion(void) {
    begin_suite("quaternion");
    float res_quat[4];
    float res_vec[3];

    float id_quat[4] = {1, 0, 0, 0};
    quat_multiply(id_quat, id_quat, res_quat);
    check_vec_close("identity * identity", res_quat, id_quat, 4, FLT_EPSILON);

    float rp_quat_1[4] = {0.59513523, 0.53820488, 0.00657071, 0.5967465};
    float rp_quat_2[4] = {0.78805728, 0.07242287, 0.48362778, 0.37393157};
    quat_multiply(rp_quat_1, rp_quat_2, res_quat);
    float expected_res_rp[4] = {0.203702175001, 0.181091486099, 0.134968324367, 0.952625236179};
    check_vec_close("multiply random pair", res_quat, expected_res_rp, 4, 1e-6f);

    quat_multiply(rp_quat_2, rp_quat_1, res_quat);
    float expected_res_rp_inv[4] = {0.203702175001, 0.753383864322, 0.451035727503, 0.432995312932};
    check_vec_close("multiply reversed pair (non-commutative)", res_quat, expected_res_rp_inv, 4, 1e-6f);

    float in_place[4] = {rp_quat_1[0], rp_quat_1[1], rp_quat_1[2], rp_quat_1[3]};
    quat_multiply(in_place, rp_quat_2, in_place);
    check_vec_close("multiply output may alias left input", in_place, expected_res_rp, 4, 1e-6f);

    float unnorm_quat[4] = {0.44245429, 0.97492992, 0.24490942, 0.74845996};
    float expected_norm_quat[4] = {0.33290518, 0.73354294, 0.18427127, 0.56314562};
    quat_norm(unnorm_quat, res_quat);
    check_vec_close("normalization", res_quat, expected_norm_quat, 4, 1e-6f);
    float zero_q[4] = {0, 0, 0, 0};
    quat_norm(zero_q, res_quat);
    check_vec_close("normalizing the zero quaternion gives identity (not NaN)", res_quat, id_quat, 4, 0);

    float quat_to_invert[4] = {7, 4, 5, 9};
    float expected_inverted_quat[4] = {0.0409357, -0.0233918, -0.0292398, -0.0526316};
    quat_inv(quat_to_invert, res_quat);
    check_vec_close("inverse of non-unit quaternion", res_quat, expected_inverted_quat, 4, 1e-6f);
    float prod[4];
    quat_multiply(quat_to_invert, res_quat, prod);
    check_vec_close("q * q^-1 = identity", prod, id_quat, 4, 1e-6f);

    float conj[4];
    quat_conj(rp_quat_1, conj);
    quat_inv(rp_quat_1, res_quat);
    check_vec_close("conjugate == inverse for unit quaternion", conj, res_quat, 4, 1e-5f);

    float expected_rp_1_rotvec[3] = {1.2502, 0.0153, 1.3862};
    float rp1_copy[4];
    memcpy(rp1_copy, rp_quat_1, sizeof(rp1_copy));
    quat2rotationvec(rp_quat_1, res_vec);
    check_vec_close("quat -> rotation vector", res_vec, expected_rp_1_rotvec, 3, 1e-4f);

    float neg_rp1[4] = {-rp_quat_1[0], -rp_quat_1[1], -rp_quat_1[2], -rp_quat_1[3]};
    float neg_copy[4];
    memcpy(neg_copy, neg_rp1, sizeof(neg_copy));
    quat2rotationvec(neg_rp1, res_vec);
    check_vec_close("quat -> rotation vector picks the short way for -q", res_vec, expected_rp_1_rotvec, 3, 1e-4f);
    check_vec_close("quat -> rotation vector does not modify its input (old version negated it)", neg_rp1,
                    neg_copy, 4, 0);

    float expected_quat_diff[4] = {0.7343, -0.66718, 0.12461, 0.012084};
    quat_diff(rp_quat_1, rp_quat_2, res_quat);
    check_vec_close("quat_diff", res_quat, expected_quat_diff, 4, 1e-4f);
    float robustness_test[4];
    quat_multiply(res_quat, rp_quat_1, robustness_test);
    check_vec_close("quat_diff(a, b) * a == b", robustness_test, rp_quat_2, 4, 1e-4f);

    float random_vec[3] = {0.96667295, 0.7002543, 0.61082435};
    float expected_quat_conv[4] = {0.78355, 0.44793, 0.32448, 0.28304};
    rotationvec2quat(random_vec, res_quat);
    check_vec_close("rotation vector -> quat", res_quat, expected_quat_conv, 4, 1e-4f);

    float expected_rotated_vec[3] = {0.1828, 0.1028, 1.3244};
    quat_apply(rp_quat_1, random_vec, res_vec);
    check_vec_close("quat_apply", res_vec, expected_rotated_vec, 3, 1e-4f);

    // 90 deg about z maps x -> y (active rotation)
    float rz90[3] = {0, 0, (float)M_PI / 2};
    float qz90[4];
    rotationvec2quat(rz90, qz90);
    float ex[3] = {1, 0, 0}, ey[3] = {0, 1, 0};
    quat_apply(qz90, ex, res_vec);
    check_vec_close("90 deg about z maps x to y", res_vec, ey, 3, 1e-6f);
    float qz90_expected[4] = {(float)M_SQRT1_2, 0, 0, (float)M_SQRT1_2};
    check_vec_close("90 deg about z quaternion value", qz90, qz90_expected, 4, 1e-6f);

    // 180 deg rotation vector round trip
    float rx180[3] = {(float)M_PI, 0, 0};
    float qx180[4];
    rotationvec2quat(rx180, qx180);
    quat2rotationvec(qx180, res_vec);
    check_close("180 deg rotation vector round trip (angle)", l2_norm(res_vec, 3), M_PI, 1e-5);

    // Randomized properties
    rng_seed(101);
    float worst_roundtrip = 0, worst_apply = 0, worst_rotm = 0, worst_norm = 0, worst_rotm_quat = 0;
    for (int trial = 0; trial < 500; trial++) {
        float axis[3];
        random_unit_vec(axis);
        float angle = (float)M_PI * uniform_01();
        float rv[3] = {axis[0] * angle, axis[1] * angle, axis[2] * angle};
        float q[4], rv_back[3];
        rotationvec2quat(rv, q);
        quat2rotationvec(q, rv_back);
        for (int i = 0; i < 3; i++) {
            worst_roundtrip = fmaxf(worst_roundtrip, fabsf(rv_back[i] - rv[i]));
        }
        worst_norm = fmaxf(worst_norm, fabsf(l2_norm(q, 4) - 1.0f));

        float v[3], got[3], expected[3];
        random_unit_vec(v);
        v[0] *= 7.0f;
        v[1] *= 7.0f;
        v[2] *= 7.0f;
        quat_apply(q, v, got);
        ref_rotate(q, v, expected);
        for (int i = 0; i < 3; i++) {
            worst_apply = fmaxf(worst_apply, fabsf(got[i] - expected[i]));
        }

        float R[9], Rv[3];
        quat2rotm(q, R);
        for (int i = 0; i < 3; i++) {
            Rv[i] = R[i * 3] * v[0] + R[i * 3 + 1] * v[1] + R[i * 3 + 2] * v[2];
            worst_rotm = fmaxf(worst_rotm, fabsf(Rv[i] - expected[i]));
        }
        float q_back[4];
        rotm_to_quat(R, q_back);
        worst_rotm_quat = fmaxf(worst_rotm_quat, quat_angle_between(q, q_back));
    }
    check_less("rotation vector -> quat -> rotation vector round trip (500 random, max abs err)", worst_roundtrip, 1e-4);
    check_less("rotationvec2quat returns unit quaternions", worst_norm, 1e-6);
    check_less("quat_apply matches explicit rotation matrix formula", worst_apply, 1e-4);
    check_less("quat2rotm matches explicit rotation matrix formula", worst_rotm, 1e-4);
    check_less("rotm_to_quat(quat2rotm(q)) == q (all four CMSIS branches, rad)", worst_rotm_quat, 1e-3);

    // Small angles: the old acos-based version returned exactly 0 below ~1e-3 rad in float32
    float tiny[3] = {2e-5f, -1e-5f, 3e-5f};
    float q_tiny[4], tiny_back[3];
    rotationvec2quat(tiny, q_tiny);
    quat2rotationvec(q_tiny, tiny_back);
    check_vec_close("tiny rotation (3.7e-5 rad) survives the round trip", tiny_back, tiny, 3, 1e-7f);
    float small[3] = {5e-4f, 0, 0};
    float q_small[4], small_back[3];
    rotationvec2quat(small, q_small);
    quat2rotationvec(q_small, small_back);
    check_vec_close("small rotation (5e-4 rad) survives the round trip", small_back, small, 3, 1e-7f);
    float sub_eps[3] = {1e-8f, 0, 0};
    rotationvec2quat(sub_eps, res_quat);
    check_close("rotation below 1e-6 rad still gives a unit quaternion", l2_norm(res_quat, 4), 1, 1e-6);

    // Composition: rotating by a then b equals rotating by b*a
    float qa[4], qb[4], qba[4], v[3] = {0.3f, -1.2f, 2.0f}, t1[3], t2[3], t3[3];
    random_quat(qa);
    random_quat(qb);
    quat_multiply(qb, qa, qba);
    quat_apply(qa, v, t1);
    quat_apply(qb, t1, t2);
    quat_apply(qba, v, t3);
    check_vec_close("apply(b, apply(a, v)) == apply(b*a, v)", t2, t3, 3, 1e-5f);

    float len_before = l2_norm(v, 3);
    quat_apply(qa, v, t1);
    check_close("rotation preserves vector length", l2_norm(t1, 3), len_before, 1e-5);
    float apply_alias[3] = {v[0], v[1], v[2]};
    quat_apply(qa, apply_alias, apply_alias);
    check_vec_close("quat_apply output may alias input", apply_alias, t1, 3, 1e-6f);
}

// =============================================================================================
// QUEST
// =============================================================================================
static float quest_error_deg(const float* q_true_b2r, const float* ref, int n, float noise_rad, bool normalize) {
    float body[3 * 8];
    float q_r2b[4];
    quat_conj(q_true_b2r, q_r2b);
    for (int i = 0; i < n; i++) {
        quat_apply(q_r2b, &ref[3 * i], &body[3 * i]);
        if (noise_rad > 0) {
            float axis[3], dq[4];
            random_unit_vec(axis);
            float ang = normal_sample(0, noise_rad);
            float rv[3] = {axis[0] * ang, axis[1] * ang, axis[2] * ang};
            rotationvec2quat(rv, dq);
            quat_apply(dq, &body[3 * i], &body[3 * i]);
        }
        if (!normalize) {
            // Scale body vectors arbitrarily: QUEST must not care
            float s = 0.5f + 3.0f * uniform_01();
            body[3 * i] *= s;
            body[3 * i + 1] *= s;
            body[3 * i + 2] *= s;
        }
    }
    float guess[4];
    if (!quest(body, ref, n, guess)) {
        return 999.0f;
    }
    return RAD2DEG(quat_angle_between(guess, q_true_b2r));
}

void test_quest(void) {
    begin_suite("QUEST");

    // Original test case
    float true_ref_to_body[4] = {0.78355, 0.44793, 0.32448, 0.28304};
    float true_body_to_ref[4];
    quat_inv(true_ref_to_body, true_body_to_ref);
    float ref_1[3] = {1, 2, 3};
    float ref_2[3] = {-4, 3, -6};
    float ref[6] = {1, 2, 3, -4, 3, -6};
    float body_1[3], body_2[3];
    quat_apply(true_ref_to_body, ref_1, body_1);
    quat_apply(true_ref_to_body, ref_2, body_2);
    float body[6] = {body_1[0], body_1[1], body_1[2], body_2[0], body_2[1], body_2[2]};
    float guess[4];
    bool ok = quest(body, ref, 2, guess);
    check_true(ok, "returns true for two good vectors");
    check_vec_close("original case (unnormalized inputs)", guess, true_body_to_ref, 4, 1e-4f);
    check_close("result is a unit quaternion", l2_norm(guess, 4), 1, 1e-6);
    check_true(guess[0] >= 0, "result has w >= 0");

    float body_copy[6], ref_copy[6];
    memcpy(body_copy, body, sizeof(body));
    memcpy(ref_copy, ref, sizeof(ref));
    quest(body, ref, 2, guess);
    check_true(memcmp(body, body_copy, sizeof(body)) == 0 && memcmp(ref, ref_copy, sizeof(ref)) == 0,
               "inputs are not modified");

    rng_seed(202);
    float worst = 0, worst_unnorm = 0, worst3 = 0;
    for (int trial = 0; trial < 300; trial++) {
        float q[4], r[9];
        random_quat(q);
        random_unit_vec(&r[0]);
        do {
            random_unit_vec(&r[3]);
        } while (vec_angle(&r[0], &r[3]) < DEG2RAD(15) || vec_angle(&r[0], &r[3]) > DEG2RAD(165));
        random_unit_vec(&r[6]);
        worst = fmaxf(worst, quest_error_deg(q, r, 2, 0, true));
        worst_unnorm = fmaxf(worst_unnorm, quest_error_deg(q, r, 2, 0, false));
        worst3 = fmaxf(worst3, quest_error_deg(q, r, 3, 0, true));
    }
    check_less("300 random attitudes, 2 exact vectors: max error (deg)", worst, 0.05);
    check_less("300 random attitudes, 2 exact vectors of arbitrary length: max error (deg)", worst_unnorm, 0.05);
    check_less("300 random attitudes, 3 exact vectors: max error (deg)", worst3, 0.05);

    // Near/at 180 degrees: the plain QUEST closed form is singular here
    float r2[6] = {0.6f, 0.8f, 0.0f, 0.0f, 0.28f, 0.96f};
    const char* names[4] = {"180 deg about x", "180 deg about y", "180 deg about z", "180 deg about random axis"};
    float axes[4][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.48f, -0.6f, 0.64f}};
    for (int k = 0; k < 4; k++) {
        float q180[4] = {0, axes[k][0], axes[k][1], axes[k][2]};
        quat_norm(q180, q180);
        char name[96];
        snprintf(name, sizeof(name), "%s: error (deg)", names[k]);
        check_less(name, quest_error_deg(q180, r2, 2, 0, true), 0.1);
    }
    float rv179[3] = {0, DEG2RAD(179.0f), 0};
    float q179[4];
    rotationvec2quat(rv179, q179);
    check_less("179 deg about y: error (deg)", quest_error_deg(q179, r2, 2, 0, true), 0.1);

    // Noisy measurements: 1 deg per vector
    rng_seed(203);
    float sum_err = 0;
    int n_trials = 200;
    for (int trial = 0; trial < n_trials; trial++) {
        float q[4], r[6];
        random_quat(q);
        random_unit_vec(&r[0]);
        do {
            random_unit_vec(&r[3]);
        } while (vec_angle(&r[0], &r[3]) < DEG2RAD(45) || vec_angle(&r[0], &r[3]) > DEG2RAD(135));
        sum_err += quest_error_deg(q, r, 2, DEG2RAD(1.0f), true);
    }
    check_less("1 deg measurement noise: mean error (deg)", sum_err / n_trials, 2.0);

    // Weighting: a garbage third vector with a tiny weight barely matters
    float q_w[4] = {0.9f, 0.1f, -0.3f, 0.2f};
    quat_norm(q_w, q_w);
    float rw[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    float bw[9], q_w_r2b[4];
    quat_conj(q_w, q_w_r2b);
    for (int i = 0; i < 3; i++) {
        quat_apply(q_w_r2b, &rw[3 * i], &bw[3 * i]);
    }
    bw[6] = 1;
    bw[7] = 0;
    bw[8] = 0; // wrong
    float weights[3] = {1.0f, 1.0f, 1e-4f};
    quest_weighted(bw, rw, weights, 3, guess);
    check_less("tiny weight on a bad vector: error (deg)", RAD2DEG(quat_angle_between(guess, q_w)), 0.1);

    // Degenerate input
    float par_ref[6] = {1, 0, 0, 2, 0, 0};
    float par_body[6] = {0, 1, 0, 0, 3, 0};
    ok = quest(par_body, par_ref, 2, guess);
    float id[4] = {1, 0, 0, 0};
    check_true(!ok, "parallel reference vectors are rejected");
    check_vec_close("parallel reference vectors give identity", guess, id, 4, 0);
    float zero_ref[6] = {0, 0, 0, 0, 1, 0};
    ok = quest(par_body, zero_ref, 2, guess);
    check_true(!ok, "zero-length vector is rejected");
    float nan_body[6] = {NAN, 0, 0, 0, 1, 0};
    ok = quest(nan_body, rw, 2, guess);
    check_true(!ok && all_finite(guess, 4), "NaN input is rejected without producing NaN");
}

// =============================================================================================
// UKF
// =============================================================================================
void test_ukf_internals(void) {
    begin_suite("UKF internals");

    float lambda = calculate_lambda(STATE_SIZE, UKF_ALPHA, UKF_KAPPA);
    check_close("lambda = alpha^2 (n + kappa) - n", lambda,
                UKF_ALPHA * UKF_ALPHA * (STATE_SIZE + UKF_KAPPA) - STATE_SIZE, 1e-6);

    float wc[NUM_SIGMAS], wm[NUM_SIGMAS];
    get_weights(lambda, STATE_SIZE, UKF_ALPHA, UKF_BETA, wc, wm);
    float sum_m = 0, sum_c = 0;
    for (int i = 0; i < NUM_SIGMAS; i++) {
        sum_m += wm[i];
        sum_c += wc[i];
    }
    check_close("mean weights sum to 1", sum_m, 1, 1e-5);
    check_close("covariance weights sum to 1 + (1 - alpha^2 + beta)", sum_c, 1 + (1 - UKF_ALPHA * UKF_ALPHA + UKF_BETA), 1e-5);

    // Old parameters (alpha = 0.01, kappa = 3 - n): the centre weight is huge and negative,
    // which is what made the float32 filter numerically fragile
    float lam_old = calculate_lambda(STATE_SIZE, 0.01f, 3.0f - STATE_SIZE);
    float wc_old[NUM_SIGMAS], wm_old[NUM_SIGMAS];
    get_weights(lam_old, STATE_SIZE, 0.01f, 2.0f, wc_old, wm_old);
    check_less("old alpha=0.01 parameters had a centre weight below -1e4 (documented reason for the change)", wm_old[0], -1e4);

    // Sigma points must reproduce the mean and a NON-DIAGONAL covariance. The old code used rows
    // of the lower Cholesky factor (plus an uninitialized upper triangle), which fails this.
    float A[36];
    rng_seed(303);
    for (int i = 0; i < 36; i++) {
        A[i] = uniform_pm1() * 0.1f;
    }
    float P[36] = {0};
    for (int r = 0; r < 6; r++) {
        for (int c = 0; c < 6; c++) {
            float s = 0;
            for (int k = 0; k < 6; k++) {
                s += A[r * 6 + k] * A[c * 6 + k];
            }
            P[r * 6 + c] = s + (r == c ? 0.01f : 0.0f);
        }
    }
    float x[6] = {0.01f, -0.02f, 0.03f, 0.001f, -0.002f, 0.0005f};
    float sigmas[NUM_SIGMAS * STATE_SIZE];
    // Poison the output so any unwritten element shows up
    for (int i = 0; i < NUM_SIGMAS * STATE_SIZE; i++) {
        sigmas[i] = NAN;
    }
    bool ok = get_sigma_points(lambda, x, P, sigmas);
    check_true(ok, "sigma points generated for a valid non-diagonal covariance");
    check_true(all_finite(sigmas, NUM_SIGMAS * STATE_SIZE), "every sigma point is written and finite");
    check_vec_close("first sigma point is the mean", sigmas, x, 6, 0);

    float mean[6] = {0}, cov[36] = {0};
    for (int i = 0; i < NUM_SIGMAS; i++) {
        for (int j = 0; j < 6; j++) {
            mean[j] += wm[i] * sigmas[i * 6 + j];
        }
    }
    for (int i = 0; i < NUM_SIGMAS; i++) {
        for (int r = 0; r < 6; r++) {
            for (int c = 0; c < 6; c++) {
                cov[r * 6 + c] += wc[i] * (sigmas[i * 6 + r] - x[r]) * (sigmas[i * 6 + c] - x[c]);
            }
        }
    }
    check_vec_close("weighted sigma mean reproduces the state", mean, x, 6, 1e-6f);
    check_vec_close("weighted sigma covariance reproduces a non-diagonal P", cov, P, 36, 1e-5f);

    float P_saved[36];
    memcpy(P_saved, P, sizeof(P));
    get_sigma_points(lambda, x, P, sigmas);
    check_vec_close("get_sigma_points does not modify P", P, P_saved, 36, 0);

    // ensure_psd
    float healthy[36];
    memcpy(healthy, P, sizeof(P));
    check_true(ensure_psd(healthy, 6), "ensure_psd accepts a healthy matrix");
    check_vec_close("ensure_psd leaves a healthy matrix unchanged (old version always added 1e-6)", healthy, P, 36, 0);

    float asym[9] = {2, 1, 0, 0.5f, 2, 0, 0, 0, 1};
    check_true(ensure_psd(asym, 3), "ensure_psd repairs an asymmetric matrix");
    check_true(matrix_is_symmetric(asym, 3, 0), "ensure_psd output is exactly symmetric");
    check_close("ensure_psd symmetrizes by averaging", asym[1], 0.75, 1e-7);

    float indefinite[9] = {1, 2, 0, 2, 1, 0, 0, 0, 1}; // eigenvalues 3, -1, 1
    check_true(ensure_psd(indefinite, 3) && matrix_is_pd(indefinite, 3), "ensure_psd makes an indefinite matrix PD");

    float neg_diag[9] = {-1e-3f, 0, 0, 0, 1, 0, 0, 0, 1};
    check_true(ensure_psd(neg_diag, 3) && matrix_is_pd(neg_diag, 3), "ensure_psd fixes a negative diagonal");

    float nan_m[9] = {1, 0, 0, 0, NAN, 0, 0, 0, 1};
    check_true(!ensure_psd(nan_m, 3), "ensure_psd rejects NaN");
    float too_big[49] = {0};
    check_true(!ensure_psd(too_big, 7), "ensure_psd rejects sizes above STATE_SIZE");
}

// Shared single-step truth simulation for the UKF behaviour tests
typedef struct {
    float q_true[4]; // body -> ref
    float omega[3];  // rad/s, body
    float bias[3];   // rad/s
    float gyro_noise;
    float vec_noise; // rad
} ukf_truth_t;

static void truth_step(ukf_truth_t* t, float dt, float* gyro_out) {
    for (int i = 0; i < 3; i++) {
        gyro_out[i] = t->omega[i] + t->bias[i] + normal_sample(0, t->gyro_noise);
    }
    float rv[3] = {t->omega[0] * dt, t->omega[1] * dt, t->omega[2] * dt};
    float dq[4];
    rotationvec2quat(rv, dq);
    quat_multiply(t->q_true, dq, t->q_true);
    quat_norm(t->q_true, t->q_true);
}

static void truth_measure(const ukf_truth_t* t, const float* ref_vec, float* body_vec) {
    ref_rotate_inv(t->q_true, ref_vec, body_vec);
    if (t->vec_noise > 0) {
        for (int i = 0; i < 3; i++) {
            body_vec[i] += normal_sample(0, t->vec_noise);
        }
    }
    normalize_vec(body_vec, 3);
}

void test_ukf_behaviour(void) {
    begin_suite("UKF behaviour");

    float mag_ref[3] = {0.3f, -0.5f, 0.81f};
    float sun_ref[3] = {0.9f, 0.4f, -0.1f};
    normalize_vec(mag_ref, 3);
    normalize_vec(sun_ref, 3);
    float ref2[6] = {mag_ref[0], mag_ref[1], mag_ref[2], sun_ref[0], sun_ref[1], sun_ref[2]};

    float Q[36], R2[36], R1[9];
    float q_diag[6] = {1e-6f, 1e-6f, 1e-6f, 1e-10f, 1e-10f, 1e-10f};
    diag(Q, 6, q_diag);
    float r2_diag[6] = {1e-4f, 1e-4f, 1e-4f, 1e-4f, 1e-4f, 1e-4f};
    diag(R2, 6, r2_diag);
    float r1_diag[3] = {1e-4f, 1e-4f, 1e-4f};
    diag(R1, 3, r1_diag);

    // --- Predict-only step integrates the gyro exactly ---
    {
        float q_true[4] = {0.5f, 0.5f, -0.5f, 0.5f};
        float q_est[4];
        memcpy(q_est, q_true, sizeof(q_est));
        float x[6] = {0};
        float P[36];
        float p_diag[6] = {1e-4f, 1e-4f, 1e-4f, 1e-6f, 1e-6f, 1e-6f};
        diag(P, 6, p_diag);
        float omega[3] = {0.02f, -0.01f, 0.03f};
        float dt = 0.1f;
        bool all_ok = true, grows = true;
        float last_var = P[0];
        for (int k = 0; k < 200; k++) {
            float rv[3] = {omega[0] * dt, omega[1] * dt, omega[2] * dt};
            float dq[4];
            rotationvec2quat(rv, dq);
            quat_multiply(q_true, dq, q_true);
            all_ok &= iterate(x, q_est, P, NULL, NULL, 0, omega, Q, NULL, dt, x, q_est, P) == UKF_OK;
            grows &= P[0] >= last_var;
            last_var = P[0];
        }
        check_true(all_ok, "predict-only (num_vecs = 0) steps succeed with NULL body/ref/R");
        check_less("predict-only tracks gyro integration over 200 steps (deg)", RAD2DEG(quat_angle_between(q_est, q_true)), 0.01);
        check_true(grows, "predict-only attitude variance never shrinks");
        // Growth = N Q_att + (N dt)^2 var(bias): the uncertain bias integrates into attitude error
        check_close("predict-only attitude variance = P0 + N Q + (N dt)^2 var(bias)", P[0],
                    1e-4 + 200 * 1e-6 + (200 * 0.1) * (200 * 0.1) * 1e-6, 3e-5);
    }

    // --- Perfect measurements starting at the truth stay at the truth ---
    {
        rng_seed(401);
        ukf_truth_t t = {{0.2812f, -0.6998f, 0.5497f, -0.3592f}, {DEG2RAD(2), DEG2RAD(0.5f), DEG2RAD(-1)}, {0, 0, 0}, 0, 0};
        quat_norm(t.q_true, t.q_true);
        float q_est[4], x[6] = {0}, P[36];
        memcpy(q_est, t.q_true, sizeof(q_est));
        float p_diag[6] = {1e-3f, 1e-3f, 1e-3f, 1e-6f, 1e-6f, 1e-6f};
        diag(P, 6, p_diag);
        float worst = 0;
        for (int k = 0; k < 300; k++) {
            float gyro[3], body[6];
            truth_step(&t, 0.5f, gyro);
            truth_measure(&t, &ref2[0], &body[0]);
            truth_measure(&t, &ref2[3], &body[3]);
            iterate(x, q_est, P, body, ref2, 2, gyro, Q, R2, 0.5f, x, q_est, P);
            worst = fmaxf(worst, quat_angle_between(q_est, t.q_true));
        }
        check_less("noise-free 2-vector run starting at truth: max error over 300 steps (deg)", RAD2DEG(worst), 0.05);
    }

    // --- 2-vector convergence with gyro bias and noise (the main in-sun mode) ---
    float bias_true[3] = {0.002f, -0.001f, 0.0015f};
    float q_after[4], x_after[6], P_after[36];
    ukf_truth_t t_after;
    {
        rng_seed(402);
        ukf_truth_t t = {{0.2812f, -0.6998f, 0.5497f, -0.3592f},
                         {DEG2RAD(2), DEG2RAD(0.5f), DEG2RAD(-1)},
                         {bias_true[0], bias_true[1], bias_true[2]},
                         0.001f,
                         0.01f};
        quat_norm(t.q_true, t.q_true);
        float err_rv[3] = {DEG2RAD(12), DEG2RAD(-10), DEG2RAD(8)}; // ~18 deg initial error
        float dq[4], q_est[4];
        rotationvec2quat(err_rv, dq);
        quat_multiply(t.q_true, dq, q_est);
        float x[6] = {0}, P[36];
        float p_diag[6] = {0.1f, 0.1f, 0.1f, 1e-4f, 1e-4f, 1e-4f};
        diag(P, 6, p_diag);
        float dt = 0.5f;
        bool all_ok = true, sym = true, pd = true;
        float tail_err = 0;
        int tail_n = 0;
        float err_at_50 = 0;
        for (int k = 0; k < 2000; k++) {
            float gyro[3], body[6];
            truth_step(&t, dt, gyro);
            truth_measure(&t, &ref2[0], &body[0]);
            truth_measure(&t, &ref2[3], &body[3]);
            all_ok &= iterate(x, q_est, P, body, ref2, 2, gyro, Q, R2, dt, x, q_est, P) == UKF_OK;
            sym &= matrix_is_symmetric(P, 6, 0);
            pd &= matrix_is_pd(P, 6);
            float e = RAD2DEG(quat_angle_between(q_est, t.q_true));
            if (k == 50) {
                err_at_50 = e;
            }
            if (k >= 1800) {
                tail_err += e;
                tail_n++;
            }
        }
        float bias_err[3] = {x[3] - bias_true[0], x[4] - bias_true[1], x[5] - bias_true[2]};
        check_true(all_ok, "2-vector: every step returns UKF_OK");
        check_true(sym, "2-vector: covariance stays exactly symmetric");
        check_true(pd, "2-vector: covariance stays positive definite");
        check_less("2-vector: converged from ~18 deg within 50 steps (deg)", err_at_50, 2.0);
        check_less("2-vector: mean error over last 200 of 2000 steps (deg)", tail_err / tail_n, 0.5);
        check_less("2-vector: gyro bias error after 2000 steps (deg/s)", RAD2DEG(l2_norm(bias_err, 3)), 0.02);
        check_true(exactly(x[0], 0) && exactly(x[1], 0) && exactly(x[2], 0), "attitude part of the returned error state is zeroed");
        memcpy(q_after, q_est, sizeof(q_after));
        memcpy(x_after, x, sizeof(x_after));
        memcpy(P_after, P, sizeof(P_after));
        t_after = t;
    }

    // Orbit-like field: the magnetic field direction turns ~two times per 90 min orbit
    float field_axis[3] = {0.3f, 0.2f, 0.93f};
    normalize_vec(field_axis, 3);
    const float field_rate = 2.0f * 2.0f * (float)M_PI / 5400.0f;

    // --- Eclipse: magnetometer only (num_vecs = 1) keeps the converged estimate bounded ---
    for (int rotating = 0; rotating <= 1; rotating++) {
        rng_seed(403);
        ukf_truth_t t = t_after;
        float q_est[4], x[6], P[36], ref_now[3];
        memcpy(q_est, q_after, sizeof(q_est));
        memcpy(x, x_after, sizeof(x));
        memcpy(P, P_after, sizeof(P));
        memcpy(ref_now, mag_ref, sizeof(ref_now));
        bool all_ok = true;
        float worst = 0;
        for (int k = 0; k < 2000; k++) { // ~17 minutes of eclipse
            float gyro[3], body[3];
            if (rotating) {
                float rv[3] = {field_axis[0] * field_rate * 0.5f, field_axis[1] * field_rate * 0.5f,
                               field_axis[2] * field_rate * 0.5f};
                float field_step[4];
                rotationvec2quat(rv, field_step);
                quat_apply(field_step, ref_now, ref_now);
            }
            truth_step(&t, 0.5f, gyro);
            truth_measure(&t, ref_now, body);
            all_ok &= iterate(x, q_est, P, body, ref_now, 1, gyro, Q, R1, 0.5f, x, q_est, P) == UKF_OK;
            worst = fmaxf(worst, RAD2DEG(quat_angle_between(q_est, t.q_true)));
        }
        if (rotating) {
            check_true(all_ok, "mag-only, orbit-like field: every step returns UKF_OK");
            // Rotation about the current field direction is only observed as the field turns,
            // so ~1-2.5 deg excursions are expected with 0.6 deg (1 sigma) magnetometer noise
            check_less("mag-only, orbit-like field: max error during 1000 s eclipse after convergence (deg)", worst, 3.0);
        } else {
            // Rotation about a fixed inertial vector is unobservable, so this only drifts with
            // the residual bias error; it must stay bounded, not stay accurate
            check_true(all_ok, "mag-only, fixed field (worst case): every step returns UKF_OK");
            check_less("mag-only, fixed field (worst case): max error during 1000 s eclipse (deg)", worst, 5.0);
        }
    }

    // --- Mag-only from a cold start converges when the reference direction changes ---
    // (Rotation about a single FIXED inertial vector is never observable, whatever the body does;
    // in orbit the field direction turns ~720 deg per orbit, which is what makes this work.)
    {
        rng_seed(404);
        ukf_truth_t t = {{1, 0, 0, 0}, {DEG2RAD(3), DEG2RAD(-2), DEG2RAD(1)}, {0, 0, 0}, 0.0005f, 0.01f};
        float err_rv[3] = {DEG2RAD(10), DEG2RAD(5), DEG2RAD(-5)};
        float q_est[4], dq[4];
        rotationvec2quat(err_rv, dq);
        quat_multiply(t.q_true, dq, q_est);
        float x[6] = {0}, P[36];
        float p_diag[6] = {0.05f, 0.05f, 0.05f, 1e-8f, 1e-8f, 1e-8f};
        diag(P, 6, p_diag);
        float tail = 0;
        float ref_now[3];
        memcpy(ref_now, mag_ref, sizeof(ref_now));
        for (int k = 0; k < 1500; k++) {
            float gyro[3], body[3];
            float rv[3] = {field_axis[0] * field_rate * 0.5f, field_axis[1] * field_rate * 0.5f,
                           field_axis[2] * field_rate * 0.5f};
            float field_step[4];
            rotationvec2quat(rv, field_step);
            quat_apply(field_step, ref_now, ref_now);
            truth_step(&t, 0.5f, gyro);
            truth_measure(&t, ref_now, body);
            iterate(x, q_est, P, body, ref_now, 1, gyro, Q, R1, 0.5f, x, q_est, P);
            if (k >= 1300) {
                tail += RAD2DEG(quat_angle_between(q_est, t.q_true));
            }
        }
        check_less("mag-only cold start, orbit-like field rotation: mean error over last 200 steps (deg)", tail / 200, 2.0);
    }

    // --- API robustness ---
    {
        float q[4] = {1, 0, 0, 0}, x[6] = {0}, P[36], body[6], gyro[3] = {0.01f, 0, 0};
        float p_diag[6] = {1e-3f, 1e-3f, 1e-3f, 1e-6f, 1e-6f, 1e-6f};
        diag(P, 6, p_diag);
        memcpy(body, ref2, sizeof(body));
        float q_in[4], x_in[6], P_in[36], body_in[6], ref_in[6], gyro_in[3], Q_in[36], R_in[36];
        memcpy(q_in, q, sizeof(q));
        memcpy(x_in, x, sizeof(x));
        memcpy(P_in, P, sizeof(P));
        memcpy(body_in, body, sizeof(body));
        memcpy(ref_in, ref2, sizeof(ref2));
        memcpy(gyro_in, gyro, sizeof(gyro));
        memcpy(Q_in, Q, sizeof(Q));
        memcpy(R_in, R2, sizeof(R2));

        float q_out[4], x_out[6], P_out[36];
        ukf_status_t st = iterate(x, q, P, body, ref2, 2, gyro, Q, R2, 0.1f, x_out, q_out, P_out);
        check_true(st == UKF_OK, "valid step succeeds");
        check_true(memcmp(q, q_in, sizeof(q)) == 0 && memcmp(x, x_in, sizeof(x)) == 0 && memcmp(P, P_in, sizeof(P)) == 0
                       && memcmp(body, body_in, sizeof(body)) == 0 && memcmp(ref2, ref_in, sizeof(ref_in)) == 0
                       && memcmp(gyro, gyro_in, sizeof(gyro)) == 0 && memcmp(Q, Q_in, sizeof(Q)) == 0
                       && memcmp(R2, R_in, sizeof(R_in)) == 0,
                   "inputs are not modified");

        float q_alias[4], x_alias[6], P_alias[36];
        memcpy(q_alias, q, sizeof(q));
        memcpy(x_alias, x, sizeof(x));
        memcpy(P_alias, P, sizeof(P));
        iterate(x_alias, q_alias, P_alias, body, ref2, 2, gyro, Q, R2, 0.1f, x_alias, q_alias, P_alias);
        check_true(memcmp(q_alias, q_out, sizeof(q_out)) == 0 && memcmp(P_alias, P_out, sizeof(P_out)) == 0,
                   "outputs may alias inputs (same result as separate buffers)");

        float nan_gyro[3] = {NAN, 0, 0};
        st = iterate(x, q, P, body, ref2, 2, nan_gyro, Q, R2, 0.1f, x_out, q_out, P_out);
        check_true(st == UKF_ERR_ARGS, "NaN gyro is rejected with UKF_ERR_ARGS");
        check_true(memcmp(q_out, q, sizeof(q)) == 0 && memcmp(P_out, P, sizeof(P)) == 0,
                   "on failure the outputs are copies of the inputs");
        float nan_body[6] = {NAN, 0, 0, 0, 1, 0};
        check_true(iterate(x, q, P, nan_body, ref2, 2, gyro, Q, R2, 0.1f, x_out, q_out, P_out) == UKF_ERR_ARGS,
                   "NaN measurement is rejected");
        check_true(iterate(x, q, P, body, ref2, 3, gyro, Q, R2, 0.1f, x_out, q_out, P_out) == UKF_ERR_ARGS,
                   "num_vecs > MAX_MSMT_VECS is rejected");
        check_true(iterate(x, q, P, body, ref2, 2, gyro, Q, R2, -0.1f, x_out, q_out, P_out) == UKF_ERR_ARGS,
                   "negative dt is rejected");
        float P_nan[36];
        memcpy(P_nan, P, sizeof(P));
        P_nan[7] = NAN;
        check_true(iterate(x, q, P_nan, body, ref2, 2, gyro, Q, R2, 0.1f, x_out, q_out, P_out) != UKF_OK,
                   "NaN covariance is rejected");
        float R_zero[36] = {0};
        float P_zero[36] = {0};
        float Q_zero[36] = {0};
        ukf_status_t st_zero = iterate(x, q, P_zero, body, ref2, 2, gyro, Q_zero, R_zero, 0.1f, x_out, q_out, P_out);
        check_true(all_finite(q_out, 4) && all_finite(P_out, 36),
                   "all-zero P, Q and R never produce NaN (either repaired or rejected)");
        (void)st_zero;
    }
}

// Legacy simulation from the original test file: a single slowly-varying reference vector plus
// its finite-difference derivative as the second "vector" (the MUKF_Dot / simulate_1vec.m idea).
void test_iteration_1vec(void) {
    begin_suite("UKF legacy: 1 vector + derivative (simulate_1vec.m)");
    rng_seed(501);
    float dt = 0.5;
    float true_body_to_ref[4] = {0.2812, -0.6998, 0.5497, -0.3592};
    float true_ref_to_body[4];
    quat_norm(true_body_to_ref, true_body_to_ref);
    quat_inv(true_body_to_ref, true_ref_to_body);

    float current_guess[4];
    for (int i = 0; i < 4; i++) {
        current_guess[i] = true_body_to_ref[i] + normal_sample(0.05f, 0.1f);
    }
    quat_norm(current_guess, current_guess);
    float initial_err = RAD2DEG(quat_angle_between(current_guess, true_body_to_ref));

    float true_omega[3] = {2, 0.5, -0.5};
    float true_bias[3] = {0.001f, 0.001f, 0.001f};
    arm_scale_f32(true_omega, M_PI / 180.0f, true_omega, 3);

    float gyro_noise = 0.001;
    float msmt_noise = 0.03;

    float state[6] = {0, 0, 0, 0, 0, 0};
    float cov[36];
    float Q[36];
    float R[36];
    eye_f32(cov, 6);
    arm_scale_f32(cov, 0.1, cov, 36);
    memset(R, 0, sizeof(float) * 36);
    R[0] = 0.03f;
    R[7] = 0.03f;
    R[14] = 0.03f;
    R[21] = 1.0f;
    R[28] = 1.0f;
    R[35] = 1.0f;
    // The original Q = 0.01 I (a 0.1 rad attitude random walk every 0.5 s step) with the bias
    // ignored gives ~23 deg here even with a correct filter. Realistic process noise plus bias
    // estimation gets ~2 deg with the same R.
    memset(Q, 0, sizeof(Q));
    Q[0] = Q[7] = Q[14] = 1e-5f;
    Q[21] = Q[28] = Q[35] = 1e-10f;
    eye_f32(cov, 6);
    arm_scale_f32(cov, 0.1, cov, 36);
    cov[21] = cov[28] = cov[35] = 1e-4f;

    float ref[3] = {40, -40, 30};
    float last_ref[3] = {40, -40, 30};
    float last_body[3];
    quat_apply(true_ref_to_body, last_ref, last_body);

    bool all_ok = true;
    float tail_err = 0;
    int tail_n = 0;
    for (int iter = 0; iter < 1000; iter++) {
        float simulated_gyro_measurement[3];
        for (int i = 0; i < 3; i++) {
            simulated_gyro_measurement[i] = true_omega[i] + normal_sample(true_bias[i], gyro_noise);
        }

        float delta_vec[3];
        memcpy(delta_vec, true_omega, sizeof(float) * 3);
        arm_scale_f32(delta_vec, dt, delta_vec, 3);
        float delta_q[4];
        rotationvec2quat(delta_vec, delta_q);
        quat_multiply(true_body_to_ref, delta_q, true_body_to_ref);
        quat_norm(true_body_to_ref, true_body_to_ref);
        float new_ref_to_body[4];
        quat_inv(true_body_to_ref, new_ref_to_body);

        float d_ref[3];
        for (int i = 0; i < 3; i++) {
            d_ref[i] = (ref[i] - last_ref[i]) / dt;
        }
        float body[3];
        quat_apply(new_ref_to_body, ref, body);
        for (int i = 0; i < 3; i++) {
            body[i] = body[i] + normal_sample(0.0f, msmt_noise);
        }
        float d_body[3];
        float c_prod[3];
        cross(simulated_gyro_measurement, body, c_prod);
        for (int i = 0; i < 3; i++) {
            d_body[i] = (body[i] - last_body[i]) / dt + c_prod[i];
        }
        float body_full[6] = {body[0], body[1], body[2], d_body[0], d_body[1], d_body[2]};
        float reference[6] = {ref[0], ref[1], ref[2], d_ref[0], d_ref[1], d_ref[2]};

        float new_error_state[6];
        float new_quat[4];
        float new_P[36];
        all_ok &= iterate(state, current_guess, cov, body_full, reference, 2, simulated_gyro_measurement, Q, R, dt,
                          new_error_state, new_quat, new_P)
                  == UKF_OK;

        memcpy(state, new_error_state, sizeof(float) * 6);
        memcpy(cov, new_P, sizeof(float) * 36);
        memcpy(current_guess, new_quat, sizeof(float) * 4);

        float err = RAD2DEG(quat_angle_between(true_body_to_ref, current_guess));
        if (iter % 200 == 0) {
            printf("  iter %4d: attitude error %.3f deg\n", iter, err);
        }
        if (iter >= 800) {
            tail_err += err;
            tail_n++;
        }
        memcpy(last_body, body, sizeof(float) * 3);
        memcpy(last_ref, ref, sizeof(float) * 3);

        // THE IMPORTANT PART: REF NEEDS TO CHANGE A DECENT AMT
        for (int i = 0; i < 3; i++) {
            ref[i] = ref[i] + normal_sample(0.0f, 0.1f);
        }
    }
    printf("  initial error %.2f deg\n", initial_err);
    check_true(all_ok, "every step returns UKF_OK");
    // Before the fixes this diverged to >100 deg
    check_less("mean error over the last 200 of 1000 steps (deg)", tail_err / tail_n, 3.0);
}

// Legacy simulation from the original test file: two fixed, unnormalized reference vectors
void test_iteration_2vec(void) {
    begin_suite("UKF legacy: 2 fixed vectors (simulate_2vec.m)");
    rng_seed(502);
    float dt = 0.5;
    float true_body_to_ref[4] = {0.2812, -0.6998, 0.5497, -0.3592};
    quat_norm(true_body_to_ref, true_body_to_ref);

    float current_guess[4];
    for (int i = 0; i < 4; i++) {
        current_guess[i] = true_body_to_ref[i] + normal_sample(0.05f, 0.1f);
    }
    quat_norm(current_guess, current_guess);

    float true_omega[3] = {2, 0.5, -1};
    float true_bias[3] = {0.001f, 0.001f, 0.001f};
    arm_scale_f32(true_omega, M_PI / 180.0f, true_omega, 3);

    float gyro_noise = 0.001;
    float msmt_noise = 0.03;

    float state[6] = {0, 0, 0, 0, 0, 0};
    float cov[36], Q[36], R[36];
    eye_f32(cov, 6);
    arm_scale_f32(cov, 0.1, cov, 36);
    eye_f32(R, 6);
    arm_scale_f32(R, 0.1, R, 36);
    eye_f32(Q, 6);
    arm_scale_f32(Q, 0.01, Q, 36);
    // The legacy Q = 0.01 I (bias random walk of 0.1 rad/s per step) makes the bias
    // unobservable; use a small bias process noise so the bias can actually be estimated
    Q[21] = Q[28] = Q[35] = 1e-10f;

    float ref[6] = {40, 0, 0, 0, 40, 0};
    float ref_1[3] = {40, 0, 0};
    float ref_2[3] = {0, 40, 0};

    bool all_ok = true;
    float tail_err = 0;
    int tail_n = 0;
    for (int iter = 0; iter < 1000; iter++) {
        float simulated_gyro_measurement[3];
        for (int i = 0; i < 3; i++) {
            simulated_gyro_measurement[i] = true_omega[i] + true_bias[i] + uniform_pm1() * gyro_noise;
        }

        float delta_vec[3];
        memcpy(delta_vec, true_omega, sizeof(float) * 3);
        arm_scale_f32(delta_vec, dt, delta_vec, 3);
        float delta_q[4];
        rotationvec2quat(delta_vec, delta_q);
        quat_multiply(true_body_to_ref, delta_q, true_body_to_ref);
        quat_norm(true_body_to_ref, true_body_to_ref);
        float new_ref_to_body[4];
        quat_inv(true_body_to_ref, new_ref_to_body);

        float body_1[3], body_2[3];
        quat_apply(new_ref_to_body, ref_1, body_1);
        quat_apply(new_ref_to_body, ref_2, body_2);
        float body[6] = {body_1[0], body_1[1], body_1[2], body_2[0], body_2[1], body_2[2]};
        for (int i = 0; i < 6; i++) {
            body[i] = body[i] + uniform_pm1() * msmt_noise;
        }

        float new_error_state[6], new_quat[4], new_P[36];
        all_ok &= iterate(state, current_guess, cov, body, ref, 2, simulated_gyro_measurement, Q, R, dt,
                          new_error_state, new_quat, new_P)
                  == UKF_OK;
        memcpy(state, new_error_state, sizeof(float) * 6);
        memcpy(cov, new_P, sizeof(float) * 36);
        memcpy(current_guess, new_quat, sizeof(float) * 4);

        float err = RAD2DEG(quat_angle_between(true_body_to_ref, current_guess));
        if (iter % 200 == 0) {
            printf("  iter %4d: attitude error %.4f deg, bias [%.5f %.5f %.5f] deg/s\n", iter, err,
                   RAD2DEG(state[3]), RAD2DEG(state[4]), RAD2DEG(state[5]));
        }
        if (iter >= 800) {
            tail_err += err;
            tail_n++;
        }
    }
    float bias_err[3] = {state[3] - true_bias[0], state[4] - true_bias[1], state[5] - true_bias[2]};
    check_true(all_ok, "every step returns UKF_OK");
    check_less("mean error over the last 200 of 1000 steps (deg)", tail_err / tail_n, 0.2);
    check_less("gyro bias error after 1000 steps (deg/s)", RAD2DEG(l2_norm(bias_err, 3)), 0.02);
}

// =============================================================================================
// Control chain
// =============================================================================================
void test_bdot(void) {
    begin_suite("B-dot");
    float B1[3] = {2e-5f, -1e-5f, 3e-5f};
    float B0[3] = {1e-5f, 1e-5f, 3.5e-5f};
    float m[3];
    Bdot(B1, B0, 2.0f, 0.5f, m);
    float expected[3] = {-2.0f * (1e-5f) / 0.5f, -2.0f * (-2e-5f) / 0.5f, -2.0f * (-0.5e-5f) / 0.5f};
    check_vec_close("m = -k (B_t - B_t-1) / dt  (bDot.m)", m, expected, 3, 1e-10f);
    float bdot[3] = {(B1[0] - B0[0]) / 0.5f, (B1[1] - B0[1]) / 0.5f, (B1[2] - B0[2]) / 0.5f};
    check_true(dot3(m, bdot) < 0, "dipole opposes dB/dt (removes rotational energy)");

    Bdot(B1, B1, 2.0f, 0.5f, m);
    float zero[3] = {0, 0, 0};
    check_vec_close("constant field -> zero dipole", m, zero, 3, 0);

    Bdot(B1, B0, 2.0f, 0.0f, m);
    check_vec_close("dt = 0 -> zero dipole instead of inf", m, zero, 3, 0);
    Bdot(B1, B0, 2.0f, -1.0f, m);
    check_vec_close("negative dt -> zero dipole", m, zero, 3, 0);

    // The magnitude scales linearly with the gain
    float m1[3], m2[3];
    Bdot(B1, B0, 1.0f, 0.5f, m1);
    Bdot(B1, B0, 3.0f, 0.5f, m2);
    check_close("dipole scales with gain", m2[0] / m1[0], 3.0, 1e-5);
}

void test_pd(void) {
    begin_suite("PD controller (PD_loop.m)");
    float id[4] = {1, 0, 0, 0}, zero[3] = {0, 0, 0}, tau[3];
    pd_loop(id, zero, tau);
    check_vec_close("no error, no rate -> no torque", tau, zero, 3, 0);

    // Attitude error only: tau = I * Kp * theta * axis
    float rv[3] = {0.1f, 0, 0};
    float q_err[4];
    rotationvec2quat(rv, q_err);
    pd_loop(q_err, zero, tau);
    float expected[3] = {PD_INERTIA[0] * PD_KP * 0.1f, PD_INERTIA[3] * PD_KP * 0.1f, PD_INERTIA[6] * PD_KP * 0.1f};
    check_vec_close("0.1 rad about x -> I * Kp * 0.1 x", tau, expected, 3, 1e-8f);

    float q_neg[4] = {-q_err[0], -q_err[1], -q_err[2], -q_err[3]};
    float tau_neg[3];
    pd_loop(q_neg, zero, tau_neg);
    check_vec_close("q and -q give the same torque (shortest rotation)", tau_neg, tau, 3, 1e-8f);

    // Rate only: tau = -I * Kd * omega
    float omega[3] = {0.01f, -0.02f, 0.03f};
    pd_loop(id, omega, tau);
    float t[3] = {-PD_KD * omega[0], -PD_KD * omega[1], -PD_KD * omega[2]};
    float expected_rate[3];
    for (int i = 0; i < 3; i++) {
        expected_rate[i] = PD_INERTIA[i * 3] * t[0] + PD_INERTIA[i * 3 + 1] * t[1] + PD_INERTIA[i * 3 + 2] * t[2];
    }
    check_vec_close("rate only -> -I * Kd * omega", tau, expected_rate, 3, 1e-8f);
    // Damping should remove energy: omega . I^-1 tau < 0 is equivalent to omega . t < 0
    check_true(dot3(omega, t) < 0, "rate term opposes the rotation");

    // Matches the MATLAB formula on a random case (hand-expanded)
    float rv2[3] = {0.3f, -0.2f, 0.5f};
    float q2[4], om2[3] = {0.05f, 0.01f, -0.04f};
    rotationvec2quat(rv2, q2);
    pd_loop(q2, om2, tau);
    float t2[3], e2[3];
    for (int i = 0; i < 3; i++) {
        t2[i] = PD_KP * rv2[i] - PD_KD * om2[i];
    }
    for (int i = 0; i < 3; i++) {
        e2[i] = PD_INERTIA[i * 3] * t2[0] + PD_INERTIA[i * 3 + 1] * t2[1] + PD_INERTIA[i * 3 + 2] * t2[2];
    }
    check_vec_close("combined error + rate matches PD_loop.m", tau, e2, 3, 1e-6f);

    float huge[3] = {1e4f, -1e4f, 1e4f};
    pd_loop(id, huge, tau);
    bool clamped = true;
    for (int i = 0; i < 3; i++) {
        clamped &= fabsf(tau[i]) <= PD_MAX_TAU;
    }
    check_true(clamped, "torque saturates at +-PD_MAX_TAU");
}

// Magnetorquer constants from the hardware team (X, Y solid core; Z air core). Duplicated from
// moments2currents.c on purpose: changing them there without updating here should fail.
static const float MTQ_K[3] = {2.6f, 2.6f, 1.02f};         // A·m² per A
static const float MTQ_IMAX[3] = {0.156f, 0.156f, 0.104f}; // A

void test_torque_and_currents(void) {
    begin_suite("torque -> moment -> current");
    rng_seed(601);
    float worst_perp = 0, worst_torque = 0;
    for (int trial = 0; trial < 200; trial++) {
        float B[3], tau[3], m[3];
        random_unit_vec(B);
        float bmag = 2e-5f + 4e-5f * uniform_01();
        B[0] *= bmag;
        B[1] *= bmag;
        B[2] *= bmag;
        random_unit_vec(tau);
        tau[0] *= 1e-6f;
        tau[1] *= 1e-6f;
        tau[2] *= 1e-6f;
        torque_2_moments(B, tau, m);
        worst_perp = fmaxf(worst_perp, fabsf(dot3(m, B)) / (l2_norm(m, 3) * bmag + 1e-30f));
        // Produced torque m x B must equal the part of tau perpendicular to B
        float produced[3];
        cross(m, B, produced);
        float bhat[3] = {B[0] / bmag, B[1] / bmag, B[2] / bmag};
        float tb = dot3(tau, bhat);
        for (int i = 0; i < 3; i++) {
            float perp = tau[i] - tb * bhat[i];
            worst_torque = fmaxf(worst_torque, fabsf(produced[i] - perp) / 1e-6f);
        }
    }
    check_less("dipole is perpendicular to B (cosine)", worst_perp, 1e-5);
    check_less("m x B reproduces the achievable (perpendicular) torque (relative)", worst_torque, 1e-4);

    float B[3] = {0, 0, 3e-5f}, tau_par[3] = {0, 0, 1e-6f}, m[3], zero[3] = {0, 0, 0};
    torque_2_moments(B, tau_par, m);
    check_vec_close("torque parallel to B needs no dipole", m, zero, 3, 1e-12f);
    float B0[3] = {0, 0, 0}, tau[3] = {1, 2, 3};
    torque_2_moments(B0, tau, m);
    check_vec_close("zero field -> zero dipole (no divide by zero)", m, zero, 3, 0);
    // Explicit sign check: tau = +y, B = +z  ->  m must be +x (x cross z = -y... so check m x B = tau)
    float Bz[3] = {0, 0, 1}, ty[3] = {0, 1, 0}, mx[3], produced[3];
    torque_2_moments(Bz, ty, mx);
    cross(mx, Bz, produced);
    check_vec_close("sign: m x B == tau (torque2moment3axis.m has the opposite sign)", produced, ty, 3, 1e-7f);

    float I[3];
    float lim[3] = {1, 1, 1};
    float m_small[3] = {0.2f, -0.5f, 0.9f};
    moment2current3axis(m_small, lim, I);
    float expected_small[3] = {m_small[0] / MTQ_K[0], m_small[1] / MTQ_K[1], m_small[2] / MTQ_K[2]};
    check_vec_close("currents = m / k (k = 2.6 solid, 1.02 air)", I, expected_small, 3, 1e-7f);
    float m_big[3] = {5, -5, 0.5f};
    moment2current3axis(m_big, lim, I);
    float expected[3] = {1, -1, m_big[2] / MTQ_K[2]};
    check_vec_close("currents saturate at +-Imax", I, expected, 3, 1e-7f);
    float no_lim[3] = {-1, 1, -1};
    moment2current3axis(m_big, no_lim, I);
    float expected2[3] = {m_big[0] / MTQ_K[0], -1, m_big[2] / MTQ_K[2]};
    check_vec_close("negative Imax means unlimited on that axis", I, expected2, 3, 1e-7f);
    moment2current3axis(m_big, NULL, I);
    float expected3[3] = {MTQ_IMAX[0], -MTQ_IMAX[1], MTQ_IMAX[2]};
    check_vec_close("NULL Imax clamps to hardware limits (0.156 solid, 0.104 air)", I, expected3, 3, 0);
}

// =============================================================================================
// Orbits / time / environment models
// =============================================================================================
static float wrap180(float deg) {
    float w = fmodf(deg + 180.0f, 360.0f);
    if (w < 0) {
        w += 360.0f;
    }
    return w - 180.0f;
}

void test_kepler(void) {
    begin_suite("Kepler propagation");
    float out[6];

    // Curtis, Orbital Mechanics for Engineering Students (the MATLAB comment's "3.3371 rad" is a
    // typo; the double-precision algorithm and the textbook both give 193.16 deg = 3.3712 rad)
    propogateOrbitalElements((9600e3f + 21000e3f) / 2, 0.37255f, 0, 0, 0, 0, 10800.0f, out);
    check_close("Curtis example: true anomaly after 3 h (deg, python ref)", wrap180(out[5] - 193.155768543f), 0, 0.05);

    propogateOrbitalElements(7.0e6f, 0.001f, 51.6f, 30.0f, 45.0f, 10.0f, 0.0f, out);
    check_close("dt = 0 leaves the anomaly unchanged (deg)", wrap180(out[5] - 10.0f), 0, 0.01);
    check_true(exactly(out[0], 7.0e6f) && exactly(out[1], 0.001f) && exactly(out[2], 51.6f) && exactly(out[3], 30.0f)
                   && exactly(out[4], 45.0f),
               "a, e, i, RAAN, argp are passed through");

    float a = 7.0e6f;
    float period = 2.0f * (float)M_PI * sqrtf(a * a * a / 3.986004418e14f);
    propogateOrbitalElements(a, 0.1f, 0, 0, 0, 77.0f, period, out);
    check_close("one full period returns to the same anomaly (deg)", wrap180(out[5] - 77.0f), 0, 0.05);

    propogateOrbitalElements(7.0e6f, 0.001f, 0, 0, 0, 10.0f, 600.0f, out);
    check_close("near-circular LEO, 600 s (deg, python ref)", wrap180(out[5] - 47.123219969f), 0, 0.02);
    propogateOrbitalElements(1.0e7f, 0.7f, 0, 0, 0, 1.0f, 3600.0f, out);
    check_close("e = 0.7, 1 h (deg, python ref)", wrap180(out[5] - 167.267351552f), 0, 0.05);
    propogateOrbitalElements(6.9e6f, 0.05f, 0, 0, 0, 250.0f, 1234.5f, out);
    check_close("e = 0.05 from nu = 250 deg (deg, python ref)", wrap180(out[5] - 330.645500650f), 0, 0.02);
    propogateOrbitalElements(1.0e7f, 0.95f, 0, 0, 0, 1.0f, 5000.0f, out);
    check_true(isfinite(out[5]), "e = 0.95 stays finite");

    // Stepping 1 s at a time for 10 minutes agrees with one 600 s step (what body.c does)
    float nu = 10.0f;
    for (int i = 0; i < 600; i++) {
        propogateOrbitalElements(7.0e6f, 0.001f, 0, 0, 0, nu, 1.0f, out);
        nu = out[5];
    }
    check_close("600 x 1 s steps == one 600 s step (deg)", wrap180(nu - 47.123219969f), 0, 0.05);

    // Many orbits in one step doesn't lose precision (mean anomaly is wrapped)
    propogateOrbitalElements(a, 0.01f, 0, 0, 0, 20.0f, 100.0f * period, out);
    check_close("100 periods in one step returns to the same anomaly (deg)", wrap180(out[5] - 20.0f), 0, 0.5);

    propogateOrbitalElements(7.0e6f, 1.2f, 0, 0, 0, 10.0f, 60.0f, out);
    check_close("hyperbolic e >= 1 is left unchanged instead of producing NaN", out[5], 10.0, 0);
}

void test_orbital_to_eci(void) {
    begin_suite("orbital elements -> ECI (orbitalToECI.m)");
    // Example from orbitalToECI.m, checked there against MATLAB's keplerian2ijk (km, km/s)
    float kep_km[6] = {7000, 0.01f, 1.1f, 0.3f, 0.3f, 0.7f};
    float r[3], v[3];
    orbital_to_eci_posvel(kep_km, 398600.4418f, r, v);
    float r_expected[3] = {2801.905474f, 3641.952685f, 5209.109530f};
    float v_expected[3] = {-6.644010f, -0.085065f, 3.698018f};
    // The old inclination quaternion [cos(i/2), 0.5 sin(i), 0, 0] fails this by ~1000 km
    check_vec_close("position matches keplerian2ijk (km)", r, r_expected, 3, 0.05f);
    check_vec_close("velocity matches keplerian2ijk (km/s)", v, v_expected, 3, 1e-4f);

    float kep_m[6] = {7.0e6f, 0.01f, 1.1f, 0.3f, 0.3f, 0.7f};
    float r_m[3];
    orbital_to_eci(kep_m, r_m);
    float r_expected_m[3] = {2801905.474f, 3641952.685f, 5209109.530f};
    check_vec_close("position-only wrapper in metres", r_m, r_expected_m, 3, 50.0f);

    // Physics checks on random orbits
    rng_seed(701);
    float worst_radius = 0, worst_energy = 0, worst_rv_h = 0;
    for (int trial = 0; trial < 100; trial++) {
        float k[6] = {6.7e6f + 3e6f * uniform_01(), 0.3f * uniform_01(), (float)M_PI * uniform_01(),
                      2 * (float)M_PI * uniform_01(), 2 * (float)M_PI * uniform_01(), 2 * (float)M_PI * uniform_01()};
        float rr[3], vv[3];
        orbital_to_eci_posvel(k, MU_EARTH_M3_S2, rr, vv);
        float rmag = l2_norm(rr, 3), vmag = l2_norm(vv, 3);
        float r_expected_mag = k[0] * (1 - k[1] * k[1]) / (1 + k[1] * cosf(k[5]));
        worst_radius = fmaxf(worst_radius, fabsf(rmag - r_expected_mag) / r_expected_mag);
        // vis-viva: v^2/2 - mu/r = -mu/(2a)
        float energy = 0.5f * vmag * vmag - MU_EARTH_M3_S2 / rmag;
        float expected_energy = -MU_EARTH_M3_S2 / (2 * k[0]);
        worst_energy = fmaxf(worst_energy, fabsf(energy - expected_energy) / fabsf(expected_energy));
        // angular momentum h = r x v points along the orbit normal: z-component = |h| cos(i)
        float h[3];
        cross(rr, vv, h);
        worst_rv_h = fmaxf(worst_rv_h, fabsf(h[2] / l2_norm(h, 3) - cosf(k[2])));
    }
    check_less("|r| = a(1-e^2)/(1+e cos nu) (relative)", worst_radius, 1e-5);
    check_less("vis-viva energy (relative)", worst_energy, 1e-4);
    check_less("orbit normal makes angle i with z (cosine)", worst_rv_h, 1e-4);

    float equatorial[6] = {7.0e6f, 0.0f, 0.0f, 0.5f, 0.2f, 1.0f};
    orbital_to_eci_posvel(equatorial, MU_EARTH_M3_S2, r, v);
    check_close("i = 0 -> z = 0", r[2], 0, 1e-2);
    check_close("i = 0 -> vz = 0", v[2], 0, 1e-5);
    check_close("i = 0: angle from x axis = RAAN + argp + nu (rad)", atan2f(r[1], r[0]), 1.7, 1e-5);
}

void test_ecef_to_eci(void) {
    begin_suite("ECEF -> ECI (eceftoeci.m)");
    check_close("unix 0 = JD 2440587.5", unix_2_jd(0), 2440587.5, 0);
    check_close("J2000 epoch (unix 946728000) = JD 2451545", unix_2_jd(946728000), 2451545.0, 1e-9);
    check_close("GMST at J2000 (deg)", jd_2_gmst_deg(2451545.0), 280.46061837, 1e-6);
    double g = jd_2_gmst_deg(2461116.25);
    check_true(g >= 0 && g < 360, "GMST is reduced to [0, 360)");

    float ex[3] = {1, 0, 0}, out[3];
    ecef_2_eci(ex, out, 946728000);
    float expected_j2000[3] = {cosf(DEG2RAD(280.46061837f)), sinf(DEG2RAD(280.46061837f)), 0};
    check_vec_close("ECEF x axis at J2000", out, expected_j2000, 3, 1e-5f);

    // Providence at a 2026 time (python ref, double). The old float32 JD was off by hours here.
    float pvd[3] = {6.3761e6f, -0.1387e6f, 0.0807e6f};
    ecef_2_eci(pvd, out, 1772814426);
    float expected[3] = {4103984.46f, 4881721.05f, 80700.0f};
    check_vec_close("Providence at unix 1772814426 (m, python ref)", out, expected, 3, 10.0f);
    check_close("rotation preserves |r| (m)", l2_norm(out, 3), l2_norm(pvd, 3), 1.0);
    check_close("rotation preserves z", out[2], pvd[2], 0);

    // Consistency with the (different) GMST polynomial inside the WMM code
    int unix_times[3] = {946728000, 1772814426, 1790000000};
    double worst = 0;
    for (int i = 0; i < 3; i++) {
        double jd = unix_2_jd(unix_times[i]);
        double a = jd_2_gmst_deg(jd) * M_PI / 180.0;
        double b = gmst_from_jd(jd);
        double d = fabs(remainder(a - b, 2 * M_PI));
        worst = fmax(worst, d);
    }
    check_less("agrees with magnetosphere's gmst_from_jd (rad)", worst, 1e-5);

    // One sidereal day later the rotation is (almost) the same; half a solar day is ~180 deg off
    float r1[3], r2[3];
    ecef_2_eci(pvd, r1, 1772814426);
    ecef_2_eci(pvd, r2, 1772814426 + 86164);
    check_less("one sidereal day later -> same ECI position (angle, deg)", RAD2DEG(vec_angle(r1, r2)), 0.01);
    ecef_2_eci(ex, r1, 1772814426);
    ecef_2_eci(ex, r2, 1772814426 + 43082);
    check_close("half a sidereal day later -> x axis flipped (deg)", RAD2DEG(vec_angle(r1, r2)), 180, 0.02);
}

void test_sun_vec(void) {
    begin_suite("sun vector (sunVectorECI.m)");
    float s[3];

    // March equinox 2024-03-20 03:06 UTC: the Sun is at the vernal equinox direction (+x)
    sun_vec(1710903960, s);
    float plus_x[3] = {1, 0, 0};
    check_less("March 2024 equinox -> +x (deg)", RAD2DEG(vec_angle(s, plus_x)), 0.05);
    // June solstice 2024-06-20 20:51 UTC: ecliptic longitude 90 deg -> (0, cos eps, sin eps)
    sun_vec(1718916660, s);
    float solstice[3] = {0, cosf(DEG2RAD(23.4393f)), sinf(DEG2RAD(23.4393f))};
    check_less("June 2024 solstice -> (0, cos eps, sin eps) (deg)", RAD2DEG(vec_angle(s, solstice)), 0.05);

    sun_vec(1767225600, s);
    float ref1[3] = {0.183377229f, -0.901947043f, -0.390978673f};
    check_vec_close("2026-01-01 (python ref)", s, ref1, 3, 2e-6f);
    sun_vec(1790000000, s);
    float ref2[3] = {-0.999706841f, 0.0222148713f, 0.0096297249f};
    check_vec_close("2026-09-21 (python ref)", s, ref2, 3, 2e-6f);

    float worst_norm = 0, worst_ecliptic = 0;
    for (int k = 0; k < 365; k++) {
        sun_vec(1767225600 + k * 86400, s);
        worst_norm = fmaxf(worst_norm, fabsf(l2_norm(s, 3) - 1));
        // Always in the ecliptic plane: s . n_ecliptic = 0, n = (0, -sin eps, cos eps)
        float n[3] = {0, -sinf(DEG2RAD(23.4393f)), cosf(DEG2RAD(23.4393f))};
        worst_ecliptic = fmaxf(worst_ecliptic, fabsf(dot3(s, n)));
    }
    check_less("unit length every day for a year", worst_norm, 1e-6);
    check_less("stays in the ecliptic plane", worst_ecliptic, 1e-4);

    float s1[3], s2[3];
    sun_vec(1767225600, s1);
    sun_vec(1767225600 + 86400, s2);
    check_close("moves ~0.99 deg per day", RAD2DEG(vec_angle(s1, s2)), 1.0, 0.05);
}

void test_magnetosphere(void) {
    begin_suite("WMM2025 magnetic field");

    check_close("jd2year(J2000) = 2000", jd2year(2451545.0), 2000.0, 1e-4);
    check_close("jd2year(J2000 + 365.25 * 25) = 2025", jd2year(2451545.0 + 365.25 * 25), 2025.0, 1e-3);
    check_close("gmst_from_jd(J2000) (rad, python ref)", gmst_from_jd(2451545.0), 4.894961213, 1e-6);
    check_close("gmst_from_jd(2461116.25) (rad, python ref)", gmst_from_jd(2461116.25), 1.471975777, 1e-6);

    // Geodetic conversion round trip against the closed-form geodetic -> ECEF map
    double a = 6378137.0, f = 1 / 298.257223563, e2 = f * (2 - f);
    double lats[7] = {-89.5, -45, -10, 0, 30, 60, 89.9};
    double alts[2] = {0, 550e3};
    float worst_lat = 0, worst_lon = 0, worst_alt = 0;
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 2; j++) {
            double lat = lats[i] * M_PI / 180, lon = (-170 + 50 * i) * M_PI / 180, h = alts[j];
            double N = a / sqrt(1 - e2 * sin(lat) * sin(lat));
            float r[3] = {(float)((N + h) * cos(lat) * cos(lon)), (float)((N + h) * cos(lat) * sin(lon)),
                          (float)((N * (1 - e2) + h) * sin(lat))};
            float lat_o, lon_o, alt_o;
            ecef_to_geodetic(r, &lat_o, &lon_o, &alt_o);
            worst_lat = fmaxf(worst_lat, (float)fabs(lat_o - lat));
            worst_lon = fmaxf(worst_lon, (float)fabs(remainder(lon_o - lon, 2 * M_PI)));
            worst_alt = fmaxf(worst_alt, (float)fabs(alt_o - h));
        }
    }
    check_less("ecef_to_geodetic latitude round trip (deg)", RAD2DEG(worst_lat), 1e-3);
    check_less("ecef_to_geodetic longitude round trip (deg)", RAD2DEG(worst_lon), 1e-3);
    check_less("ecef_to_geodetic altitude round trip (m)", worst_alt, 5.0);
    float pole[3] = {0, 0, 6356752.3f + 1000.0f}, lat_p, lon_p, alt_p;
    ecef_to_geodetic(pole, &lat_p, &lon_p, &alt_p);
    check_close("north pole latitude", lat_p, M_PI / 2, 1e-6);
    check_close("north pole altitude (m)", alt_p, 1000, 1.0);

    // Field vs an independent reference: B = -grad(V) of the spherical harmonic potential,
    // evaluated in double by tools/reference_values.py (no geodetic conversion, no NED frame)
    struct {
        float r[3];
        int jd_int;
        float jd_frac;
        float b[3];
        const char* name;
    } cases[5] = {
        {{6878137, 0, 0}, 2461116, 0.25f, {-7.94777569e-06f, -4.26611433e-07f, 2.14833334e-05f}, "equatorial, 500 km"},
        {{0, -4200000, 5400000}, 2461116, 0.25f, {1.41981506e-06f, 3.57564235e-05f, -1.84210689e-05f}, "northern mid-latitude"},
        {{3100000, 2200000, -5600000}, 2461300, 0.7f, {1.91990743e-05f, 1.15574319e-05f, -1.01150619e-05f}, "southern mid-latitude"},
        {{-4000000, 1000000, 5300000}, 2460900, 0.1f, {3.59249385e-05f, -7.92033649e-06f, -2.02947473e-05f}, "2025 epoch, north"},
        {{1000, 2000, 6900000}, 2461000, 0.5f, {-5.74409126e-07f, -8.64492328e-07f, -4.55131861e-05f}, "near the north pole"},
    };
    for (int i = 0; i < 5; i++) {
        float b[3];
        wmm_eci_embedded_v2(cases[i].r, cases[i].jd_int, cases[i].jd_frac, b);
        float diff[3] = {b[0] - cases[i].b[0], b[1] - cases[i].b[1], b[2] - cases[i].b[2]};
        char name[128];
        snprintf(name, sizeof(name), "%s: |B - B_ref| / |B_ref|", cases[i].name);
        check_less(name, l2_norm(diff, 3) / l2_norm(cases[i].b, 3), 2e-4);
        snprintf(name, sizeof(name), "%s: direction error (deg)", cases[i].name);
        check_less(name, RAD2DEG(vec_angle(b, cases[i].b)), 0.01);
    }

    // Physical sanity: strength at LEO and ~r^-3 fall-off
    float r_leo[3] = {6878137, 0, 0}, r_hi[3] = {2 * 6878137.0f, 0, 0}, b_leo[3], b_hi[3];
    wmm_eci_embedded_v2(r_leo, 2461116, 0.25f, b_leo);
    wmm_eci_embedded_v2(r_hi, 2461116, 0.25f, b_hi);
    float mag = l2_norm(b_leo, 3);
    check_true(mag > 1.5e-5f && mag < 6.5e-5f, "field strength at 500 km is 15-65 uT");
    check_close("doubling the radius divides |B| by ~8 (dipole)", l2_norm(b_leo, 3) / l2_norm(b_hi, 3), 8.0, 1.0);
}

// =============================================================================================
// Sensors
// =============================================================================================
static void ideal_photodiode_readings(const float* sun_body, float scale, float noise, float* readings) {
    for (int i = 0; i < NUM_DIODES; i++) {
        float d = dot3(PHOTODIODES[i], sun_body);
        float v = scale * (d > 0 ? d : 0);
        if (noise > 0) {
            v += normal_sample(0, noise);
        }
        readings[i] = v > 0 ? v : 0;
    }
}

void test_photodiodes(void) {
    begin_suite("photodiode sun vector");

    bool unit = true, pairs = true;
    for (int i = 0; i < NUM_DIODES; i++) {
        unit &= fabsf(l2_norm(PHOTODIODES[i], 3) - 1) < 1e-5f;
        if (i % 2 == 0) {
            float s[3] = {PHOTODIODES[i][0] + PHOTODIODES[i + 1][0], PHOTODIODES[i][1] + PHOTODIODES[i + 1][1],
                          PHOTODIODES[i][2] + PHOTODIODES[i + 1][2]};
            pairs &= l2_norm(s, 3) < 1e-6f;
        }
    }
    check_true(unit, "diode normals are unit vectors");
    check_true(pairs, "diodes 2k and 2k+1 point in opposite directions");

    rng_seed(801);
    float worst = 0;
    int in_sun = 0, n_trials = 1000;
    for (int trial = 0; trial < n_trials; trial++) {
        float s[3], readings[NUM_DIODES], est[3];
        random_unit_vec(s);
        ideal_photodiode_readings(s, MAX_READING, 0, readings);
        if (get_vec_from_photodiode_readings(readings, est)) {
            in_sun++;
            worst = fmaxf(worst, RAD2DEG(vec_angle(est, s)));
        }
    }
    check_true(in_sun == n_trials, "ideal readings: always reports in-sun");
    check_less("ideal readings: max error over 1000 directions (deg)", worst, 0.05);

    rng_seed(802);
    float sum = 0;
    int n_ok = 0;
    for (int trial = 0; trial < n_trials; trial++) {
        float s[3], readings[NUM_DIODES], est[3];
        random_unit_vec(s);
        ideal_photodiode_readings(s, MAX_READING, 0.025f * MAX_READING, readings);
        if (get_vec_from_photodiode_readings(readings, est)) {
            n_ok++;
            sum += RAD2DEG(vec_angle(est, s));
        }
    }
    check_true(n_ok > n_trials * 0.99, "2.5% noise: still reports in-sun");
    check_less("2.5% noise: mean error (deg)", sum / (n_ok ? n_ok : 1), 3.0);

    // Output is a unit vector regardless of the brightness scale
    float s[3] = {0.2f, 0.9f, -0.3f}, readings[NUM_DIODES], est[3];
    normalize_vec(s, 3);
    ideal_photodiode_readings(s, 0.8f * MAX_READING, 0, readings);
    get_vec_from_photodiode_readings(readings, est);
    check_close("dimmer sun still gives a unit vector", l2_norm(est, 3), 1, 1e-5);
    check_less("dimmer sun: direction still correct (deg)", RAD2DEG(vec_angle(est, s)), 0.05);

    float dark[NUM_DIODES] = {0};
    check_true(!get_vec_from_photodiode_readings(dark, est), "eclipse (all zero) -> not in sun");
    float two[NUM_DIODES] = {0};
    two[0] = 1.0f;
    two[4] = 1.0f;
    check_true(!get_vec_from_photodiode_readings(two, est), "only two lit pairs -> not in sun");
    float nan_r[NUM_DIODES];
    ideal_photodiode_readings(s, MAX_READING, 0, nan_r);
    nan_r[3] = NAN;
    check_true(!get_vec_from_photodiode_readings(nan_r, est), "NaN reading -> rejected");
}

// =============================================================================================
// Pointing
// =============================================================================================
void test_pointing(void) {
    begin_suite("pointing (pointing_error.m / down_quat)");
    float ez[3] = {0, 0, 1}, ey[3] = {0, 1, 0}, ex[3] = {1, 0, 0};

    float r[3] = {7.0e6f, 0, 0}, v[3] = {0, 7500, 0}, target[3] = {6.37e6f, 0, 0};
    float id[4] = {1, 0, 0, 0}, q_err[4], z_want[3];
    float prev[4] = {1, 0, 0, 0};
    bool ok = pointing_error(r, v, id, target, prev, q_err, z_want);
    check_true(ok, "normal geometry succeeds");
    float nadir[3] = {-1, 0, 0};
    check_vec_close("z_want points at the target", z_want, nadir, 3, 1e-6f);

    // With current attitude = identity, q_err == q_want
    float bz[3], by[3], bx[3];
    quat_apply(q_err, ez, bz);
    quat_apply(q_err, ey, by);
    quat_apply(q_err, ex, bx);
    check_vec_close("desired body z points at the target", bz, z_want, 3, 1e-5f);
    check_close("desired body y is perpendicular to velocity", dot3(by, v) / l2_norm(v, 3), 0, 1e-5);
    float yz[3];
    cross(by, bz, yz);
    check_vec_close("desired frame is right-handed (x = y cross z)", bx, yz, 3, 1e-5f);

    // Already at the desired attitude -> identity error, and the PD produces no torque
    float q_want[4];
    memcpy(q_want, q_err, sizeof(q_want));
    ok = pointing_error(r, v, q_want, target, prev, q_err, NULL);
    float rv[3];
    quat2rotationvec(q_err, rv);
    check_less("at the goal the error rotation is ~0 (deg)", RAD2DEG(l2_norm(rv, 3)), 1e-3);

    // Error is in BODY coordinates: rotating it into ECI with the current attitude reproduces
    // q_want (q_b2eci * q_err == q_want)
    float q_cur[4] = {0.9f, 0.2f, -0.3f, 0.1f};
    quat_norm(q_cur, q_cur);
    pointing_error(r, v, q_cur, target, prev, q_err, NULL);
    float recon[4];
    quat_multiply(q_cur, q_err, recon);
    check_less("q_b2eci * q_err == q_want (body-frame error)", RAD2DEG(quat_angle_between(recon, q_want)), 1e-3);

    // Hemisphere continuity via q_want_prev
    float prev_neg[4] = {-q_want[0], -q_want[1], -q_want[2], -q_want[3]};
    float q_err2[4];
    pointing_error(r, v, id, target, prev_neg, q_err2, NULL);
    check_true(dot3(&q_err2[1], &prev_neg[1]) + q_err2[0] * prev_neg[0] > 0,
               "q_want is kept in the same hemisphere as q_want_prev");

    float v_par[3] = {-1, 0, 0};
    ok = pointing_error(r, v_par, id, target, prev, q_err, NULL);
    check_true(!ok, "velocity along the line of sight is rejected");
    check_vec_close("degenerate geometry gives identity error", q_err, id, 4, 0);
    ok = pointing_error(r, v, id, r, prev, q_err, NULL);
    check_true(!ok, "target at the satellite position is rejected");

    // down_quat: body z toward target, right-handed, minimal roll
    float goal[4];
    down_quat(r, target, id, goal);
    quat_apply(goal, ez, bz);
    check_vec_close("down_quat: body z points at the target (old version pointed away)", bz, nadir, 3, 1e-5f);
    float Rg[9];
    quat2rotm(goal, Rg);
    float det = Rg[0] * (Rg[4] * Rg[8] - Rg[5] * Rg[7]) - Rg[1] * (Rg[3] * Rg[8] - Rg[5] * Rg[6])
                + Rg[2] * (Rg[3] * Rg[7] - Rg[4] * Rg[6]);
    check_close("down_quat: proper rotation (det = +1; old version built a det = -1 frame)", det, 1, 1e-5);
    quat_apply(goal, ey, by);
    check_vec_close("down_quat: keeps the current body y axis when possible", by, ey, 3, 1e-5f);
    float target_along_y[3] = {7.0e6f, 1.0e6f, 0};
    down_quat(r, target_along_y, id, goal);
    quat_apply(goal, ez, bz);
    check_vec_close("down_quat: degenerate (y parallel to target) still points z at the target", bz, ey, 3, 1e-5f);
    down_quat(r, r, id, goal);
    check_vec_close("down_quat: coincident points give identity", goal, id, 4, 0);
}

// =============================================================================================
// Filters (filters.c was not compiled before)
// =============================================================================================
void test_filters(void) {
    begin_suite("signal filters");
    float src[20], dst[20];
    for (int i = 0; i < 20; i++) {
        src[i] = 3.0f;
    }
    gaussian_smooth_1d(src, dst, 20, 1.5f);
    check_vec_close("gaussian: constant signal is preserved", dst, src, 20, 1e-5f);

    for (int i = 0; i < 20; i++) {
        src[i] = 0.5f * i - 2;
    }
    gaussian_smooth_1d(src, dst, 20, 1.0f);
    check_vec_close("gaussian: linear ramp preserved away from the edges", &dst[4], &src[4], 12, 1e-4f);

    float impulse[21] = {0};
    impulse[10] = 1;
    float out21[21];
    gaussian_smooth_1d(impulse, out21, 21, 1.0f);
    float s = 0;
    for (int i = 0; i < 21; i++) {
        s += out21[i];
    }
    check_close("gaussian: kernel sums to 1", s, 1, 1e-5);
    check_true(out21[10] > out21[9] && out21[9] > out21[8] && fabsf(out21[9] - out21[11]) < 1e-7f,
               "gaussian: impulse response is symmetric and peaked");

    // Signals shorter than the kernel radius used to index outside the array
    float tiny_src[2] = {1.0f, 1.0f}, tiny_dst[2];
    gaussian_smooth_1d(tiny_src, tiny_dst, 2, 3.0f);
    check_vec_close("gaussian: n = 2 with a wide kernel stays in bounds", tiny_dst, tiny_src, 2, 1e-6f);
    float one_src[1] = {4.0f}, one_dst[1];
    gaussian_smooth_1d(one_src, one_dst, 1, 3.0f);
    check_close("gaussian: n = 1 is a copy", one_dst[0], 4.0, 0);

    for (int i = 0; i < 20; i++) {
        src[i] = (float)(i % 3);
    }
    exponential_filter(src, dst, 20, 1.0f);
    check_vec_close("exponential: alpha = 1 is a copy", dst, src, 20, 0);
    exponential_filter(src, dst, 20, 0.0f);
    bool all_first = true;
    for (int i = 0; i < 20; i++) {
        all_first &= exactly(dst[i], src[0]);
    }
    check_true(all_first, "exponential: alpha = 0 holds the first sample");
    float step[20];
    for (int i = 0; i < 20; i++) {
        step[i] = i < 5 ? 0.0f : 1.0f;
    }
    exponential_filter(step, dst, 20, 0.3f);
    check_close("exponential: step response after 15 samples", dst[19], 1 - powf(0.7f, 15), 1e-5);

    for (int i = 0; i < 20; i++) {
        src[i] = 2.0f + 0.25f * i;
    }
    holtz_double_exp_filter(src, dst, 20, 0.4f, 0.3f);
    check_vec_close("holt: tracks a linear trend exactly", dst, src, 20, 1e-5f);
    // n <= 0 must not touch dst
    dst[0] = 42;
    holtz_double_exp_filter(src, dst, 0, 0.4f, 0.3f);
    exponential_filter(src, dst, 0, 0.4f);
    gaussian_smooth_1d(src, dst, 0, 1.0f);
    check_close("n = 0 writes nothing", dst[0], 42, 0);
}

// =============================================================================================
// Full ADCS step (body.c)
// =============================================================================================
typedef struct {
    float q_true[4]; // body -> ECI
    float omega[3];
    float gyro_bias[3];
    float tle_mean[6]; // n (rev/day), e, i, RAAN, argp, M (rad)
    double tle_epoch_jd;
    float tle_bstar;
    bool tle_sent;
    int unix_time;
    int jd_int;
    double jd_frac;
    float last_mag[3];
} body_sim_t;

static void body_sim_init(body_sim_t* s) {
    memset(s, 0, sizeof(*s));
    s->q_true[0] = 0.7f;
    s->q_true[1] = 0.1f;
    s->q_true[2] = -0.5f;
    s->q_true[3] = 0.5f;
    quat_norm(s->q_true, s->q_true);
    double mean_motion = sqrt(3.986004418e14 / (6.9e6 * 6.9e6 * 6.9e6)) * 86400.0
                         / (2.0 * M_PI);
    float tle[6] = {(float)mean_motion, 0.001f, DEG2RAD(51.6f), DEG2RAD(40.0f),
                    DEG2RAD(10.0f), DEG2RAD(30.0f)};
    memcpy(s->tle_mean, tle, sizeof(tle));
    s->unix_time = 1767225600; // 2026-01-01 00:00:00 UTC
    double jd = unix_2_jd(s->unix_time);
    s->jd_int = (int)floor(jd);
    s->jd_frac = jd - s->jd_int;
    s->tle_epoch_jd = jd;
    s->tle_bstar = 0.0;
}

// Sensor readings consistent with the truth state
static void body_sim_sensors(body_sim_t* s, bool sunlit, float* mag, float* gyro, float* diodes) {
    float r[3], b_eci[3], sun_eci[3], sun_body[3];
    float r_km[3], v_km_s[3], oe_osc[7];
    double elapsed = ((double)s->jd_int + s->jd_frac - s->tle_epoch_jd) * 86400.0;
    int orbit_status = sgp4_step(s->tle_mean, s->tle_bstar, s->tle_epoch_jd, (float)elapsed,
                                 r_km, v_km_s, oe_osc);
    if (orbit_status != 0) {
        memset(mag, 0, sizeof(float) * 3);
        memset(gyro, 0, sizeof(float) * 3);
        memset(diodes, 0, sizeof(float) * NUM_DIODES);
        return;
    }
    for (int i = 0; i < 3; i++) r[i] = (float)(1000.0 * r_km[i]);
    wmm_eci_embedded_v2(r, s->jd_int, (float)s->jd_frac, b_eci);
    ref_rotate_inv(s->q_true, b_eci, mag);
    for (int i = 0; i < 3; i++) {
        mag[i] += normal_sample(0, 5e-8f); // 50 nT noise
        gyro[i] = s->omega[i] + s->gyro_bias[i] + normal_sample(0, 2e-4f);
    }
    sun_vec(s->unix_time, sun_eci);
    ref_rotate_inv(s->q_true, sun_eci, sun_body);
    if (sunlit) {
        ideal_photodiode_readings(sun_body, MAX_READING, 0.01f * MAX_READING, diodes);
    } else {
        memset(diodes, 0, sizeof(float) * NUM_DIODES);
    }
}

static void body_sim_advance(body_sim_t* s, int dt_s) {
    float rv[3] = {s->omega[0] * dt_s, s->omega[1] * dt_s, s->omega[2] * dt_s};
    float dq[4];
    rotationvec2quat(rv, dq);
    quat_multiply(s->q_true, dq, s->q_true);
    quat_norm(s->q_true, s->q_true);
    s->unix_time += dt_s;
    s->jd_frac += dt_s / 86400.0;
    while (s->jd_frac >= 1.0) {
        s->jd_frac -= 1.0;
        s->jd_int++;
    }
}

// One body() call with SGP4-consistent sensors; upload the TLE only once.
static void body_sim_step(body_sim_t* s, bool sunlit, int dt_s, float* currents) {
    body_sim_advance(s, dt_s);
    float mag[3], gyro[3], diodes[NUM_DIODES];
    body_sim_sensors(s, sunlit, mag, gyro, diodes);
    const float* update = NULLPTR;
    if (!s->tle_sent) {
        body_set_tle_metadata(s->tle_bstar, s->tle_epoch_jd);
        update = s->tle_mean;
        s->tle_sent = true;
    }
    body(s->last_mag, mag, gyro, diodes, update, s->unix_time, s->jd_int,
         (float)s->jd_frac, (float)dt_s, currents);
    memcpy(s->last_mag, mag, sizeof(mag));
}

static float body_attitude_error_deg(const body_sim_t* s) {
    float q[4];
    body_get_attitude(q);
    return RAD2DEG(quat_angle_between(q, s->q_true));
}

void test_body(void) {
    begin_suite("ADCS step (body.c)");
    float currents[3], zero[3] = {0, 0, 0};
    float zero_moment_currents[3];
    moment2current3axis(zero, Imax, zero_moment_currents);
    float diodes_dark[NUM_DIODES] = {0};
    float mag[3] = {2e-5f, -1e-5f, 3e-5f}, mag_prev[3] = {2.1e-5f, -1.2e-5f, 2.9e-5f}, gyro[3] = {0.3f, -0.2f, 0.1f};

    body_reset();
    currents[0] = currents[1] = currents[2] = 123.0f;
    body(mag_prev, mag, gyro, diodes_dark, NULLPTR, 1767225600, 2461041, 0.5f, 1.0f, currents);
    check_vec_close("no position yet -> zero currents (and currents are always written)", currents, zero, 3, 0);
    check_true(!body_is_pointing(), "no position yet -> not pointing");
    {
        float tle_out[6], r_out[3], bstar_out;
        double epoch_out;
        check_true(!body_get_tle(tle_out, &bstar_out, &epoch_out), "no TLE yet -> body_get_tle returns false");
        check_true(!body_get_position(r_out), "no position yet -> body_get_position returns false");

        body_sim_t getter_sim;
        body_sim_init(&getter_sim);
        body_reset();
        body_set_tle_metadata(getter_sim.tle_bstar, getter_sim.tle_epoch_jd);
        body(mag_prev, mag, gyro, diodes_dark, getter_sim.tle_mean, getter_sim.unix_time, getter_sim.jd_int,
             (float)getter_sim.jd_frac, 0.5f, currents);
        check_true(body_get_tle(tle_out, &bstar_out, &epoch_out), "committed TLE -> body_get_tle returns true");
        check_true(memcmp(tle_out, getter_sim.tle_mean, sizeof(tle_out)) == 0, "body_get_tle returns the committed elements");
        check_true(memcmp(&epoch_out, &getter_sim.tle_epoch_jd, sizeof(epoch_out)) == 0, "body_get_tle returns the committed epoch");
        check_true(body_get_position(r_out) && l2_norm(r_out, 3) > 6.0e6f, "body_get_position returns an orbital radius (m)");
        body_reset();
        check_true(!body_get_tle(tle_out, &bstar_out, &epoch_out), "body_reset clears the TLE");
        check_true(!body_get_position(r_out), "body_reset clears the position");
    }

    // Detumbling: fast spin, eclipse -> B-dot
    body_sim_t s;
    body_sim_init(&s);
    body_reset();
    body(mag_prev, mag, gyro, diodes_dark, s.tle_mean, s.unix_time, s.jd_int, (float)s.jd_frac, 0.5f, currents);
    check_vec_close("TLE without metadata (epoch, B*) is rejected -> zero currents", currents, zero, 3, 0);
    body_set_tle_metadata(s.tle_bstar, s.tle_epoch_jd);
    body(mag_prev, mag, gyro, diodes_dark, s.tle_mean, s.unix_time, s.jd_int, (float)s.jd_frac, 0.5f, currents);
    float dipole[3], expected[3];
    for (int i = 0; i < 3; i++) {
        dipole[i] = -BDOT_GAIN * (mag[i] - mag_prev[i]) / 0.5f;
    }
    moment2current3axis(dipole, Imax, expected);
    check_true(!body_is_pointing(), "tumbling in eclipse -> stays in detumble mode");
    check_vec_close("detumble output is the B-dot dipole (-k dB/dt) through the current model", currents, expected, 3, 1e-10f);

    // Garbage in never crashes and never produces garbage out
    float nan_mag[3] = {NAN, 0, 0};
    body(mag_prev, nan_mag, gyro, diodes_dark, s.tle_mean, s.unix_time, s.jd_int, (float)s.jd_frac, 0.5f, currents);
    check_vec_close("NaN magnetometer -> zero currents", currents, zero, 3, 0);
    body(mag_prev, mag, gyro, diodes_dark, s.tle_mean, s.unix_time, s.jd_int, (float)s.jd_frac, 0.0f, currents);
    check_vec_close("dt = 0 -> zero currents", currents, zero, 3, 0);
    body(mag_prev, mag, gyro, diodes_dark, s.tle_mean, s.unix_time, s.jd_int, (float)s.jd_frac, NAN, currents);
    check_vec_close("dt = NaN -> zero currents", currents, zero, 3, 0);
    body(NULLPTR, mag, gyro, diodes_dark, s.tle_mean, s.unix_time, s.jd_int, (float)s.jd_frac, 0.5f, currents);
    check_vec_close("NULL previous magnetometer -> zero currents", currents, zero, 3, 0);

    // Slow rotation in sunlight -> QUEST initialization and pointing
    rng_seed(901);
    body_sim_init(&s);
    s.omega[0] = 0.002f;
    s.omega[1] = -0.001f;
    s.omega[2] = 0.0015f;
    s.gyro_bias[0] = 0.001f;
    s.gyro_bias[1] = -0.0005f;
    s.gyro_bias[2] = 0.0008f;
    body_reset();
    {
        float m0[3], g0[3], d0[NUM_DIODES];
        body_sim_sensors(&s, true, m0, g0, d0);
        memcpy(s.last_mag, m0, sizeof(m0));
    }
    body_sim_step(&s, true, 1, currents);
    check_true(body_is_pointing(), "low rates + sun visible -> switches to pointing");
    float q_dummy[4];
    check_true(body_get_attitude(q_dummy), "attitude is initialized");
    check_less("QUEST initial attitude error (deg)", body_attitude_error_deg(&s), 5.0);
    check_vec_close("transition step commands zero moment", currents, zero_moment_currents, 3, 0);

    // Sunlit tracking for 20 minutes
    bool finite = true, bounded = true, any_nonzero = false, quiet_window_ok = true;
    float prior_bias[3], max_predict_bias_step = 0.0f, max_update_bias_step = 0.0f;
    body_get_gyro_bias(prior_bias);
    float tail = 0;
    for (int k = 0; k < 1200; k++) {
        body_sim_step(&s, true, 1, currents);
        unsigned int cycle_tick = (unsigned int)k % BODY_MEASUREMENT_UPDATE_PERIOD_TICKS + 1u;
        finite &= all_finite(currents, 3);
        for (int i = 0; i < 3; i++) {
            bounded &= fabsf(currents[i]) <= MTQ_IMAX[i];
            any_nonzero |= !exactly(currents[i], 0);
            if (cycle_tick >= BODY_MAGNETORQUER_QUIET_START_TICK) {
                quiet_window_ok &= exactly(currents[i], zero_moment_currents[i]);
            }
        }
        float current_bias[3], bias_step[3];
        body_get_gyro_bias(current_bias);
        for (int i = 0; i < 3; i++) bias_step[i] = current_bias[i] - prior_bias[i];
        float bias_change = l2_norm(bias_step, 3);
        if (cycle_tick == BODY_MEASUREMENT_UPDATE_PERIOD_TICKS) {
            max_update_bias_step = fmaxf(max_update_bias_step, bias_change);
        } else {
            max_predict_bias_step = fmaxf(max_predict_bias_step, bias_change);
        }
        memcpy(prior_bias, current_bias, sizeof(prior_bias));
        if (k >= 1000) {
            tail += body_attitude_error_deg(&s);
        }
    }
    float bias[3];
    body_get_gyro_bias(bias);
    float bias_err[3] = {bias[0] - s.gyro_bias[0], bias[1] - s.gyro_bias[1], bias[2] - s.gyro_bias[2]};
    check_true(body_is_pointing(), "sunlit: still pointing after 20 minutes");
    check_true(body_get_filter_resets() == 0, "sunlit: no filter resets");
    check_less("sunlit: mean attitude error over the last 200 s (deg)", tail / 200, 1.5);
    check_less("sunlit: gyro bias error (deg/s)", RAD2DEG(l2_norm(bias_err, 3)), 0.02);
    check_true(finite, "sunlit: currents always finite");
    check_true(bounded, "sunlit: currents always within +-Imax");
    check_true(any_nonzero, "sunlit: the controller actually commands something");
    check_true(quiet_window_ok, "ticks n through K command zero magnetic moment");
    check_less("predict-only ticks leave gyro bias unchanged", max_predict_bias_step, 1e-5);
    check_true(max_update_bias_step > max_predict_bias_step,
               "only measurement ticks make a meaningful gyro-bias correction");

    // Eclipse for 15 minutes: magnetometer-only filtering keeps us close
    float worst = 0;
    for (int k = 0; k < 900; k++) {
        body_sim_step(&s, false, 1, currents);
        worst = fmaxf(worst, body_attitude_error_deg(&s));
    }
    check_true(body_is_pointing(), "eclipse: still pointing after 15 minutes");
    check_true(body_get_filter_resets() == 0, "eclipse: no filter resets");
    check_less("eclipse: max attitude error (deg)", worst, 5.0);

    // A bad magnetometer forces detumbling on both prediction and measurement ticks.
    const char* bad_mag_cases[] = {"zero field", "NaN field", "NULL field",
                                   "NULL previous field", "NaN previous field",
                                   "zero field on measurement tick"};
    float zero_mag[3] = {0}, bad_prev[3] = {NAN, 0, 0};
    for (int bad = 0; bad < 6; bad++) {
        if (bad == 5) {
            for (unsigned int k = 0; k < BODY_MEASUREMENT_UPDATE_PERIOD_TICKS - 1u; k++) {
                body_sim_step(&s, true, 1, currents);
            }
        }
        body_sim_advance(&s, 1);
        float m1[3], g1[3], d1[NUM_DIODES];
        body_sim_sensors(&s, true, m1, g1, d1);
        const float* current = m1;
        const float* previous = s.last_mag;
        if (bad == 0 || bad == 5) current = zero_mag;
        if (bad == 1) current = nan_mag;
        if (bad == 2) current = NULLPTR;
        if (bad == 3) previous = NULLPTR;
        if (bad == 4) previous = bad_prev;
        body(previous, current, g1, d1, NULLPTR, s.unix_time, s.jd_int,
             (float)s.jd_frac, 1.0f, currents);
        report(!body_is_pointing() && !body_get_attitude(q_dummy),
               "bad magnetometer enters detumbling", "%s", bad_mag_cases[bad]);
        body(previous, current, g1, d1, NULLPTR, s.unix_time, s.jd_int,
             (float)s.jd_frac, 1.0f, currents);
        report(!body_is_pointing(), "bad magnetometer stays in detumbling",
               "%s", bad_mag_cases[bad]);
        memcpy(s.last_mag, m1, sizeof(m1));
        body_sim_step(&s, true, 1, currents);
        report(body_is_pointing(), "valid magnetometer can resume pointing",
               "%s", bad_mag_cases[bad]);
    }

    // A large rate while pointing is treated as a new tumble or gyro fault.
    {
        float m1[3], g1[3], d1[NUM_DIODES];
        body_sim_advance(&s, 1);
        body_sim_sensors(&s, true, m1, g1, d1);
        g1[0] += 3.0f; // 3 rad for one second
        body(s.last_mag, m1, g1, d1, NULLPTR, s.unix_time, s.jd_int, (float)s.jd_frac, 1.0f, currents);
        memcpy(s.last_mag, m1, sizeof(m1));
    }
    int steps_to_reset = -1;
    for (int k = 0; k < 20; k++) {
        body_sim_step(&s, true, 1, currents);
        if (body_get_filter_resets() > 0) {
            steps_to_reset = k;
            break;
        }
    }
    check_true(steps_to_reset >= 0, "gyro glitch: divergence detected and the filter is reset");
    for (int k = 0; k < 5; k++) {
        body_sim_step(&s, true, 1, currents);
    }
    check_true(body_is_pointing(), "gyro glitch: re-initialized with QUEST and pointing again");
    check_less("gyro glitch: attitude error after re-initialization (deg)", body_attitude_error_deg(&s), 5.0);

    // Position is propagated internally when no update is given
    body_reset();
    body_sim_init(&s);
    body_set_tle_metadata(s.tle_bstar, s.tle_epoch_jd);
    body(mag_prev, mag, gyro, diodes_dark, s.tle_mean, s.unix_time, s.jd_int, (float)s.jd_frac, 1.0f, currents);
    body(mag_prev, mag, gyro, diodes_dark, NULLPTR, s.unix_time + 1, s.jd_int, (float)s.jd_frac, 1.0f, currents);
    check_true(all_finite(currents, 3) && !body_is_pointing(), "propagating without a position update works");
}

// =============================================================================================
// MATLAB parity: replay the inputs recorded by tools/matlab/generate_matlab_reference.m through
// the C code and compare with what the team's MATLAB produced. Known MATLAB-side bugs are
// printed as [NOTE]s; the C code is then compared with the MATLAB after that single fix.
// =============================================================================================

static void to_float(const double* in, float* out, int n) {
    for (int i = 0; i < n; i++) {
        out[i] = (float)in[i];
    }
}

// |got - expected| / max(|expected|, floor), 2-norm
static double rel_diff(const float* got, const double* expected, int n, double floor) {
    double d = 0, e = 0;
    for (int i = 0; i < n; i++) {
        d += ((double)got[i] - expected[i]) * ((double)got[i] - expected[i]);
        e += expected[i] * expected[i];
    }
    double scale = sqrt(e) > floor ? sqrt(e) : floor;
    return isfinite(d) ? sqrt(d) / scale : INFINITY;
}

static double max_abs_diff(const float* got, const double* expected, int n) {
    double worst = 0;
    for (int i = 0; i < n; i++) {
        double d = fabs((double)got[i] - expected[i]);
        worst = (d > worst || !isfinite(d)) ? d : worst;
    }
    return worst;
}

static float quat_angle_deg_d(const float* q, const double* qd) {
    float qf[4];
    to_float(qd, qf, 4);
    return RAD2DEG(quat_angle_between(q, qf));
}

static float vec_angle_deg_d(const float* v, const double* vd) {
    float vf[3];
    to_float(vd, vf, 3);
    return RAD2DEG(vec_angle(v, vf));
}

static void parity_quaternion(void) {
    double w_mul = 0, w_norm = 0, w_inv = 0, w_apply = 0, w_q2rv = 0, w_rv2q = 0, w_diff = 0;
    for (int k = 0; k < ML_QUAT_IN_N; k++) {
        float q1[4], q2[4], v[3], rv[3], out[4], vec[3], q1n[4], q2n[4];
        to_float(&ML_QUAT_IN[k][0], q1, 4);
        to_float(&ML_QUAT_IN[k][4], q2, 4);
        to_float(&ML_QUAT_IN[k][8], v, 3);
        to_float(&ML_QUAT_IN[k][11], rv, 3);
        const double* e = ML_QUAT_OUT[k];

        quat_multiply(q1, q2, out);
        w_mul = fmax(w_mul, rel_diff(out, &e[0], 4, 1));
        quat_norm(q1, q1n);
        w_norm = fmax(w_norm, rel_diff(q1n, &e[4], 4, 1));
        quat_inv(q1, out);
        w_inv = fmax(w_inv, rel_diff(out, &e[8], 4, 1e-3));
        quat_apply(q1, v, vec);
        w_apply = fmax(w_apply, rel_diff(vec, &e[12], 3, 1));
        quat2rotationvec(q1n, vec);
        w_q2rv = fmax(w_q2rv, max_abs_diff(vec, &e[15], 3));
        rotationvec2quat(rv, out);
        w_rv2q = fmax(w_rv2q, max_abs_diff(out, &e[18], 4));
        quat_norm(q2, q2n);
        float d[4];
        quat_diff(q1n, q2n, d);
        quat2rotationvec(d, vec);
        w_diff = fmax(w_diff, max_abs_diff(vec, &e[22], 3));
    }
    check_less("quat_multiply == Quaternion.quaternion_multiply (relative)", w_mul, 1e-6);
    check_less("quat_norm == quaternion_normalize (relative)", w_norm, 1e-6);
    check_less("quat_inv == quaternion_inverse (relative)", w_inv, 1e-5);
    check_less("quat_apply == apply_rotation (relative)", w_apply, 1e-5);
    check_less("quat2rotationvec == quaternion2rotation_vec (rad)", w_q2rv, 1e-5);
    check_less("rotationvec2quat == rotation_vec2quaternion (incl. tiny and ~pi angles)", w_rv2q, 1e-6);
    check_less("quat_diff + quat2rotationvec == Quaternion.quat_diff (rad)", w_diff, 1e-5);
}

static void parity_control(void) {
    double worst = 0;
    for (int k = 0; k < ML_BDOT_IN_N; k++) {
        float Mt[3], Mt1[3], m[3];
        to_float(&ML_BDOT_IN[k][0], Mt, 3);
        to_float(&ML_BDOT_IN[k][3], Mt1, 3);
        Bdot(Mt, Mt1, (float)ML_BDOT_IN[k][6], (float)ML_BDOT_IN[k][7], m);
        worst = fmax(worst, rel_diff(m, ML_BDOT_OUT[k], 3, 1e-30));
    }
    check_less("Bdot == bDot.m (relative)", worst, 1e-4);

    worst = 0;
    int saturated = 0;
    for (int k = 0; k < ML_PD_IN_N; k++) {
        float omega[3], q[4], tau[3];
        for (int i = 0; i < 3; i++) {
            omega[i] = (float)(ML_PD_IN[k][i] * M_PI / 180.0); // MATLAB takes deg/s
        }
        to_float(&ML_PD_IN[k][3], q, 4);
        pd_loop(q, omega, tau);
        worst = fmax(worst, max_abs_diff(tau, ML_PD_OUT[k], 3) / (1e-3 + fabs(ML_PD_OUT[k][0]) + fabs(ML_PD_OUT[k][1]) + fabs(ML_PD_OUT[k][2])));
        for (int i = 0; i < 3; i++) {
            saturated += fabs(fabs(ML_PD_OUT[k][i]) - PD_MAX_TAU) < 1e-9;
        }
    }
    check_less("pd_loop == PD_loop.m (relative, incl. saturated cases)", worst, 1e-4);
    check_true(saturated > 0, "PD comparison includes saturated outputs");

    // torque2moment3axis.m solves B x m = tau; the C code solves m x B = tau (see torque2moments.c)
    double worst_neg = 0, worst_same = 0;
    int ml_wrong_way = 0;
    for (int k = 0; k < ML_T2M_IN_N; k++) {
        float tau[3], B[3], m[3];
        to_float(&ML_T2M_IN[k][0], tau, 3);
        to_float(&ML_T2M_IN[k][3], B, 3);
        torque_2_moments(B, tau, m);
        double neg[3] = {-ML_T2M_OUT[k][0], -ML_T2M_OUT[k][1], -ML_T2M_OUT[k][2]};
        worst_neg = fmax(worst_neg, rel_diff(m, neg, 3, 1e-30));
        worst_same = fmax(worst_same, rel_diff(m, ML_T2M_OUT[k], 3, 1e-30));
        float m_ml[3], produced[3];
        to_float(ML_T2M_OUT[k], m_ml, 3);
        cross(m_ml, B, produced);
        ml_wrong_way += dot3(produced, tau) < 0;
    }
    if (worst_same < 1e-4) {
        check_less("torque_2_moments == torque2moment3axis.m (relative)", worst_same, 1e-4);
    } else {
        check_less("torque_2_moments == -torque2moment3axis.m (relative; same magnitude, opposite sign)", worst_neg, 1e-4);
        note("torque2moment3axis.m: for %d/%d cases its dipole makes a torque (m x B) OPPOSITE to the "
             "request. pinv(skew(B)) solves B x m = tau; a magnetorquer makes m x B. The C code keeps the "
             "physical sign (verified in 'torque -> moment -> current').",
             ml_wrong_way, ML_T2M_IN_N);
    }

    if (strlen(ML_M2C_RAW_ERROR) > 0) {
        note("moment2current3axis.m errors as written: \"%s\". Compared against it with `mu_r = [1 1 1]` "
             "(air core, as in the C code) %s.",
             ML_M2C_RAW_ERROR, ML_M2C_FIX_APPLIED ? "added" : "(fix could not be applied!)");
    }
    // The .m still uses unit coil parameters (I = m). Scale m by k so the C side divides back to the
    // same unclamped current; this checks the clamp / no-limit behavior against MATLAB.
    worst = 0;
    for (int k = 0; k < ML_M2C_IN_N; k++) {
        float m[3], lim[3], I[3];
        to_float(&ML_M2C_IN[k][0], m, 3);
        to_float(&ML_M2C_IN[k][3], lim, 3);
        for (int i = 0; i < 3; i++) m[i] *= MTQ_K[i];
        moment2current3axis(m, lim, I);
        worst = fmax(worst, max_abs_diff(I, ML_M2C_OUT[k], 3));
    }
    check_less("moment2current3axis == moment2current3axis.m (A, incl. clamped and unlimited)", worst, 1e-6);
}

static void parity_orbits_time(void) {
    double wr = 0, wv = 0;
    for (int k = 0; k < ML_O2E_IN_N; k++) {
        float el[6], r[3], v[3];
        to_float(ML_O2E_IN[k], el, 6);
        orbital_to_eci_posvel(el, 398600.4418f, r, v);
        wr = fmax(wr, rel_diff(r, &ML_O2E_OUT[k][0], 3, 1));
        wv = fmax(wv, rel_diff(v, &ML_O2E_OUT[k][3], 3, 1e-3));
    }
    check_less("orbital_to_eci_posvel position == orbitalToECI.m (relative)", wr, 1e-5);
    check_less("orbital_to_eci_posvel velocity == orbitalToECI.m (relative)", wv, 1e-5);

    double wpos = 0, wg = 0;
    for (int k = 0; k < ML_E2E_IN_N; k++) {
        float r[3], out[3];
        to_float(ML_E2E_IN[k], r, 3);
        int t = (int)ML_E2E_IN[k][3];
        ecef_2_eci(r, out, t);
        wpos = fmax(wpos, max_abs_diff(out, ML_E2E_OUT[k], 3));
        wg = fmax(wg, fabs(remainder(jd_2_gmst_deg(unix_2_jd(t)) - ML_E2E_OUT[k][3], 360.0)));
    }
    check_less("ecef_2_eci == eceftoeci.m (m, |r| ~ 6.4e6)", wpos, 5.0);
    check_less("jd_2_gmst_deg == jd_to_gmst (deg)", wg, 1e-8);

    double ws = 0;
    for (int k = 0; k < ML_SUN_IN_N; k++) {
        float s[3];
        sun_vec((int)ML_SUN_IN[k][0], s);
        ws = fmax(ws, max_abs_diff(s, ML_SUN_OUT[k], 3));
    }
    check_less("sun_vec == sunVectorECI.m", ws, 1e-6);

    double wnu = 0, wpass = 0;
    for (int k = 0; k < ML_KEP_IN_N; k++) {
        const double* in = ML_KEP_IN[k];
        float out[6];
        propogateOrbitalElements((float)in[0], (float)in[1], (float)in[2], (float)in[3], (float)in[4], (float)in[5],
                                 (float)in[6], out);
        wnu = fmax(wnu, fabs(remainder(out[5] - ML_KEP_OUT[k][5], 360.0)));
        for (int i = 0; i < 5; i++) {
            wpass = fmax(wpass, fabs(out[i] - ML_KEP_OUT[k][i]) / (1 + fabs(ML_KEP_OUT[k][i])));
        }
    }
    check_less("propogateOrbitalElements true anomaly == Two Body Propogation.m (deg)", wnu, 0.02);
    check_less("propogateOrbitalElements passes a, e, i, RAAN, argp through like MATLAB", wpass, 1e-6);
}

static void parity_magnetosphere(void) {
    double w_year = 0, w_gmst = 0, w_lat = 0, w_lon = 0, w_alt = 0;
    double w_ned_fix = 0, w_ne_raw = 0, raw_bn_diff = 0;
    double w_b_fix_rel = 0, w_b_fix_ang = 0, raw_ang = 0, raw_rel = 0;
    for (int k = 0; k < ML_WMM_IN_N; k++) {
        const double* in = ML_WMM_IN[k];
        const double* e = ML_WMM_OUT[k];
        double jd = in[3] + in[4];
        w_year = fmax(w_year, fabs(jd2year(jd) - e[6]));
        w_gmst = fmax(w_gmst, fabs(remainder(gmst_from_jd(jd) - e[7], 2 * M_PI)));

        // Same ECI -> ECEF rotation MATLAB used, so geodetic conversion is compared on equal inputs
        double th = e[7];
        float r_ecef[3] = {(float)(cos(th) * in[0] + sin(th) * in[1]), (float)(-sin(th) * in[0] + cos(th) * in[1]),
                           (float)in[2]};
        float lat, lon, alt;
        ecef_to_geodetic(r_ecef, &lat, &lon, &alt);
        w_lat = fmax(w_lat, fabs(lat - e[8]));
        w_lon = fmax(w_lon, fabs(remainder(lon - e[9], 2 * M_PI)));
        w_alt = fmax(w_alt, fabs(alt - e[10]));

        // Field synthesis in NED on MATLAB's own geodetic coordinates
        static float g[WMM_DIM][WMM_DIM], h[WMM_DIM][WMM_DIM];
        load_wmm2025(jd, g, h);
        float ned[3];
        synthesize_mag_field((float)e[8], (float)e[9], (float)e[10], g, h, &ned[0], &ned[1], &ned[2]);
        w_ned_fix = fmax(w_ned_fix, max_abs_diff(ned, &e[14], 3));
        double ed_raw[2] = {e[12], e[13]};
        w_ne_raw = fmax(w_ne_raw, max_abs_diff(&ned[1], ed_raw, 2));
        raw_bn_diff = fmax(raw_bn_diff, fabs(ned[0] - e[11]));

        // Full ECI field (C returns Tesla, MATLAB nT)
        float r[3], B[3];
        to_float(in, r, 3);
        wmm_eci_embedded_v2(r, (int32_t)in[3], (float)in[4], B);
        for (int i = 0; i < 3; i++) {
            B[i] *= 1e9f;
        }
        w_b_fix_rel = fmax(w_b_fix_rel, rel_diff(B, &e[3], 3, 1));
        w_b_fix_ang = fmax(w_b_fix_ang, vec_angle_deg_d(B, &e[3]));
        raw_ang = fmax(raw_ang, vec_angle_deg_d(B, &e[0]));
        raw_rel = fmax(raw_rel, rel_diff(B, &e[0], 3, 1));
    }
    check_less("jd2year == magnetosphere.m jd2year (years)", w_year, 1e-4);
    check_less("gmst_from_jd == magnetosphere.m gmst_from_jd (rad)", w_gmst, 1e-6);
    check_less("ecef_to_geodetic latitude == ecef2geodetic (rad)", w_lat, 1e-6);
    check_less("ecef_to_geodetic longitude == ecef2geodetic (rad)", w_lon, 1e-6);
    check_less("ecef_to_geodetic altitude == ecef2geodetic (m)", w_alt, 2.0);
    check_less("synthesize_mag_field East/Down == synthesizeMagField_manual (nT)", w_ne_raw, 1.0);
    check_less("synthesize_mag_field N/E/D == synthesizeMagField_manual with the B_N sign fix (nT)", w_ned_fix, 1.0);
    check_less("wmm_eci_embedded_v2 == magnetosphere.m with both fixes (relative)", w_b_fix_rel, 1e-4);
    check_less("wmm_eci_embedded_v2 == magnetosphere.m with both fixes (direction, deg)", w_b_fix_ang, 0.01);
    if (!ML_WMM_FIXES_APPLIED) {
        note("the magnetosphere.m fixes no longer match its text (changed upstream?); 'fixed' == raw here");
    }
    note("magnetosphere.m synthesizeMagField_manual: B_N = -X'cos(psi) - Z'sin(psi) should be + Z'sin(psi); "
         "up to %.0f nT North error at mid-latitudes (East/Down agree with C).",
         raw_bn_diff);
    note("magnetosphere.m computeWMMfieldFromCoeffs_manual multiplies NED by R_ned2ecef, which is the "
         "ECEF->NED matrix (needs a transpose). Raw MATLAB field vs C: up to %.1f deg direction and %.0f%% "
         "magnitude difference. Both fixed MATLAB and C agree with an independent -grad(V) reference.",
         raw_ang, 100 * raw_rel);
}

static void parity_pointing(void) {
    if (strlen(ML_USED_SHIMS) > 0) {
        note("pointing_error.m ran with shims for missing toolbox functions: %s (tools/matlab/shims)", ML_USED_SHIMS);
    }
    double wq = 0, wz = 0;
    for (int k = 0; k < ML_POINT_IN_N; k++) {
        float r[3], v[3], qb[4], target[3], q_err[4], z[3], prev[4] = {1, 0, 0, 0};
        to_float(&ML_POINT_IN[k][0], r, 3);
        to_float(&ML_POINT_IN[k][3], v, 3);
        to_float(&ML_POINT_IN[k][6], qb, 4);
        to_float(&ML_POINT_IN[k][10], target, 3);
        pointing_error(r, v, qb, target, prev, q_err, z);
        wq = fmax(wq, quat_angle_deg_d(q_err, &ML_POINT_OUT[k][0]));
        wz = fmax(wz, max_abs_diff(z, &ML_POINT_OUT[k][4], 3));
    }
    check_less("pointing_error q_err == pointing_error.m (deg)", wq, 1e-3);
    check_less("pointing_error z_want == pointing_error.m", wz, 1e-6);
}

// C local-frame covariance -> MATLAB global frame: attitude rows/cols rotated by R(q)
static void local_to_global_cov(const float* P, const float* q, double* out) {
    float Rm[9];
    quat2rotm(q, Rm);
    double T[36] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            T[i * 6 + j] = Rm[i * 3 + j];
        }
        T[(i + 3) * 6 + (i + 3)] = 1;
    }
    double TP[36] = {0};
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 6; k++) {
                TP[i * 6 + j] += T[i * 6 + k] * P[k * 6 + j];
            }
        }
    }
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            double acc = 0;
            for (int k = 0; k < 6; k++) {
                acc += TP[i * 6 + k] * T[j * 6 + k];
            }
            out[i * 6 + j] = acc;
        }
    }
}

typedef struct {
    bool all_ok;
    double att_diff_deg, bias_diff, cov_rel_diff, cov_rel_diff_predict, max_correction_deg;
    double trace_ratio_min, trace_ratio_max, trace_ratio_mean; // C / MATLAB attitude variance, update steps
} ukf_step_stats_t;

static ukf_step_stats_t run_ukf_steps(const double (*in_tab)[27], const double (*out_tab)[46], int n, float p_att,
                                      float q_scale) {
    ukf_step_stats_t st = {true, 0, 0, 0, 0, 0, 1e30, 0, 0};
    int n_updates = 0;
    float Q[36], R2[36], R1[9], P0[36];
    float qd[6] = {4e-6f * q_scale, 4e-6f * q_scale, 4e-6f * q_scale, 1e-10f * q_scale, 1e-10f * q_scale, 1e-10f * q_scale};
    diag(Q, 6, qd);
    float r2d[6] = {2e-5f, 2e-5f, 2e-5f, 1e-4f, 1e-4f, 1e-4f}, r1d[3] = {2e-5f, 2e-5f, 2e-5f};
    diag(R2, 6, r2d);
    diag(R1, 3, r1d);
    float p0d[6] = {p_att, p_att, p_att, 1e-6f, 1e-6f, 1e-6f};
    diag(P0, 6, p0d);
    for (int k = 0; k < n; k++) {
        const double* in = in_tab[k];
        const double* e = out_tab[k];
        float q[4], x[6], gyro[3], ref[6], body[6];
        to_float(&in[0], q, 4);
        to_float(&in[4], x, 6);
        to_float(&in[10], gyro, 3);
        to_float(&in[13], ref, 6);
        to_float(&in[19], body, 6);
        int nv = (int)in[25];
        float q_new[4], x_new[6], P_new[36];
        st.all_ok &= iterate(x, q, P0, body, ref, nv, gyro, Q, nv == 2 ? R2 : R1, (float)in[26], x_new, q_new, P_new)
                     == UKF_OK;
        st.att_diff_deg = fmax(st.att_diff_deg, quat_angle_deg_d(q_new, &e[0]));
        st.bias_diff = fmax(st.bias_diff, max_abs_diff(&x_new[3], &e[4], 3));
        st.max_correction_deg = fmax(st.max_correction_deg, quat_angle_deg_d(q, &e[0]));
        // MATLAB's covariance is about its mean before the final correction: mean = rv2q(x_hat)^-1 * q_new
        float x_ml[3], dq_ml[4], dq_inv[4], q_ml[4], mean_ml[4];
        to_float(&e[43], x_ml, 3);
        to_float(&e[0], q_ml, 4);
        rotationvec2quat(x_ml, dq_ml);
        quat_conj(dq_ml, dq_inv);
        quat_multiply(dq_inv, q_ml, mean_ml);
        double Pg[36], dmax = 0, pmax = 0;
        local_to_global_cov(P_new, mean_ml, Pg);
        for (int i = 0; i < 36; i++) {
            dmax = fmax(dmax, fabs(Pg[i] - e[7 + i]));
            pmax = fmax(pmax, fabs(e[7 + i]));
        }
        if (nv == 0) {
            st.cov_rel_diff_predict = fmax(st.cov_rel_diff_predict, dmax / pmax);
        }
        st.cov_rel_diff = fmax(st.cov_rel_diff, dmax / pmax);
        double tr_c = Pg[0] + Pg[7] + Pg[14], tr_ml = e[7] + e[7 + 7] + e[7 + 14];
        if (nv > 0) {
            double ratio = tr_c / tr_ml;
            st.trace_ratio_min = fmin(st.trace_ratio_min, ratio);
            st.trace_ratio_max = fmax(st.trace_ratio_max, ratio);
            st.trace_ratio_mean += ratio;
            n_updates++;
        }
    }
    st.trace_ratio_mean /= n_updates > 0 ? n_updates : 1;
    return st;
}

static void parity_ukf(void) {
    // Building blocks
    float lam = calculate_lambda(STATE_SIZE, UKF_ALPHA, UKF_KAPPA);
    check_close("calculate_lambda == simulink.m calculate_lambda", lam, ML_UKF_LAMBDA, 1e-7);
    float wc[NUM_SIGMAS], wm[NUM_SIGMAS];
    get_weights(lam, STATE_SIZE, UKF_ALPHA, UKF_BETA, wc, wm);
    check_less("get_weights covariance weights == simulink.m", max_abs_diff(wc, ML_UKF_WEIGHTS[0], NUM_SIGMAS), 1e-7);
    check_less("get_weights mean weights == simulink.m", max_abs_diff(wm, ML_UKF_WEIGHTS[1], NUM_SIGMAS), 1e-7);
    float x[6], P[36], sig[NUM_SIGMAS * STATE_SIZE];
    to_float(ML_UKF_SIGMA_X[0], x, 6);
    for (int i = 0; i < 6; i++) {
        to_float(ML_UKF_SIGMA_P[i], &P[i * 6], 6);
    }
    get_sigma_points(lam, x, P, sig);
    double wsig = 0;
    for (int i = 0; i < NUM_SIGMAS; i++) {
        wsig = fmax(wsig, max_abs_diff(&sig[i * 6], ML_UKF_SIGMAS[i], 6));
    }
    check_less("get_sigma_points == simulink.m get_sigma_points (same order, non-diagonal P)", wsig, 1e-5);

    // Single steps from identical states. The C filter uses body-frame (local) attitude errors and
    // adds Q before generating sigma points; simulink.m uses ECI-frame (global) errors and adds Q
    // after propagation, so its measurement update never sees Q (its posterior keeps all of Q).
    //  - Q = 0, small prior: the two are the same algorithm to second order -> tight agreement.
    //  - Q != 0 (the simulink.m tuning): estimates stay close and C's posterior covariance is
    //    smaller, exactly as that difference predicts.
    ukf_step_stats_t q0 = run_ukf_steps(ML_UKF_STEPQ0_IN, ML_UKF_STEPQ0_OUT, ML_UKF_STEPQ0_IN_N, 1e-4f, 0.0f);
    check_true(q0.all_ok, "single steps, Q = 0: every C step returns UKF_OK");
    check_less("single steps, Q = 0: updated attitude C vs simulink.m (deg)", q0.att_diff_deg, 0.002);
    check_less("single steps, Q = 0: updated gyro bias C vs simulink.m (rad/s)", q0.bias_diff, 1e-6);
    check_less("single steps, Q = 0: covariance C (rotated to ECI frame) vs simulink.m (relative)", q0.cov_rel_diff, 0.01);

    ukf_step_stats_t qs = run_ukf_steps(ML_UKF_STEP_IN, ML_UKF_STEP_OUT, ML_UKF_STEP_IN_N, 0.004f, 1.0f);
    check_true(qs.all_ok, "single steps, simulink.m tuning: every C step returns UKF_OK");
    report_value(true, "single steps, simulink.m tuning: (context) largest attitude change in one MATLAB step", "%.3f deg",
                 qs.max_correction_deg);
    check_less("single steps, simulink.m tuning: updated attitude C vs simulink.m (deg)", qs.att_diff_deg, 0.1);
    check_less("single steps, simulink.m tuning: updated gyro bias C vs simulink.m (rad/s)", qs.bias_diff, 2e-5);
    check_less("single steps, simulink.m tuning: predict-only covariance C vs simulink.m (relative)", qs.cov_rel_diff_predict, 0.01);
    check_close("single steps, Q = 0: posterior attitude variance ratio C / simulink.m", q0.trace_ratio_mean, 1.0, 1e-3);
    // MATLAB's update ignores Q, so it keeps all of Q in the posterior: C should be smaller on
    // average. Nonlinearity at ~3 deg corrections is the same order as Q here, so individual
    // steps can go either way by a little.
    check_less("single steps, simulink.m tuning: mean posterior attitude variance ratio C / simulink.m (< 1)", qs.trace_ratio_mean, 1.0);
    check_less("single steps, simulink.m tuning: max posterior attitude variance ratio C / simulink.m", qs.trace_ratio_max, 1.02);
    check_true(qs.trace_ratio_min > 0.8, "single steps, simulink.m tuning: min posterior attitude variance ratio C / simulink.m > 0.8");

    // Full trajectory: 2 min at 10 Hz, updates every 5th step, magnetometer-only for the second half
    const double* setup = ML_TRAJ_SETUP[0];
    float q_est[4], xs[6] = {0}, Ps[36], Qt[36], R2t[36], R1t[9], sun[3];
    to_float(&setup[4], q_est, 4);
    float d6[6];
    to_float(&setup[8], d6, 6);
    diag(Ps, 6, d6);
    to_float(&setup[14], d6, 6);
    diag(Qt, 6, d6);
    to_float(&setup[20], d6, 6);
    diag(R2t, 6, d6);
    float d3[3];
    to_float(&setup[26], d3, 3);
    diag(R1t, 3, d3);
    to_float(&setup[29], sun, 3);
    double w_vs_ml_sun = 0, w_vs_ml_ecl = 0, w_bias = 0;
    double c_sum[2] = {0, 0}, ml_sum[2] = {0, 0};
    int n_sum[2] = {0, 0};
    bool all_ok = true;
    int check_idx = 0;
    printf("  step | C err (deg) | MATLAB err (deg) | C vs MATLAB (deg)\n");
    for (int k = 0; k < ML_TRAJ_STEPS; k++) {
        const double* s = ML_TRAJ_STEP[k];
        float gyro[3], body[6], ref[6];
        to_float(&s[0], gyro, 3);
        to_float(&s[3], body, 6);
        to_float(&s[9], ref, 3);
        memcpy(&ref[3], sun, sizeof(sun));
        bool update = k % ML_TRAJ_UPDATE_EVERY == 0;
        int nv = !update ? 0 : (k < ML_TRAJ_SWITCH_STEP ? 2 : 1);
        all_ok &= iterate(xs, q_est, Ps, body, ref, nv, gyro, Qt, nv == 2 ? R2t : R1t, (float)ML_TRAJ_DT, xs, q_est, Ps)
                  == UKF_OK;
        if ((k + 1) % 10 == 0 && check_idx < ML_TRAJ_CHECK_N) {
            const double* c = ML_TRAJ_CHECK[check_idx++];
            float c_vs_ml = quat_angle_deg_d(q_est, &c[1]);
            float c_err = quat_angle_deg_d(q_est, &c[8]);
            float qml[4];
            to_float(&c[1], qml, 4);
            float ml_err = quat_angle_deg_d(qml, &c[8]);
            if (k >= 100) { // skip the initial convergence transient
                if (k < ML_TRAJ_SWITCH_STEP) {
                    w_vs_ml_sun = fmax(w_vs_ml_sun, c_vs_ml);
                } else {
                    w_vs_ml_ecl = fmax(w_vs_ml_ecl, c_vs_ml);
                }
                w_bias = fmax(w_bias, max_abs_diff(&xs[3], &c[5], 3));
            }
            if ((k + 1) % 200 == 0) {
                printf("  %4d | %11.4f | %16.4f | %.4f\n", k + 1, c_err, ml_err, c_vs_ml);
            }
            if (k >= 100) {
                int phase = k < ML_TRAJ_SWITCH_STEP ? 0 : 1;
                c_sum[phase] += c_err;
                ml_sum[phase] += ml_err;
                n_sum[phase]++;
            }
        }
    }
    check_true(all_ok, "trajectory: every C step returns UKF_OK");
    check_less("trajectory, sun + mag: max C vs simulink.m attitude difference after 10 s (deg)", w_vs_ml_sun, 0.05);
    // Magnetometer only: rotation about the field is only weakly observable, so the two filters'
    // errors wander independently along it (both stay ~0.5 deg from the truth, checked below)
    check_less("trajectory, mag only: max C vs simulink.m attitude difference (deg)", w_vs_ml_ecl, 1.0);
    check_less("trajectory: max C vs simulink.m gyro bias difference (rad/s)", w_bias, 2e-4);
    const char* phase_name[2] = {"sun + mag", "mag only"};
    for (int ph = 0; ph < 2; ph++) {
        double c_mean = c_sum[ph] / n_sum[ph], ml_mean = ml_sum[ph] / n_sum[ph];
        char name[160];
        snprintf(name, sizeof(name), "trajectory, %s: C mean error minus MATLAB mean error (deg; C %.3f, MATLAB %.3f)",
                 phase_name[ph], c_mean, ml_mean);
        check_less(name, c_mean - ml_mean, 0.1);
    }
}

void test_matlab_parity(void) {
    begin_suite("MATLAB parity");
    parity_quaternion();
    parity_control();
    parity_orbits_time();
    parity_magnetosphere();
    parity_pointing();
    parity_ukf();
}
