#include "include/laextension.h"
#include "include/quat.h"
#include <string.h>

#define EPS 1e-6f

static const float IDENTITY_Q[4] = {1.0f, 0.0f, 0.0f, 0.0f};

// Row-major rotation matrix with the given columns (body axes expressed in ECI)
static void rotm_from_columns(const float* x, const float* y, const float* z, float* R) {
    R[0] = x[0];
    R[1] = y[0];
    R[2] = z[0];
    R[3] = x[1];
    R[4] = y[1];
    R[5] = z[1];
    R[6] = x[2];
    R[7] = y[2];
    R[8] = z[2];
}
/**
 * \fn pointing_error
 *
 * \brief Port of pointing_error.m (see header).
 */
bool pointing_error(const float* r_eci, const float* v_eci, const float* q_b2eci,
                    const float* target_eci, float* q_want_prev, float* q_err) {
    float z[3] = {target_eci[0] - r_eci[0], target_eci[1] - r_eci[1], target_eci[2] - r_eci[2]};
    if (!normalize_vec(z, 3)) {
        memcpy(q_err, IDENTITY_Q, sizeof(IDENTITY_Q));
        return false;
    }
  
    float y[3];
    cross(z, v_eci, y);
    if (!(l2_norm(y, 3) > EPS * l2_norm(v_eci, 3)) || !normalize_vec(y, 3)) {
        // Velocity (nearly) along the line of sight; MATLAB would divide by zero here
        memcpy(q_err, IDENTITY_Q, sizeof(IDENTITY_Q));
        return false;
    }

    float x[3];
    cross(y, z, x);
    normalize_vec(x, 3);

    float R_want[9];
    rotm_from_columns(x, y, z, R_want);
    float q_want[4];
    rotm_to_quat(R_want, q_want);

    if (q_want_prev != NULL) {
        float d = q_want[0] * q_want_prev[0] + q_want[1] * q_want_prev[1]
                  + q_want[2] * q_want_prev[2] + q_want[3] * q_want_prev[3];
        if (d < 0.0f) {
            for (int i = 0; i < 4; i++) {
                q_want[i] = -q_want[i];
            }
        }
        memcpy(q_want_prev, q_want, sizeof(q_want));
    }

    float q_inv[4];
    quat_inv(q_b2eci, q_inv);
    quat_multiply(q_inv, q_want, q_err);
    return true;
}
