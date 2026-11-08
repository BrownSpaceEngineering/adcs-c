#include "include/newpointing.h"
#include "include/laextension.h"
#include "include/quat.h"
#include <math.h>
#include <string.h>

const float NEW_PD_MAX_DIPOLE[3] = {0.32f, 0.32f, 0.11f};

// out = R^T * v, with R row-major 3x3 (ECI->body when R is body->ECI)
static void mat3_transpose_mul(const float* R, const float* v, float* out) {
    for (int i = 0; i < 3; i++) {
        out[i] = R[0 * 3 + i] * v[0] + R[1 * 3 + i] * v[1] + R[2 * 3 + i] * v[2];
    }
}

static void zero_outputs(float* m_out) {
    if (m_out != NULL) memset(m_out, 0, sizeof(float) * 3);
}

bool new_pd(const float* r_eci, const float* v_eci, const float* q_b2eci, const float* b_eci,
            const float* w_eci2body_in_body, float* m_out) {
    zero_outputs(m_out);
    if (r_eci == NULL || v_eci == NULL || q_b2eci == NULL || b_eci == NULL
        || w_eci2body_in_body == NULL || !all_finite(r_eci, 3) || !all_finite(v_eci, 3)
        || !all_finite(q_b2eci, 4) || !all_finite(b_eci, 3) || !all_finite(w_eci2body_in_body, 3)) {
        return false;
    }

    // nadir_eci = -r_eci / norm(r_eci)
    float nadir_eci[3] = {-r_eci[0], -r_eci[1], -r_eci[2]};
    if (!normalize_vec(nadir_eci, 3)) return false;

    // R_body_to_eci = quat2rotm(q_b2eci); R_eci_to_body = R_body_to_eci.'
    // (quat2rotm normalizes q first, like MATLAB's)
    float R_body_to_eci[9];
    quat2rotm(q_b2eci, R_body_to_eci);

    // w_eci2orbital_in_eci = cross(r_eci, v_eci) / dot(r_eci, r_eci)
    float r_dot_r = dot3(r_eci, r_eci);
    if (!(r_dot_r > 0.0f) || !isfinite(r_dot_r)) return false;
    float w_eci2orbital_in_eci[3];
    cross(r_eci, v_eci, w_eci2orbital_in_eci);
    for (int i = 0; i < 3; i++) w_eci2orbital_in_eci[i] /= r_dot_r;

    // Rotate into the body frame
    float w_eci2orbital_in_body[3], nadir_body[3], b_body[3];
    mat3_transpose_mul(R_body_to_eci, w_eci2orbital_in_eci, w_eci2orbital_in_body);
    mat3_transpose_mul(R_body_to_eci, nadir_eci, nadir_body);
    mat3_transpose_mul(R_body_to_eci, b_eci, b_body);

    // w_orbital2body_in_body = w_eci2body_in_body - w_eci2orbital_in_body
    float w_orbital2body_in_body[3];
    for (int i = 0; i < 3; i++) {
        w_orbital2body_in_body[i] = w_eci2body_in_body[i] - w_eci2orbital_in_body[i];
    }

    // torque_cmd = kp * cross(z_current_body, nadir_body) - kd * w_orbital2body_in_body
    const float z_current_body[3] = {0.0f, 0.0f, 1.0f};
    float z_cross_nadir[3];
    cross(z_current_body, nadir_body, z_cross_nadir);
    float torque_cmd[3];
    for (int i = 0; i < 3; i++) {
        torque_cmd[i] = NEW_PD_KP * z_cross_nadir[i] - NEW_PD_KD * w_orbital2body_in_body[i];
    }

    // m = cross(b_body, torque_cmd); m = max(-maxdipole, min(m, maxdipole))
    float m[3];
    cross(b_body, torque_cmd, m);
    for (int i = 0; i < 3; i++) {
        m[i] = fmaxf(-NEW_PD_MAX_DIPOLE[i], fminf(m[i], NEW_PD_MAX_DIPOLE[i]));
    }

    if (m_out != NULL) memcpy(m_out, m, sizeof(m));
    return true;
}
