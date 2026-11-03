#ifndef DOWN_QUAT_H
#define DOWN_QUAT_H
#include <stdbool.h>

// Port of pointing_error.m. Desired attitude: body +Z toward target, body Y = normalize(Z x v),
// X = Y x Z. Writes the body-frame error quaternion q_err = q_b2eci^-1 * q_want (feed this to
// pd_loop) and, if non-NULL, z_want (unit vector to target, ECI).
// q_want_prev (in/out, may be NULL) replaces MATLAB's persistent variable: it keeps q_want in
// the same hemisphere between calls. Initialize it to {1, 0, 0, 0}.
// Returns false (q_err = identity) if the geometry is degenerate.
bool pointing_error(const float* r_eci, const float* v_eci, const float* q_b2eci,
                    const float* target_eci, float* q_want_prev, float* q_err);
#endif
