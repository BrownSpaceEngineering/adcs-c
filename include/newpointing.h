#ifndef NEWPOINTING_H
#define NEWPOINTING_H
#include <stdbool.h>

// Gains / limits from the MATLAB nadir-pointing fcn
#define NEW_PD_KP 2e4f
#define NEW_PD_KD 1.5e5f
extern const float NEW_PD_MAX_DIPOLE[3]; // A*m^2, per body axis

/**
 * Nadir-pointing PD (port of the MATLAB fcn): drives body +Z toward nadir and damps the body rate
 * relative to the orbital frame, then maps the torque to a saturated magnetic dipole.
 *
 *   r_eci, v_eci:        position (m) and velocity (m/s), ECI
 *   q_b2eci:             body->ECI quaternion, ACTIVE, scalar first
 *   b_eci:               magnetic field (T), ECI
 *   w_eci2body_in_body:  body rate relative to ECI, body coordinates, rad/s
 *
 * Outputs (any may be NULL):
 *   m_out:      commanded dipole (A*m^2), clamped to +-NEW_PD_MAX_DIPOLE
 *   torque_out: torque that dipole produces, m x b_body
 *   crs_out:    |z_body x nadir_body| = sin(pointing error)
 *
 * Returns false (outputs zeroed) on non-finite inputs or a zero position vector.
 */
bool new_pd(const float* r_eci, const float* v_eci, const float* q_b2eci, const float* b_eci,
            const float* w_eci2body_in_body, float* m_out);
#endif
