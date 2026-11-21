#ifndef NEWPOINTING_H
#define NEWPOINTING_H
#include <stdbool.h>

// Gains / limits. The MATLAB nadir-pointing fcn used kp = 2e4, kd = 1.5e5. Chosen from a host Monte
// Carlo over the full uncertainty envelope (tools/mc: inertia +-20% / 10 deg, orbit, tip-off, sensor,
// coil errors; 200 runs each, 6 orbits), settled nadir p95 median / P90 and runs with p95 <= 10 deg:
//   kp 2e4 kd 1.5e5: 14.7 / 38.2 deg, 46/187      kp 4e4 kd 3e5: 9.8 / 28.7 deg, 96/187
//   kp 8e4 kd 6e5  : 11.1 / 59.5 deg, 87/187 (best with only sensor noise varied, but an instability
//   tail for some inertia / orbit combinations).
#ifndef NEW_PD_KP
#define NEW_PD_KP 4e4f
#endif
#ifndef NEW_PD_KD
#define NEW_PD_KD 3e5f
#endif
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
