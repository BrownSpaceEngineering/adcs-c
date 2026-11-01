#ifndef CLOSEST_APPROACH_H
#define CLOSEST_APPROACH_H
#include <stdbool.h>

#define CLOSEST_APPROACH_DT_COARSE_S 10.0 // coarse search step (s)
#define CLOSEST_APPROACH_WINDOW_S 86400.0 // search window length (s)

/*
 * Port of closest_approach.m. Searches [t0, t0 + CLOSEST_APPROACH_WINDOW_S] for the time the
 * satellite is closest to target_ecef: a coarse pass every CLOSEST_APPROACH_DT_COARSE_S seconds,
 * then a 1 s pass around the coarse minimum.
 *
 * Unlike the MATLAB, this has no persistent cache: the caller runs it once per uplinked TLE and
 * keeps the result. The elements are SGP4 mean elements in body()'s posn_update order, so the
 * MATLAB's (a, nu) -> (n, M) conversion isn't needed.
 *
 * Runs ~8660 SGP4 propagations; don't call it every tick.
 *
 *   oe_epoch, bstar, epoch_jd: as for sgp4_step
 *   t0:          window start, s since epoch (normally (current UTC - epoch_jd) * 86400)
 *   target_ecef: target position, Earth-fixed (m)
 *   u:           unit vector from the satellite to the target at closest approach, ECI
 *   t_star:      time of closest approach, s since epoch
 *   rho_min:     range at closest approach (m)
 * Returns false (u zeroed) if t0 isn't finite or no SGP4 sample succeeded.
 */
bool closest_approach(const float* oe_epoch, float bstar, double epoch_jd, double t0,
                      const float* target_ecef, float* u, double* t_star, float* rho_min);
#endif
