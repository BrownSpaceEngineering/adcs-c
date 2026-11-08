#ifndef TLE_FIT_H
#define TLE_FIT_H
/*
 * Ground-side helper for the Simulink harness (NOT flight code).
 * Fits SGP4 mean elements so that the adcs-unix-testing SGP4 (sgp4_init/sgp4_propagate)
 * reproduces a given GCRF state vector at the epoch -- i.e. "make the TLE we would uplink".
 *
 *   r_gcrf_m, v_gcrf_mps : state at epoch, GCRF/ICRF, m and m/s
 *   epoch_jd             : UTC Julian date of the state
 *   bstar                : B* to fit with (0 when drag is off)
 *   tle_out              : [n (rev/day), e, i, RAAN, argp, M (rad)], body() posn_update order
 * Returns the remaining position residual at epoch in metres (negative on failure).
 */
double adcs_tle_fit(const double* r_gcrf_m, const double* v_gcrf_mps, double epoch_jd,
                    float bstar, float* tle_out);
#endif
