#include "include/closest_approach.h"
#include "include/ecef2eci.h"
#include "include/sgp4_step.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N_COARSE_SAMPLES ((int)(CLOSEST_APPROACH_WINDOW_S / CLOSEST_APPROACH_DT_COARSE_S) + 1)

/*
 * Port of rangeAt. Range (m) and unit vector (ECI) from the satellite to target_ecef at t seconds
 * past epoch. Returns false if SGP4 fails or the range is degenerate.
 */
static bool range_at(double t, const float* oe_epoch, float bstar, double epoch_jd,
                     const float* target_ecef, float* rho, float* u) {
    float r_km[3], v_km_s[3], oe_osc[7];
    if (sgp4_step(oe_epoch, bstar, epoch_jd, (float)t, r_km, v_km_s, oe_osc) != 0) {
        return false;
    }

    // Rotate the target into ECI by GMST at that instant (same rotation as ecef_2_eci).
    // jd_2_gmst_deg (src/ecef2eci.c) is the IAU GMST formula, same as gmst() in
    // closest_approach.m: Julian date -> Earth rotation angle in degrees, [0, 360). Kept in
    // double because a float JD only resolves ~0.25 days.
    double theta = jd_2_gmst_deg(epoch_jd + t / 86400.0) * (M_PI / 180.0);
    double c = cos(theta), s = sin(theta);
    double target_eci[3] = {c * target_ecef[0] - s * target_ecef[1],
                            s * target_ecef[0] + c * target_ecef[1], target_ecef[2]};

    double d[3];
    for (int i = 0; i < 3; i++) {
        d[i] = target_eci[i] - (double)r_km[i] * 1000.0;
    }
    double norm = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (!isfinite(norm) || !(norm > 0.0)) {
        return false;
    }
    *rho = (float)norm;
    if (u != NULL) {
        for (int i = 0; i < 3; i++) {
            u[i] = (float)(d[i] / norm);
        }
    }
    return true;
}

bool closest_approach(const float* oe_epoch, float bstar, double epoch_jd, double t0,
                      const float* target_ecef, float* u, double* t_star, float* rho_min) {
    memset(u, 0, sizeof(float) * 3);
    *t_star = t0;
    *rho_min = INFINITY;
    if (!isfinite(t0)) {
        return false;
    }

    // Coarse pass: running minimum of the range every DT_COARSE seconds
    float best_rho = INFINITY;
    double best_t = t0;
    bool found = false;
    for (int k = 0; k < N_COARSE_SAMPLES; k++) {
        double t = t0 + k * CLOSEST_APPROACH_DT_COARSE_S;
        float rho;
        if (range_at(t, oe_epoch, bstar, epoch_jd, target_ecef, &rho, NULL) && rho < best_rho) {
            best_rho = rho;
            best_t = t;
            found = true;
        }
    }
    if (!found) {
        return false;
    }

    // Refine at 1 s around the coarse minimum, not before the window start
    double ta = fmax(best_t - CLOSEST_APPROACH_DT_COARSE_S, t0);
    double tb = best_t + CLOSEST_APPROACH_DT_COARSE_S;
    for (int j = 0; ta + j <= tb; j++) {
        double t = ta + j;
        float rho;
        if (range_at(t, oe_epoch, bstar, epoch_jd, target_ecef, &rho, NULL) && rho < best_rho) {
            best_rho = rho;
            best_t = t;
        }
    }

    // Line of sight at the refined time
    float rho;
    if (!range_at(best_t, oe_epoch, bstar, epoch_jd, target_ecef, &rho, u)) {
        memset(u, 0, sizeof(float) * 3);
        return false;
    }
    *t_star = best_t;
    *rho_min = best_rho;
    return true;
}
