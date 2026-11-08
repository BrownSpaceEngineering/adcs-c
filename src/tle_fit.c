/* Ground-side SGP4 mean-element fit (NOT flight code). Copied from
 * FINAL_SIMULINK/CubesatModel/C files/adcs_tle_fit.c for plant.c; keep the two in sync. */
#include "include/tle_fit.h"
#include "include/sgp4_init.h"
#include "include/sgp4_propagate.h"
#include <math.h>

#define MU_KM 398600.8 /* WGS-72, same as GRAV_MODEL 72 in sgp4_init.c */
#define TWO_PI 6.28318530717958647692

/* Nonsingular elements: [a (km), ex = e cos w, ey = e sin w, i, RAAN, lambda = w + M] */
static void rv2ns(const double* r, const double* v, double* x) {
    double h[3] = {r[1] * v[2] - r[2] * v[1], r[2] * v[0] - r[0] * v[2], r[0] * v[1] - r[1] * v[0]};
    double rn = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    double hn = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
    double rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
    double a = 1.0 / (2.0 / rn - v2 / MU_KM);
    double inc = acos(h[2] / hn);
    double raan = atan2(h[0], -h[1]);
    double ev[3];
    for (int k = 0; k < 3; k++) ev[k] = ((v2 - MU_KM / rn) * r[k] - rv * v[k]) / MU_KM;
    /* node line and the in-plane axis 90 deg ahead of it */
    double nx = cos(raan), ny = sin(raan);
    double px = -sin(raan) * cos(inc), py = cos(raan) * cos(inc), pz = sin(inc);
    double ex = ev[0] * nx + ev[1] * ny;
    double ey = ev[0] * px + ev[1] * py + ev[2] * pz;
    double u = atan2(r[0] * px + r[1] * py + r[2] * pz, r[0] * nx + r[1] * ny); /* arg of latitude */
    double e = sqrt(ex * ex + ey * ey);
    double w = atan2(ey, ex);
    double nu = u - w;
    double E = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(nu / 2.0));
    double M = E - e * sin(E);
    x[0] = a; x[1] = ex; x[2] = ey; x[3] = inc; x[4] = raan; x[5] = w + M;
}

static double wrap(double a) { return a - TWO_PI * floor((a + TWO_PI / 2.0) / TWO_PI); }

static void ns2tle(const double* x, float* tle) {
    double e = sqrt(x[1] * x[1] + x[2] * x[2]);
    double w = atan2(x[2], x[1]);
    tle[0] = (float)(sqrt(MU_KM / (x[0] * x[0] * x[0])) * 86400.0 / TWO_PI);
    tle[1] = (float)e;
    tle[2] = (float)x[3];
    tle[3] = (float)fmod(x[4] + 2.0 * TWO_PI, TWO_PI);
    tle[4] = (float)fmod(w + 2.0 * TWO_PI, TWO_PI);
    tle[5] = (float)fmod(x[5] - w + 4.0 * TWO_PI, TWO_PI);
}

double adcs_tle_fit(const double* r_gcrf_m, const double* v_gcrf_mps, double epoch_jd,
                    float bstar, float* tle_out) {
    sgp4_sat_t sat;
    double r_t[3], v_t[3], target[6], x[6];
    float tle[6], rg[3], vg[3], osc[7], mean[6], rt[3], vt[3];

    /* first pass: osculating GCRF elements, only to get the TEME rotation at this epoch */
    double rk[3], vk[3];
    for (int k = 0; k < 3; k++) { rk[k] = r_gcrf_m[k] / 1000.0; vk[k] = v_gcrf_mps[k] / 1000.0; }
    rv2ns(rk, vk, x);
    ns2tle(x, tle);
    sgp4_init(tle, bstar, epoch_jd, &sat);
    if (sat.initErr != 0) return -1.0;
    /* target state in TEME: R_teme2gcrf is row-major, so TEME = R^T * GCRF */
    for (int i = 0; i < 3; i++) {
        r_t[i] = 0.0; v_t[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            r_t[i] += sat.R_teme2gcrf[3 * j + i] * rk[j];
            v_t[i] += sat.R_teme2gcrf[3 * j + i] * vk[j];
        }
    }
    rv2ns(r_t, v_t, target);
    for (int k = 0; k < 6; k++) x[k] = target[k];

    /* fixed point: mean += (osculating target - osculating SGP4 output at dt = 0) */
    double resid = -1.0;
    for (int it = 0; it < 40; it++) {
        ns2tle(x, tle);
        sgp4_init(tle, bstar, epoch_jd, &sat);
        if (sat.initErr != 0) return -1.0;
        if (sgp4_propagate(&sat, 0.0f, rg, vg, osc, mean, rt, vt) != 0) return -1.0;
        double rr[3] = {rt[0], rt[1], rt[2]}, vv[3] = {vt[0], vt[1], vt[2]}, xo[6];
        rv2ns(rr, vv, xo);
        double d[6];
        for (int k = 0; k < 6; k++) d[k] = target[k] - xo[k];
        d[4] = wrap(d[4]); d[5] = wrap(d[5]);
        for (int k = 0; k < 6; k++) x[k] += d[k];
        double dr = 0.0;
        for (int k = 0; k < 3; k++) dr += (rr[k] - r_t[k]) * (rr[k] - r_t[k]);
        resid = sqrt(dr) * 1000.0;
        if (resid < 1.0) break;
    }
    ns2tle(x, tle_out);
    return resid;
}
