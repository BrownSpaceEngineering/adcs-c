/*
 * Host closed loop: plant.c truth + sensors -> body() -> coil currents -> plant.c.
 * This is the PIL loop with the target replaced by a direct body() call.
 *
 *   ./build/adcs-plant-sim [orbits=3] [worst_case=0] [seed=1] [w0_degps=0] [csv=plant_sim.csv]
 *
 * Initial conditions follow the Simulink model (asbCubeSat): 2025-01-01 00:00 UTC, a = 7021 km,
 * e = 1e-5, i = 60 deg, RAAN = 137 deg, argp = 349.1886 deg, nu = 0, attitude [.7071 0 .7071 0],
 * rate = w0 * [1 -1 1]/sqrt(3) (Simulink uses adcs_worst_case*[10 -10 10] deg/s).
 */
#include "include/body.h"
#include "include/plant.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEG (M_PI / 180.0)

static double angle_between(const double* a, const double* b) {
    double c = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
               / (sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) * sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]));
    return acos(c > 1 ? 1 : (c < -1 ? -1 : c));
}

/* Per-tick logging / scoring, driven by plant_run()'s on_tick callback */
typedef struct {
    FILE* csv;
    double t_orbit;
    long n_ticks;
    double orbit_sum, orbit_max;
    long orbit_n, uplinks, pointing_ticks;
    int orbit_idx;
} logger_t;

static void log_tick(long k, const plant_body_inputs_t* in, const float* I, const plant_truth_t* tr, void* ctx) {
    logger_t* L = (logger_t*)ctx;
    if (in->has_posn_update) L->uplinks++;
    /* truth nadir error: body +z vs -r */
    const double* q = tr->q_body_to_eci;
    double zb[3] = {2 * (q[1] * q[3] + q[0] * q[2]), 2 * (q[2] * q[3] - q[0] * q[1]),
                    q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]};
    double nadir[3] = {-tr->r_eci_m[0], -tr->r_eci_m[1], -tr->r_eci_m[2]};
    double nadir_err = angle_between(zb, nadir) / DEG;
    float qe[4];
    double est_err = NAN;
    if (body_get_attitude(qe)) {
        double d = fabs(qe[0] * q[0] + qe[1] * q[1] + qe[2] * q[2] + qe[3] * q[3])
                   / sqrt(qe[0] * qe[0] + qe[1] * qe[1] + qe[2] * qe[2] + qe[3] * qe[3]);
        est_err = 2.0 * acos(d > 1 ? 1 : d) / DEG;
    }
    float rp[3];
    double pos_err = NAN;
    if (body_get_position(rp)) {
        double dr[3] = {rp[0] - tr->r_eci_m[0], rp[1] - tr->r_eci_m[1], rp[2] - tr->r_eci_m[2]};
        pos_err = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]) / 1000.0;
    }
    const double* w = tr->w_body_rad_s;
    double wn = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]) / DEG;
    bool pt = body_is_pointing();
    if (pt) {
        L->pointing_ticks++;
        L->orbit_sum += nadir_err;
        L->orbit_max = fmax(L->orbit_max, nadir_err);
        L->orbit_n++;
    }
    if (L->csv && k % 100 == 0) {
        fprintf(L->csv, "%.1f,%d,%d,%.4f,%.4f,%.5f,%.4f,%d,%.5f,%.5f,%.5f\n", tr->t_s, pt, body_get_filter_resets(),
                nadir_err, est_err, wn, pos_err, tr->sunlit, I[0], I[1], I[2]);
    }
    if (tr->t_s >= L->orbit_idx * L->t_orbit || k == L->n_ticks - 1) {
        printf("%5d | %8.1f | %15.2f / %-22.2f | %13.3f | %11.4f | %12.3f | %d\n", L->orbit_idx,
               100.0 * (double)L->orbit_n / (L->t_orbit / PLANT_DT), L->orbit_n ? L->orbit_sum / (double)L->orbit_n : NAN,
               L->orbit_n ? L->orbit_max : NAN, est_err, wn, pos_err, body_get_filter_resets());
        L->orbit_sum = L->orbit_max = 0;
        L->orbit_n = 0;
        L->orbit_idx++;
    }
}

int main(int argc, char** argv) {
    double orbits = argc > 1 ? atof(argv[1]) : 3.0;
    int worst = argc > 2 ? atoi(argv[2]) : 0;
    unsigned long long seed = argc > 3 ? strtoull(argv[3], NULL, 10) : 1ULL;
    double w0 = argc > 4 ? atof(argv[4]) : 0.0;
    const char* csv_path = argc > 5 ? argv[5] : "plant_sim.csv";

    /* Simulink orbit as SGP4 mean elements (near-circular, so osculating ~ mean to start) */
    const double a_km = 7021.0, mu_km = 398600.8;
    plant_config_t cfg;
    plant_default_config(&cfg);
    cfg.tle_mean[0] = (float)(sqrt(mu_km / (a_km * a_km * a_km)) * 86400.0 / (2.0 * M_PI));
    cfg.tle_mean[1] = 1e-5f;
    cfg.tle_mean[2] = (float)(60.0 * DEG);
    cfg.tle_mean[3] = (float)(137.0 * DEG);
    cfg.tle_mean[4] = (float)(349.1886 * DEG);
    cfg.tle_mean[5] = 0.0f;
    cfg.bstar = 0.0f;
    cfg.tle_epoch_jd = 2460676.5; /* 2025-01-01 00:00 UTC */
    cfg.start_jd = cfg.tle_epoch_jd;
    double q0[4] = {0.7071, 0.0, 0.7071, 0.0};
    memcpy(cfg.q0_body_to_eci, q0, sizeof(q0));
    for (int i = 0; i < 3; i++) cfg.w0_body_rad_s[i] = w0 * DEG * (i == 1 ? -1.0 : 1.0) / sqrt(3.0);
    cfg.worst_case = worst;
    cfg.seed = seed;

    const double t_orbit = 2.0 * M_PI * sqrt(a_km * a_km * a_km / mu_km);
    logger_t L = {0};
    L.t_orbit = t_orbit;
    L.n_ticks = (long)(orbits * t_orbit / PLANT_DT);
    L.orbit_idx = 1;
    L.csv = fopen(csv_path, "w");
    if (L.csv) fprintf(L.csv, "t_s,pointing,filter_resets,nadir_err_deg,est_err_deg,w_degps,pos_err_km,sunlit,"
                              "I_x,I_y,I_z\n");
    printf("orbit | pointing%% | nadir err mean/max (deg, while pointing) | est err (deg) | |w| (deg/s) | "
           "pos err (km) | resets\n");

    /* The closed loop: plant -> body() -> plant, every 0.1 s.
     * For PIL, replace plant_fsw_step_host with a function that talks to the target. */
    body_reset();
    plant_t plant;
    long done = plant_run(&plant, &cfg, L.n_ticks, plant_fsw_step_host, NULL, log_tick, &L);
    if (L.csv) fclose(L.csv);
    if (done < 0) {
        fprintf(stderr, "plant_init failed (bad TLE?)\n");
        return 1;
    }
    printf("uplinks: %ld, last TLE fit residual %.3f m, pointing %.1f%% of the run, csv: %s\n", L.uplinks,
           plant.tle_fit_resid_m, 100.0 * (double)L.pointing_ticks / (double)L.n_ticks, csv_path);
    return done == L.n_ticks ? 0 : 1;
}
