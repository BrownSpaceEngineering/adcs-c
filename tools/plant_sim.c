/*
 * Host closed loop: plant.c truth + sensors -> body() -> coil currents -> plant.c.
 * This is the PIL loop with the target replaced by a direct body() call.
 *
 *   ./build/adcs-plant-sim [orbits=3] [worst_case=0] [seed=1] [w0_degps=0] [csv=plant_sim.csv]
 *                          [coil_coupling_uT_per_A=0] [key=value ...]
 *
 * The first (up to 6) arguments may be positional, as above. Everything can also be given as
 * key=value (vectors as x,y,z), which is what the Monte Carlo driver (tools/mc/) uses:
 *   run:      orbits, worst, seed, csv (path, or "none"), quiet (1 = no per-orbit table)
 *   orbit:    alt_km, ecc, inc_deg, raan_deg, argp_deg, ma_deg, epoch_offset_days (sun geometry)
 *   initial:  w0 (deg/s along [1,-1,1]/sqrt3) or w0_vec=x,y,z (deg/s), q0=w,x,y,z (body->ECI)
 *   truth:    J=Jxx,Jyy,Jzz,Jxy,Jxz,Jyz (kg m^2), gg (0/1), nan (0/1), coupling (uT/A)
 *   errors:   sigma_mag, mag_bias, gyro_bias0, sigma_rrw, sigma_gyro, sigma_pd, pd_spread, albedo,
 *             sf_err, m_res, sigma_D, uplink_days, sigma_at   (worst=... first, then overrides)
 *   faults:   spike_t0, spike_t1, mag_spike_p, gyro_spike_p, mag_bias_step_t, mag_bias_step,
 *             coil_fault_t, coil_gain, kick_t, kick_dw (rad/s), gyro_dead_t0, gyro_dead_t1 (NaN gyro),
 *             clock_jump_t0, clock_jump_t1, clock_jump_days (flight clock offset),
 *             time_invalid_t0, time_invalid_t1 (flight software reports time_valid = false),
 *             detumble_only_t0, detumble_only_t1 (ground commands detumble-only mode)
 * Output: per-orbit table (unless quiet), the reject counters, and one line
 *   SUMMARY key=value ...   (see print_summary) for scripts.
 *
 * Default initial conditions follow the Simulink model (asbCubeSat): 2025-01-01 00:00 UTC,
 * a = 7021 km, e = 1e-5, i = 60 deg, RAAN = 137 deg, argp = 349.1886 deg, M = 0,
 * attitude [.7071 0 .7071 0].
 */
#include "include/body.h"
#include "include/plant.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEG (M_PI / 180.0)
#define MU_KM 398600.8

static double angle_between(const double* a, const double* b) {
    double c = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
               / (sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) * sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]));
    return acos(c > 1 ? 1 : (c < -1 ? -1 : c));
}

/* ---------------- per-tick logging / scoring (plant_run's on_tick) ---------------- */
typedef struct {
    FILE* csv;
    bool quiet;
    double t_orbit;
    long n_ticks;
    double orbit_sum, orbit_max;
    long orbit_n, uplinks, pointing_ticks;
    int orbit_idx;
    /* 1 Hz history for the summary */
    long n_hist;
    float *nadir, *est, *rate;
    unsigned char* pointing;
    double t_handoff;
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
    if (pt && isnan(L->t_handoff)) L->t_handoff = tr->t_s;
    if (pt) {
        L->pointing_ticks++;
        L->orbit_sum += nadir_err;
        L->orbit_max = fmax(L->orbit_max, nadir_err);
        L->orbit_n++;
    }
    if (k % 10 == 0 && L->n_hist < L->n_ticks / 10 + 1) {
        long j = L->n_hist++;
        L->nadir[j] = (float)nadir_err;
        L->est[j] = (float)est_err;
        L->rate[j] = (float)wn;
        L->pointing[j] = pt;
    }
    if (L->csv && k % 100 == 0) {
        fprintf(L->csv, "%.1f,%d,%d,%.4f,%.4f,%.5f,%.4f,%d,%.5f,%.5f,%.5f\n", tr->t_s, pt, body_get_filter_resets(),
                nadir_err, est_err, wn, pos_err, tr->sunlit, I[0], I[1], I[2]);
    }
    if (tr->t_s >= L->orbit_idx * L->t_orbit || k == L->n_ticks - 1) {
        if (!L->quiet) {
            printf("%5d | %8.1f | %15.2f / %-22.2f | %13.3f | %11.4f | %12.3f | %d\n", L->orbit_idx,
                   100.0 * (double)L->orbit_n / (L->t_orbit / PLANT_DT),
                   L->orbit_n ? L->orbit_sum / (double)L->orbit_n : NAN, L->orbit_n ? L->orbit_max : NAN, est_err, wn,
                   pos_err, body_get_filter_resets());
        }
        L->orbit_sum = L->orbit_max = 0;
        L->orbit_n = 0;
        L->orbit_idx++;
    }
}

static int cmp_float(const void* a, const void* b) {
    float x = *(const float*)a, y = *(const float*)b;
    return (x > y) - (x < y);
}

/* RMS / p95 / max of x over 1 Hz samples [i0, n), only where pointing if need_pt, finite only */
static void stats(const logger_t* L, const float* x, long i0, bool need_pt, double* rms, double* p95, double* mx) {
    float* buf = malloc(sizeof(float) * (size_t)(L->n_hist + 1));
    long m = 0;
    double ss = 0.0;
    for (long j = i0; j < L->n_hist; j++) {
        if ((need_pt && !L->pointing[j]) || !isfinite(x[j])) continue;
        buf[m++] = x[j];
        ss += (double)x[j] * x[j];
    }
    if (m == 0) {
        *rms = *p95 = *mx = NAN;
    } else {
        qsort(buf, (size_t)m, sizeof(float), cmp_float);
        *rms = sqrt(ss / (double)m);
        *p95 = buf[(long)floor(0.95 * (double)(m - 1))];
        *mx = buf[m - 1];
    }
    free(buf);
}


/* Photodiode rejects: total over all diodes, and the per-diode list "n0,n1,...,n17" (no spaces) */
static unsigned long pd_rejects_total(const body_sensor_rejects_t* rj) {
    unsigned long n = 0;
    for (int k = 0; k < NUM_DIODES; k++) n += rj->pd_rejected[k];
    return n;
}
static void pd_rejects_list(const body_sensor_rejects_t* rj, char* buf, size_t len) {
    size_t used = 0;
    buf[0] = '\0';
    for (int k = 0; k < NUM_DIODES && used < len; k++) {
        int w = snprintf(buf + used, len - used, k ? ",%u" : "%u", (unsigned)rj->pd_rejected[k]);
        if (w < 0) break;
        used += (size_t)w;
    }
}

static void print_summary(const logger_t* L, const plant_t* p, double t_settle) {
    /* settled window: from handoff + t_settle to the end (1 Hz samples, index = seconds) */
    long i_set = isnan(L->t_handoff) ? L->n_hist : (long)(L->t_handoff + t_settle);
    if (i_set > L->n_hist) i_set = L->n_hist;
    double nad_rms, nad_p95, nad_max, est_rms, est_p95, est_max, w_rms, w_p95, w_max;
    stats(L, L->nadir, i_set, false, &nad_rms, &nad_p95, &nad_max);
    stats(L, L->est, i_set, true, &est_rms, &est_p95, &est_max);
    stats(L, L->rate, i_set, false, &w_rms, &w_p95, &w_max);
    long n_pt = 0, n_set = 0;
    for (long j = i_set; j < L->n_hist; j++) {
        n_pt += L->pointing[j];
        n_set++;
    }
    /* acquisition: time from the first handoff until pointing holds with nadir < 10 deg for 600 s */
    double t_acq = NAN;
    if (!isnan(L->t_handoff)) {
        long run = 0;
        for (long j = (long)L->t_handoff; j < L->n_hist; j++) {
            run = (L->pointing[j] && L->nadir[j] < 10.0f) ? run + 1 : 0;
            if (run >= 600) {
                t_acq = (double)(j - 599) - L->t_handoff;
                break;
            }
        }
    }
    body_sensor_rejects_t rj;
    body_get_sensor_rejects(&rj);
    char pd_list[NUM_DIODES * 11 + 1];
    pd_rejects_list(&rj, pd_list, sizeof(pd_list));
    float mb[3];
    bool mb_ok = body_get_mag_bias(mb);
    printf("SUMMARY magcal_ready=%d magcal_uT=%.3f,%.3f,%.3f ", mb_ok, mb[0] * 1e6f, mb[1] * 1e6f, mb[2] * 1e6f);
    float csign[3];
    bool cen[3];
    body_get_coil_config(csign, cen);
    printf("coil_sign=%.0f,%.0f,%.0f coil_on=%d,%d,%d ", csign[0], csign[1], csign[2], cen[0], cen[1], cen[2]);
    printf("handoff_s=%.1f acquire_s=%.1f point_frac=%.4f point_frac_settled=%.4f nadir_rms=%.3f "
           "nadir_p95=%.3f nadir_max=%.2f est_rms=%.3f est_p95=%.3f rate_rms_degps=%.4f resets=%d uplinks=%ld "
           "tle_resid_m=%.2f rej_gyro_stuck=%u rej_gyro_spike=%u rej_mag_stuck=%u rej_mag_norm=%u rej_mag_nan=%u rej_pd=%lu "
           "rej_pd_diodes=%s rej_gated=%u tle_accepted=%u tle_rejected=%u tle_last_reject=%u orbit_fail_streak=%u orbit_fail_reason=%u\n",
           L->t_handoff, t_acq, (double)L->pointing_ticks / (double)L->n_ticks,
           n_set ? (double)n_pt / (double)n_set : NAN, nad_rms, nad_p95, nad_max, est_rms, est_p95, w_rms,
           body_get_filter_resets(), L->uplinks, p->tle_fit_resid_m, rj.gyro_stuck, rj.gyro_spike, rj.mag_stuck,
           rj.mag_norm, rj.mag_nan, pd_rejects_total(&rj), pd_list, rj.innovation_gated, rj.tle_accepted, rj.tle_rejected, rj.tle_last_reject, rj.orbit_fail_streak,
           rj.orbit_fail_reason);
}

/* ---------------- argument parsing ---------------- */
static int parse_vec(const char* s, double* out, int n) {
    char* end;
    for (int i = 0; i < n; i++) {
        out[i] = strtod(s, &end);
        if (end == s) return -1;
        s = end;
        if (i < n - 1) {
            if (*s != ',') return -1;
            s++;
        }
    }
    return 0;
}

typedef struct {
    double orbits, w0, coupling, alt_km, ecc, inc, raan, argp, ma, epoch_offset_days;
    int worst, quiet, gg, nan;
    unsigned long long seed;
    const char* csv;
    bool have_w0_vec, have_q0, have_J;
    double w0_vec[3], q0[4], J6[6];
} opts_t;

static const char* RUN_KEYS[] = {"orbits",   "worst",   "seed",  "w0",       "csv",    "coupling", "quiet",
                                 "gg",       "nan",     "alt_km", "ecc",     "inc_deg", "raan_deg", "argp_deg",
                                 "ma_deg",   "epoch_offset_days", "w0_vec", "q0", "J", NULL};

static bool is_run_key(const char* key) {
    for (int i = 0; RUN_KEYS[i] != NULL; i++) {
        if (!strcmp(key, RUN_KEYS[i])) return true;
    }
    return false;
}

/* Pass 1: run / orbit / initial-condition keys */
static int apply_run_kv(const char* key, const char* val, opts_t* o) {
    if (!strcmp(key, "orbits")) o->orbits = atof(val);
    else if (!strcmp(key, "worst")) o->worst = atoi(val);
    else if (!strcmp(key, "seed")) o->seed = strtoull(val, NULL, 10);
    else if (!strcmp(key, "w0")) o->w0 = atof(val);
    else if (!strcmp(key, "csv")) o->csv = val;
    else if (!strcmp(key, "coupling")) o->coupling = atof(val);
    else if (!strcmp(key, "quiet")) o->quiet = atoi(val);
    else if (!strcmp(key, "gg")) o->gg = atoi(val);
    else if (!strcmp(key, "nan")) o->nan = atoi(val);
    else if (!strcmp(key, "alt_km")) o->alt_km = atof(val);
    else if (!strcmp(key, "ecc")) o->ecc = atof(val);
    else if (!strcmp(key, "inc_deg")) o->inc = atof(val);
    else if (!strcmp(key, "raan_deg")) o->raan = atof(val);
    else if (!strcmp(key, "argp_deg")) o->argp = atof(val);
    else if (!strcmp(key, "ma_deg")) o->ma = atof(val);
    else if (!strcmp(key, "epoch_offset_days")) o->epoch_offset_days = atof(val);
    else if (!strcmp(key, "w0_vec")) { if (parse_vec(val, o->w0_vec, 3)) return -1; o->have_w0_vec = true; }
    else if (!strcmp(key, "q0")) { if (parse_vec(val, o->q0, 4)) return -1; o->have_q0 = true; }
    else if (!strcmp(key, "J")) { if (parse_vec(val, o->J6, 6)) return -1; o->have_J = true; }
    return 0;
}

/* Pass 2 (after the worst-case error model is loaded): error-model / fault overrides */
static int apply_truth_kv(const char* key, const char* val, plant_config_t* cfg) {
    plant_errors_t* e = &cfg->err;
    plant_faults_t* f = &cfg->faults;
    double* v3 = NULL;
    if (!strcmp(key, "sigma_mag")) e->sigma_mag = atof(val);
    else if (!strcmp(key, "mag_bias")) v3 = e->mag_bias;
    else if (!strcmp(key, "gyro_bias0")) v3 = e->gyro_bias0;
    else if (!strcmp(key, "sigma_rrw")) e->sigma_rrw = atof(val);
    else if (!strcmp(key, "sigma_gyro")) v3 = e->sigma_gyro;
    else if (!strcmp(key, "sigma_pd")) e->sigma_pd = atof(val);
    else if (!strcmp(key, "pd_spread")) e->pd_spread = atof(val);
    else if (!strcmp(key, "albedo")) e->albedo = atof(val);
    else if (!strcmp(key, "sf_err")) v3 = e->sf_err;
    else if (!strcmp(key, "m_res")) v3 = e->m_res;
    else if (!strcmp(key, "sigma_D")) e->sigma_D = atof(val);
    else if (!strcmp(key, "uplink_days")) e->uplink_period_s = atof(val) * 86400.0;
    else if (!strcmp(key, "sigma_at")) e->sigma_at_s = atof(val);
    else if (!strcmp(key, "spike_t0")) f->spike_t0 = atof(val);
    else if (!strcmp(key, "spike_t1")) f->spike_t1 = atof(val);
    else if (!strcmp(key, "mag_spike_p")) f->mag_spike_p = atof(val);
    else if (!strcmp(key, "gyro_spike_p")) f->gyro_spike_p = atof(val);
    else if (!strcmp(key, "mag_bias_step_t")) f->mag_bias_step_t = atof(val);
    else if (!strcmp(key, "mag_bias_step")) v3 = f->mag_bias_step;
    else if (!strcmp(key, "coil_fault_t")) f->coil_fault_t = atof(val);
    else if (!strcmp(key, "gyro_dead_t0")) f->gyro_dead_t0 = atof(val);
    else if (!strcmp(key, "gyro_dead_t1")) f->gyro_dead_t1 = atof(val);
    else if (!strcmp(key, "clock_jump_t0")) f->clock_jump_t0 = atof(val);
    else if (!strcmp(key, "clock_jump_t1")) f->clock_jump_t1 = atof(val);
    else if (!strcmp(key, "clock_jump_days")) f->clock_jump_days = atof(val);
    else if (!strcmp(key, "time_invalid_t0")) f->time_invalid_t0 = atof(val);
    else if (!strcmp(key, "time_invalid_t1")) f->time_invalid_t1 = atof(val);
    else if (!strcmp(key, "detumble_only_t0")) f->detumble_only_t0 = atof(val);
    else if (!strcmp(key, "detumble_only_t1")) f->detumble_only_t1 = atof(val);
    else if (!strcmp(key, "coil_gain")) v3 = f->coil_gain;
    else if (!strcmp(key, "kick_t")) { cfg->kick_on = true; cfg->kick_t_s = atof(val); }
    else if (!strcmp(key, "kick_dw")) v3 = cfg->kick_dw;
    else if (!is_run_key(key)) {
        fprintf(stderr, "unknown option %s\n", key);
        return -1;
    }
    if (v3 != NULL && parse_vec(val, v3, 3)) return -1;
    return 0;
}

int main(int argc, char** argv) {
    opts_t o = {.orbits = 3.0, .seed = 1ULL, .csv = "plant_sim.csv", .gg = 1, .nan = 1, .alt_km = 7021.0 - 6378.137,
                .ecc = 1e-5, .inc = 60.0, .raan = 137.0, .argp = 349.1886, .ma = 0.0};
    plant_config_t cfg;
    plant_default_config(&cfg);
    /* pass 1: positional + run/orbit/initial keys */
    const char* positional[6] = {"orbits", "worst", "seed", "w0", "csv", "coupling"};
    int npos = 0;
    for (int a = 1; a < argc; a++) {
        char key[64];
        const char* eq = strchr(argv[a], '=');
        const char* val;
        if (eq == NULL) {
            if (npos >= 6) {
                fprintf(stderr, "too many positional arguments\n");
                return 2;
            }
            snprintf(key, sizeof(key), "%s", positional[npos++]);
            val = argv[a];
        } else {
            snprintf(key, sizeof(key), "%.*s", (int)(eq - argv[a]), argv[a]);
            val = eq + 1;
        }
        if (apply_run_kv(key, val, &o)) {
            fprintf(stderr, "bad value for %s\n", key);
            return 2;
        }
    }
    /* error model for worst, then pass 2: error-model / fault overrides */
    cfg.worst_case = o.worst;
    plant_error_model(o.worst, &cfg.err);
    for (int a = 1; a < argc; a++) {
        const char* eq = strchr(argv[a], '=');
        if (eq == NULL) continue;
        char key[64];
        snprintf(key, sizeof(key), "%.*s", (int)(eq - argv[a]), argv[a]);
        if (apply_truth_kv(key, eq + 1, &cfg)) {
            fprintf(stderr, "bad value for %s\n", key);
            return 2;
        }
    }

    /* orbit as SGP4 mean elements (near-circular, so osculating ~ mean to start) */
    double a_km = 6378.137 + o.alt_km;
    cfg.tle_mean[0] = (float)(sqrt(MU_KM / (a_km * a_km * a_km)) * 86400.0 / (2.0 * M_PI));
    cfg.tle_mean[1] = (float)o.ecc;
    cfg.tle_mean[2] = (float)(o.inc * DEG);
    cfg.tle_mean[3] = (float)(o.raan * DEG);
    cfg.tle_mean[4] = (float)(o.argp * DEG);
    cfg.tle_mean[5] = (float)(o.ma * DEG);
    cfg.bstar = 0.0f;
    cfg.tle_epoch_jd = 2460676.5 + o.epoch_offset_days; /* 2025-01-01 00:00 UTC + offset */
    cfg.start_jd = cfg.tle_epoch_jd;
    double q0[4] = {0.7071, 0.0, 0.7071, 0.0};
    if (o.have_q0) memcpy(q0, o.q0, sizeof(q0));
    memcpy(cfg.q0_body_to_eci, q0, sizeof(q0));
    for (int i = 0; i < 3; i++) {
        cfg.w0_body_rad_s[i] = o.have_w0_vec ? o.w0_vec[i] * DEG : o.w0 * DEG * (i == 1 ? -1.0 : 1.0) / sqrt(3.0);
    }
    cfg.seed = o.seed;
    cfg.gravity_gradient = o.gg != 0;
    cfg.nan_dropouts = o.nan != 0;
    if (o.have_J) {
        const double* j = o.J6;
        double J[3][3] = {{j[0], j[3], j[4]}, {j[3], j[1], j[5]}, {j[4], j[5], j[2]}};
        memcpy(cfg.J, J, sizeof(J));
    }
    /* coupling: mostly along the coil axis, with some cross-talk (geometry-dependent) */
    const double shape[3][3] = {{1.0, 0.3, -0.2}, {0.2, 1.0, 0.3}, {-0.3, 0.2, 1.0}};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) cfg.mag_coil_coupling[i][j] = o.coupling * 1e-6 * shape[i][j];

    const double t_orbit = 2.0 * M_PI * sqrt(a_km * a_km * a_km / MU_KM);
    logger_t L = {0};
    L.quiet = o.quiet != 0;
    L.t_orbit = t_orbit;
    L.n_ticks = (long)(o.orbits * t_orbit / PLANT_DT);
    L.orbit_idx = 1;
    L.t_handoff = NAN;
    size_t nh = (size_t)(L.n_ticks / 10 + 2);
    L.nadir = malloc(nh * sizeof(float));
    L.est = malloc(nh * sizeof(float));
    L.rate = malloc(nh * sizeof(float));
    L.pointing = malloc(nh);
    if (strcmp(o.csv, "none") != 0) {
        L.csv = fopen(o.csv, "w");
        if (L.csv) {
            fprintf(L.csv, "t_s,pointing,filter_resets,nadir_err_deg,est_err_deg,w_degps,pos_err_km,sunlit,I_x,I_y,I_z\n");
        }
    }
    if (!L.quiet) {
        printf("orbit | pointing%% | nadir err mean/max (deg, while pointing) | est err (deg) | |w| (deg/s) | "
               "pos err (km) | resets\n");
    }

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
    if (!L.quiet) {
        printf("uplinks: %ld, last TLE fit residual %.3f m, pointing %.1f%% of the run, csv: %s\n", L.uplinks,
               plant.tle_fit_resid_m, 100.0 * (double)L.pointing_ticks / (double)L.n_ticks, o.csv);
        body_sensor_rejects_t rj;
        body_get_sensor_rejects(&rj);
        char pd_list[NUM_DIODES * 11 + 1];
        pd_rejects_list(&rj, pd_list, sizeof(pd_list));
        printf("sensor rejects: gyro stuck %u, gyro spike %u, mag stuck %u, mag |B| %u, mag NaN %u, photodiode %lu "
               "(per diode %s), innovation gated %u, TLE accepted %u / rejected %u\n",
               rj.gyro_stuck, rj.gyro_spike, rj.mag_stuck, rj.mag_norm, rj.mag_nan, pd_rejects_total(&rj), pd_list,
               rj.innovation_gated, rj.tle_accepted, rj.tle_rejected);
    }
    print_summary(&L, &plant, t_orbit); /* settled = one orbit after the hand-off */
    free(L.nadir);
    free(L.est);
    free(L.rate);
    free(L.pointing);
    return done == L.n_ticks ? 0 : 1;
}
