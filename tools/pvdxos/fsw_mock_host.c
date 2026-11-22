/* Host replica of PVDXos PLANT_SIM: adcs_mock.c plant config + the new adcs_task.c control/TLE glue.
 * usage: fsw_mock_host <seconds> <no_utc_until_s> <gyro_dead_until_s> <csv|none> */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "include/body.h"
#include "include/plant.h"
#include "include/bdot.h"
#include "include/moments2currents.h"
#include "include/sensor_frames.h"
extern const float BDOT_GAIN; extern const float* Imax;
static int OLD = 0; static int attempts = 0;

#define MOCK_PI 3.14159265358979323846
#define RAD2DEG (180.0 / MOCK_PI)

static float pend_mean[6]; static int tle_pending = 0; static int dropped = 0, saved = 0;

static void resolve(const body_sensor_rejects_t* b, const body_sensor_rejects_t* a) {
    if (!tle_pending) return;
    if (a->tle_accepted != b->tle_accepted) { tle_pending = 0; saved++; }
    else if (a->tle_rejected != b->tle_rejected) { tle_pending = 0; dropped++; printf("  TLE dropped, reason %u\n", a->tle_last_reject); }
}

static double nadir_err(const plant_truth_t* t) {
    const double* q = t->q_body_to_eci;
    double z[3] = {2*(q[1]*q[3]+q[0]*q[2]), 2*(q[2]*q[3]-q[0]*q[1]), q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3]};
    double n[3] = {-t->r_eci_m[0], -t->r_eci_m[1], -t->r_eci_m[2]};
    double c = (z[0]*n[0]+z[1]*n[1]+z[2]*n[2]) / (sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]) * sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]));
    return acos(c > 1 ? 1 : (c < -1 ? -1 : c)) * RAD2DEG;
}

int main(int argc, char** argv) {
    double T = argc > 1 ? atof(argv[1]) : 20000, no_utc = argc > 2 ? atof(argv[2]) : 0, gdead = argc > 3 ? atof(argv[3]) : 0;
    const char* csvp = argc > 4 ? argv[4] : "none";
    OLD = argc > 5 && !strcmp(argv[5], "old");
    plant_config_t cfg; plant_default_config(&cfg);
    const double a = 7021.0;   /* adcs_mock.c */
    cfg.tle_mean[0] = (float)(sqrt(398600.8 / (a*a*a)) * 86400.0 / (2.0*MOCK_PI));
    cfg.tle_mean[1] = 1e-5f; cfg.tle_mean[2] = (float)(60.0/RAD2DEG); cfg.tle_mean[3] = (float)(137.0/RAD2DEG);
    cfg.tle_mean[4] = (float)(349.1886/RAD2DEG); cfg.tle_mean[5] = 0.0f; cfg.bstar = 0.0f;
    cfg.tle_epoch_jd = 2461345.5; cfg.start_jd = 2461345.5;
    double q0[4] = {0.7071, 0.0, 0.7071, 0.0}; memcpy(cfg.q0_body_to_eci, q0, sizeof q0);
    for (int i = 0; i < 3; i++) cfg.w0_body_rad_s[i] = 5.0/RAD2DEG * (i == 1 ? -1.0 : 1.0) / sqrt(3.0);
    cfg.seed = 1u;
    plant_t p; if (plant_init(&p, &cfg) != 0) { printf("plant_init failed\n"); return 1; }
    body_reset();
    FILE* csv = strcmp(csvp, "none") ? fopen(csvp, "w") : NULL;
    if (csv) fprintf(csv, "t_s,pointing,resets,w_degps,nadir_deg\n");
    plant_body_inputs_t in; plant_truth_t tr;
    long n = 0, npt = 0, nset = 0; double handoff = NAN, se = 0; int uplinks = 0;
    for (double t = 0; t < T; t += PLANT_DT, n++) {
        plant_sense(&p, &in);
        if (in.has_posn_update) { body_set_tle_metadata(in.bstar, in.tle_epoch_jd); memcpy(pend_mean, in.posn_update, sizeof pend_mean); tle_pending = 1; uplinks++; }
        if (t < gdead) for (int i = 0; i < 3; i++) in.gyro[i] = NAN;
        int have_time = !(t < no_utc);
        float currents[3];
        body_sensor_rejects_t b, aft; body_get_sensor_rejects(&b);
        body(in.last_magnetometer, in.magnetometer, in.gyro, in.photodiodes, (have_time && tle_pending) ? pend_mean : NULL,
             in.unix_time, in.jd_scalar, in.jd_frac, have_time, in.dt, currents);
        body_get_sensor_rejects(&aft);
        if (!OLD) resolve(&b, &aft);
        else if (tle_pending) {
            float m[6], bs; double ep;
            int gf = isfinite(in.gyro[0]) && isfinite(in.gyro[1]) && isfinite(in.gyro[2]);
            if (body_get_tle(m, &bs, &ep) && !memcmp(m, pend_mean, sizeof m)) { tle_pending = 0; saved++; attempts = 0; }
            else if (have_time && gf && ++attempts >= 3) { tle_pending = 0; dropped++; attempts = 0; }
        }
        if (OLD && !body_is_pointing()) {   /* adcs_detumble() override */
            float lb[3], mb[3], mom[3];
            int ok = 1; for (int i = 0; i < 3; i++) ok &= isfinite(in.last_magnetometer[i]) && isfinite(in.magnetometer[i]);
            if (!ok) { currents[0] = currents[1] = currents[2] = 0; }
            else { mag_sensor_to_body(in.last_magnetometer, lb); mag_sensor_to_body(in.magnetometer, mb); Bdot(mb, lb, BDOT_GAIN, in.dt, mom); moment2current3axis(mom, Imax, currents); }
        }
        plant_step(&p, currents);
        plant_get_truth(&p, &tr);
        int pt = body_is_pointing(); npt += pt;
        if (pt && isnan(handoff)) handoff = t;
        double ne = nadir_err(&tr);
        if (!isnan(handoff) && t > handoff + 5800) { se += ne*ne; nset++; }
        if (csv && n % 100 == 0) {
            const double* w = tr.w_body_rad_s;
            fprintf(csv, "%.1f,%d,%d,%.4f,%.3f\n", t, pt, body_get_filter_resets(), sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*RAD2DEG, ne);
        }
    }
    body_sensor_rejects_t r; body_get_sensor_rejects(&r);
    printf("handoff_s=%.0f point_frac=%.3f nadir_rms_settled=%.2f resets=%d uplinks=%d tle_saved=%d tle_dropped=%d tle_rejected=%u orbit_fail_reason=%u gyro_fail_streak=%u\n",
           handoff, (double)npt/n, nset ? sqrt(se/nset) : NAN, body_get_filter_resets(), uplinks, saved, dropped, r.tle_rejected, r.orbit_fail_reason, r.gyro_fail_streak);
    if (csv) fclose(csv);
    return 0;
}
