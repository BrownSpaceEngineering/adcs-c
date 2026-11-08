/* Truth plant for PIL / host closed-loop testing of body(). NOT flight code. See include/plant.h. */
#include "include/plant.h"
#include "include/body.h"
#include "include/magnetosphere.h"
#include "include/photodiode_determination.h"
#include "include/sgp4_init.h"
#include "include/sgp4_propagate.h"
#include "include/sunvec.h"
#include "include/tle_fit.h"
#include <math.h>
#include <string.h>

/* ---------------- constants ---------------- */
#define MU_E 3.986004418e14  /* m^3/s^2 */
#define R_E 6378137.0        /* m (eclipse / albedo geometry, as in the Simulink sensor block) */
#define J2_E 1.08262668e-3
#define JD_UNIX_EPOCH 2440587.5

/* Spacecraft Dynamics inertia (kg m^2), symmetric part of the block value */
static const double J_SC[3][3] = {{0.029730, 24.4275e-6, 3672.061e-6},
                                  {24.4275e-6, 36414.034e-6, -285.496e-6},
                                  {3672.061e-6, -285.496e-6, 12429.899e-6}};
/* Coil gains body() assumes (moments2currents.c), A m^2 per A */
static const double K_COIL[3] = {2.6, 2.6, 1.02};
#define PD_MAX 1.7 /* MAX_READING in photodiode_determination.c */

/* ---------------- error models (Simulink "Sensors (noise)" / "ADCS_Actuators") ---------------- */
typedef struct {
    double sigma_mag, mag_bias[3];
    double gyro_bias0[3], sigma_rrw, sigma_gyro[3];
    double sigma_pd, pd_spread, albedo;
    double sf_err[3], m_res[3], sigma_D;
    double uplink_period_s, sigma_at_s;
} err_model_t;

static void error_model(int worst_case, err_model_t* e) {
    if (worst_case == 0) {
        *e = (err_model_t){.sigma_mag = 30e-9,
                           .mag_bias = {0, 0, 0},
                           .gyro_bias0 = {0.002, -0.001, 0.0015},
                           .sigma_rrw = 1.4e-7,
                           .sigma_gyro = {2.5e-5, 2.5e-5, 3.78e-5},
                           .sigma_pd = 0.01,
                           .pd_spread = 0.0,
                           .albedo = 0.0,
                           .sf_err = {0.08, -0.06, 0.10},
                           .m_res = {4e-3, -3e-3, 5e-3},
                           .sigma_D = 1e-3,
                           .uplink_period_s = 86400.0,
                           .sigma_at_s = 0.07};
    } else {
        *e = (err_model_t){.sigma_mag = 50e-9,
                           .mag_bias = {350e-9, -300e-9, 450e-9},
                           .gyro_bias0 = {0.005, -0.0025, 0.00375},
                           .sigma_rrw = 5e-7,
                           .sigma_gyro = {7.5e-5, 7.5e-5, 1.13e-4},
                           .sigma_pd = 0.02,
                           .pd_spread = 0.05,
                           .albedo = 0.4,
                           .sf_err = {0.20, -0.20, 0.20},
                           .m_res = {15e-3, -12e-3, 15e-3},
                           .sigma_D = 0.0055,
                           .uplink_period_s = 3.0 * 86400.0,
                           .sigma_at_s = 0.27};
    }
}

/* NaN dropouts, per tick: {p_glitch, p_long, long_min, long_max} (ticks) */
static const double MAG_NAN[4] = {1e-3, 2e-5, 10, 50};
static const double GYRO_NAN[4] = {5e-4, 1e-5, 10, 30};
static const double PD_NAN[4] = {2e-4, 2e-6, 10, 60};

/* ---------------- RNG: xoshiro256** + Box-Muller (reproducible, no libc rand) ---------------- */
static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static uint64_t rng_next(plant_t* p) {
    uint64_t* s = p->rng;
    uint64_t result = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

static void rng_seed(plant_t* p, uint64_t seed) {
    for (int i = 0; i < 4; i++) { /* splitmix64 */
        uint64_t z = (seed += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        p->rng[i] = z ^ (z >> 31);
    }
    p->have_spare_normal = false;
}

static double rand_u(plant_t* p) { return (double)(rng_next(p) >> 11) * 0x1.0p-53; } /* [0,1) */

static double randn(plant_t* p) {
    if (p->have_spare_normal) {
        p->have_spare_normal = false;
        return p->spare_normal;
    }
    double u1 = rand_u(p), u2 = rand_u(p);
    if (u1 < 1e-300) u1 = 1e-300;
    double r = sqrt(-2.0 * log(u1));
    p->spare_normal = r * sin(2.0 * M_PI * u2);
    p->have_spare_normal = true;
    return r * cos(2.0 * M_PI * u2);
}

/* ---------------- small vector / quaternion helpers (double) ---------------- */
static double dot(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static double norm3(const double* a) { return sqrt(dot(a, a)); }
static void crossd(const double* a, const double* b, double* c) {
    double t[3] = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    memcpy(c, t, sizeof(t));
}
static void matvec(const double M[3][3], const double* x, double* y) {
    for (int i = 0; i < 3; i++) y[i] = M[i][0] * x[0] + M[i][1] * x[1] + M[i][2] * x[2];
}
/* R(q) for active body->ECI q: v_eci = R v_body */
static void quat_to_R(const double* q, double R[3][3]) {
    double w = q[0], x = q[1], y = q[2], z = q[3];
    R[0][0] = w * w + x * x - y * y - z * z;
    R[0][1] = 2 * (x * y - w * z);
    R[0][2] = 2 * (x * z + w * y);
    R[1][0] = 2 * (x * y + w * z);
    R[1][1] = w * w - x * x + y * y - z * z;
    R[1][2] = 2 * (y * z - w * x);
    R[2][0] = 2 * (x * z - w * y);
    R[2][1] = 2 * (y * z + w * x);
    R[2][2] = w * w - x * x - y * y + z * z;
}
static void eci_to_body(const double* q, const double* v_eci, double* v_body) {
    double R[3][3];
    quat_to_R(q, R);
    for (int i = 0; i < 3; i++) v_body[i] = R[0][i] * v_eci[0] + R[1][i] * v_eci[1] + R[2][i] * v_eci[2];
}
static void normalize_q(double* q) {
    double n = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    for (int i = 0; i < 4; i++) q[i] /= n;
}

/* ---------------- environment ---------------- */
static double jd_now(const plant_t* p) { return p->cfg.start_jd + p->t / 86400.0; }

static void field_eci(const plant_t* p, const double* r, double* b) {
    double jd = jd_now(p);
    double jd_int = floor(jd);
    float rf[3] = {(float)r[0], (float)r[1], (float)r[2]}, bf[3];
    wmm_eci_embedded_v2(rf, (int32_t)jd_int, (float)(jd - jd_int), bf);
    for (int i = 0; i < 3; i++) b[i] = bf[i];
}

static int unix_time_of(double jd) { return (int)lround((jd - JD_UNIX_EPOCH) * 86400.0); }

static void sun_eci(const plant_t* p, double* s) {
    float sf[3];
    sun_vec(unix_time_of(jd_now(p)), sf);
    for (int i = 0; i < 3; i++) s[i] = sf[i];
    double n = norm3(s);
    for (int i = 0; i < 3; i++) s[i] /= n;
}

/* Same cone test as SunInView / the Simulink sensor block */
static bool sunlit(const double* r, const double* s) {
    double rn = norm3(r);
    double c = -dot(r, s) / rn;
    c = c > 1 ? 1 : (c < -1 ? -1 : c);
    double rho = asin(fmin(1.0, R_E / rn));
    return acos(c) > rho;
}

/* ---------------- dynamics ---------------- */
static void orbit_accel(const double* r, double* a) {
    double rn = norm3(r), rn2 = rn * rn;
    double k = -MU_E / (rn2 * rn);
    double f = 1.5 * J2_E * (R_E * R_E) / rn2;
    double z2 = r[2] * r[2] / rn2;
    a[0] = k * r[0] * (1.0 + f * (1.0 - 5.0 * z2));
    a[1] = k * r[1] * (1.0 + f * (1.0 - 5.0 * z2));
    a[2] = k * r[2] * (1.0 + f * (3.0 - 5.0 * z2));
}

/* RK4 on the orbit, one step of h */
static void orbit_rk4(double* r, double* v, double h) {
    double k1r[3], k1v[3], k2r[3], k2v[3], k3r[3], k3v[3], k4r[3], k4v[3], rt[3], vt[3];
    memcpy(k1r, v, sizeof(k1r));
    orbit_accel(r, k1v);
    for (int i = 0; i < 3; i++) { rt[i] = r[i] + 0.5 * h * k1r[i]; vt[i] = v[i] + 0.5 * h * k1v[i]; }
    memcpy(k2r, vt, sizeof(k2r));
    orbit_accel(rt, k2v);
    for (int i = 0; i < 3; i++) { rt[i] = r[i] + 0.5 * h * k2r[i]; vt[i] = v[i] + 0.5 * h * k2v[i]; }
    memcpy(k3r, vt, sizeof(k3r));
    orbit_accel(rt, k3v);
    for (int i = 0; i < 3; i++) { rt[i] = r[i] + h * k3r[i]; vt[i] = v[i] + h * k3v[i]; }
    memcpy(k4r, vt, sizeof(k4r));
    orbit_accel(rt, k4v);
    for (int i = 0; i < 3; i++) {
        r[i] += h / 6.0 * (k1r[i] + 2 * k2r[i] + 2 * k3r[i] + k4r[i]);
        v[i] += h / 6.0 * (k1v[i] + 2 * k2v[i] + 2 * k3v[i] + k4v[i]);
    }
}

/* J^-1 by adjugate: inv[i][j] = cofactor[j][i] / det */
static void inertia_inverse(double inv[3][3]) {
    const double(*J)[3] = J_SC;
    double cof[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int i1 = (i + 1) % 3, i2 = (i + 2) % 3, j1 = (j + 1) % 3, j2 = (j + 2) % 3;
            cof[i][j] = J[i1][j1] * J[i2][j2] - J[i1][j2] * J[i2][j1];
        }
    }
    double det = J[0][0] * cof[0][0] + J[0][1] * cof[0][1] + J[0][2] * cof[0][2];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) inv[i][j] = cof[j][i] / det;
}

/* Attitude derivative. ext = magnetorquer (+ kick) torque, held over the step;
 * gravity gradient is re-evaluated from the attitude at each stage. */
static void att_deriv(const double* q, const double* w, const double* ext, const double* r_hat_eci,
                      double n2, bool gg, double* dq, double* dw) {
    dq[0] = -0.5 * (q[1] * w[0] + q[2] * w[1] + q[3] * w[2]);
    dq[1] = 0.5 * (q[0] * w[0] + q[2] * w[2] - q[3] * w[1]);
    dq[2] = 0.5 * (q[0] * w[1] + q[3] * w[0] - q[1] * w[2]);
    dq[3] = 0.5 * (q[0] * w[2] + q[1] * w[1] - q[2] * w[0]);
    double T[3] = {ext[0], ext[1], ext[2]};
    if (gg) {
        double rb[3], Jr[3], tg[3];
        eci_to_body(q, r_hat_eci, rb);
        matvec(J_SC, rb, Jr);
        crossd(rb, Jr, tg);
        for (int i = 0; i < 3; i++) T[i] += 3.0 * n2 * tg[i];
    }
    double Jw[3], wJw[3], rhs[3];
    matvec(J_SC, w, Jw);
    crossd(w, Jw, wJw);
    for (int i = 0; i < 3; i++) rhs[i] = T[i] - wJw[i];
    double Jinv[3][3];
    inertia_inverse(Jinv);
    matvec(Jinv, rhs, dw);
}

/* ---------------- sensors ---------------- */
/* n = ticks left in the current outage; returns whether this tick's sample is bad */
static bool outage(plant_t* p, int* n, const double* cfg) {
    if (*n == 0) {
        double r = rand_u(p);
        if (r < cfg[1]) {
            *n = (int)cfg[2] + (int)floor(rand_u(p) * (cfg[3] - cfg[2] + 1.0));
        } else if (r < cfg[1] + cfg[0]) {
            *n = 1;
        }
    }
    bool bad = *n > 0;
    if (bad) (*n)--;
    return bad;
}

/* Ground uplink: latch truth, slide it along its orbit by an OD timing error, fit SGP4 mean elements */
static int do_uplink(plant_t* p, const err_model_t* e, plant_body_inputs_t* in) {
    double dt_err = e->sigma_at_s * randn(p);
    if (p->n_uplinks == 0 && p->cfg.worst_case == 0) dt_err = 0.0; /* nominal: first TLE exact */
    double a[3], r0[3], v0[3];
    orbit_accel(p->r, a);
    for (int i = 0; i < 3; i++) {
        r0[i] = p->r[i] + p->v[i] * dt_err + 0.5 * a[i] * dt_err * dt_err;
        v0[i] = p->v[i] + a[i] * dt_err;
    }
    double jd = jd_now(p);
    float tle[6];
    float bstar = 0.0f; /* drag off in the plant */
    double resid = adcs_tle_fit(r0, v0, jd, bstar, tle);
    if (!(resid >= 0.0)) return -1;
    memcpy(in->posn_update, tle, sizeof(tle));
    in->bstar = bstar;
    in->tle_epoch_jd = jd;
    in->has_posn_update = true;
    p->tle_fit_resid_m = resid;
    p->n_uplinks++;
    p->next_uplink_jd = jd + e->uplink_period_s / 86400.0;
    return 0;
}

/* ---------------- public API ---------------- */
void plant_default_config(plant_config_t* cfg) {
    memset(cfg, 0, sizeof(*cfg));
    cfg->q0_body_to_eci[0] = 1.0;
    cfg->seed = 1;
    cfg->nan_dropouts = true;
    cfg->gravity_gradient = true;
    cfg->kick_on = false;
    cfg->kick_t_s = 20000.0;
    cfg->kick_len_s = 1.0;
    cfg->kick_dw[0] = 0.20;
    cfg->kick_dw[1] = -0.15;
    cfg->kick_dw[2] = 0.15;
}

int plant_init(plant_t* p, const plant_config_t* cfg) {
    memset(p, 0, sizeof(*p));
    p->cfg = *cfg;
    rng_seed(p, cfg->seed);
    err_model_t e;
    error_model(cfg->worst_case, &e);

    /* Truth orbit at t = 0 from the TLE (flight SGP4, then integrated in double from here on) */
    sgp4_sat_t sat;
    sgp4_init(cfg->tle_mean, cfg->bstar, cfg->tle_epoch_jd, &sat);
    if (sat.initErr != 0) return -1;
    float rg[3], vg[3], osc[7], mean[6], rt[3], vt[3];
    float dt0 = (float)((cfg->start_jd - cfg->tle_epoch_jd) * 86400.0);
    if (sgp4_propagate(&sat, dt0, rg, vg, osc, mean, rt, vt) != 0) return -1;
    for (int i = 0; i < 3; i++) {
        p->r[i] = rg[i] * 1000.0;
        p->v[i] = vg[i] * 1000.0;
    }

    memcpy(p->q, cfg->q0_body_to_eci, sizeof(p->q));
    normalize_q(p->q);
    memcpy(p->w, cfg->w0_body_rad_s, sizeof(p->w));
    memcpy(p->gyro_bias, e.gyro_bias0, sizeof(p->gyro_bias));
    for (int i = 0; i < NUM_DIODES; i++) p->pd_gain[i] = 1.0 + e.pd_spread * (2.0 * rand_u(p) - 1.0);
    p->next_uplink_jd = cfg->start_jd; /* uplink on the first tick */
    return 0;
}

void plant_sense(plant_t* p, plant_body_inputs_t* in) {
    err_model_t e;
    error_model(p->cfg.worst_case, &e);
    memset(in, 0, sizeof(*in));

    /* ---- clock (Orbit + Time block) ---- */
    double jd = jd_now(p);
    double jd_int = floor(jd);
    in->jd_scalar = (int)jd_int;
    in->jd_frac = (float)(jd - jd_int);
    in->unix_time = unix_time_of(jd);
    in->dt = (float)PLANT_DT;

    /* ---- TLE uplink ---- */
    if (p->n_uplinks == 0 || jd >= p->next_uplink_jd) {
        if (do_uplink(p, &e, in) != 0) in->has_posn_update = false;
    }

    /* ---- magnetometer: body frame, T, bias + white noise ---- */
    double b_eci[3], b_body[3];
    field_eci(p, p->r, b_eci);
    eci_to_body(p->q, b_eci, b_body);
    double mag[3];
    for (int i = 0; i < 3; i++) mag[i] = b_body[i] + e.mag_bias[i] + e.sigma_mag * randn(p);

    /* ---- gyro: bias random walk + white noise ---- */
    double gyro[3];
    for (int i = 0; i < 3; i++) {
        p->gyro_bias[i] += e.sigma_rrw * sqrt(PLANT_DT) * randn(p);
        gyro[i] = p->w[i] + p->gyro_bias[i] + e.sigma_gyro[i] * randn(p);
    }

    /* ---- photodiodes: cosine law, eclipse, albedo, gain spread, noise ---- */
    double s[3], s_body[3], r_hat[3], nadir_body[3];
    sun_eci(p, s);
    eci_to_body(p->q, s, s_body);
    double rn = norm3(p->r);
    for (int i = 0; i < 3; i++) r_hat[i] = p->r[i] / rn;
    double m_rhat[3] = {-r_hat[0], -r_hat[1], -r_hat[2]};
    eci_to_body(p->q, m_rhat, nadir_body);
    double lit = sunlit(p->r, s) ? 1.0 : 0.0;
    double alb_scale = e.albedo * (R_E / rn) * (R_E / rn) * fmax(0.0, dot(s, r_hat));
    double pd[NUM_DIODES];
    for (int k = 0; k < NUM_DIODES; k++) {
        double nk[3] = {PHOTODIODES[k][0], PHOTODIODES[k][1], PHOTODIODES[k][2]};
        double direct = fmax(0.0, dot(nk, s_body)) * lit;
        double alb = alb_scale * fmax(0.0, dot(nk, nadir_body));
        pd[k] = PD_MAX * p->pd_gain[k] * (direct + alb) + e.sigma_pd * randn(p);
        if (pd[k] < 0.0) pd[k] = 0.0;
    }

    /* ---- NaN dropouts ---- */
    bool mag_bad = false, gyro_bad = false, pd_bad[NUM_DIODES] = {false};
    if (p->cfg.nan_dropouts) {
        mag_bad = outage(p, &p->out_mag, MAG_NAN);
        gyro_bad = outage(p, &p->out_gyro, GYRO_NAN);
        for (int k = 0; k < NUM_DIODES; k++) pd_bad[k] = outage(p, &p->out_pd[k], PD_NAN);
    }
    for (int i = 0; i < 3; i++) {
        in->magnetometer[i] = mag_bad ? NAN : (float)mag[i];
        in->gyro[i] = gyro_bad ? NAN : (float)gyro[i];
    }
    for (int k = 0; k < NUM_DIODES; k++) in->photodiodes[k] = pd_bad[k] ? NAN : (float)pd[k];

    /* last_mag = previous tick's sample (carries a dropped NaN forward, like the flight buffer) */
    if (!p->have_prev_mag) {
        memcpy(p->prev_mag, in->magnetometer, sizeof(p->prev_mag));
        p->have_prev_mag = true;
    }
    memcpy(in->last_magnetometer, p->prev_mag, sizeof(p->prev_mag));
    memcpy(p->prev_mag, in->magnetometer, sizeof(p->prev_mag));
}

void plant_step(plant_t* p, const float* currents) {
    err_model_t e;
    error_model(p->cfg.worst_case, &e);

    /* ---- actuators: currents -> dipole (scale factor, residual dipole, noise) -> torque ---- */
    double b_eci[3], b_body[3], dip[3], tq[3];
    field_eci(p, p->r, b_eci);
    eci_to_body(p->q, b_eci, b_body);
    for (int i = 0; i < 3; i++) {
        double I = (currents != NULL && isfinite(currents[i])) ? currents[i] : 0.0;
        dip[i] = (1.0 + e.sf_err[i]) * K_COIL[i] * I + e.m_res[i] + e.sigma_D * randn(p);
    }
    crossd(dip, b_body, tq);
    if (p->cfg.kick_on && p->t >= p->cfg.kick_t_s && p->t < p->cfg.kick_t_s + p->cfg.kick_len_s) {
        double Jdw[3];
        matvec(J_SC, p->cfg.kick_dw, Jdw);
        for (int i = 0; i < 3; i++) tq[i] += Jdw[i] / p->cfg.kick_len_s;
    }

    /* ---- attitude RK4 (torque held over the tick), 2 substeps ---- */
    double rn = norm3(p->r);
    double r_hat[3] = {p->r[0] / rn, p->r[1] / rn, p->r[2] / rn};
    double n2 = MU_E / (rn * rn * rn);
    const int nsub = 2;
    double h = PLANT_DT / nsub;
    for (int s = 0; s < nsub; s++) {
        double k1q[4], k1w[3], k2q[4], k2w[3], k3q[4], k3w[3], k4q[4], k4w[3], qt[4], wt[3];
        att_deriv(p->q, p->w, tq, r_hat, n2, p->cfg.gravity_gradient, k1q, k1w);
        for (int i = 0; i < 4; i++) qt[i] = p->q[i] + 0.5 * h * k1q[i];
        for (int i = 0; i < 3; i++) wt[i] = p->w[i] + 0.5 * h * k1w[i];
        att_deriv(qt, wt, tq, r_hat, n2, p->cfg.gravity_gradient, k2q, k2w);
        for (int i = 0; i < 4; i++) qt[i] = p->q[i] + 0.5 * h * k2q[i];
        for (int i = 0; i < 3; i++) wt[i] = p->w[i] + 0.5 * h * k2w[i];
        att_deriv(qt, wt, tq, r_hat, n2, p->cfg.gravity_gradient, k3q, k3w);
        for (int i = 0; i < 4; i++) qt[i] = p->q[i] + h * k3q[i];
        for (int i = 0; i < 3; i++) wt[i] = p->w[i] + h * k3w[i];
        att_deriv(qt, wt, tq, r_hat, n2, p->cfg.gravity_gradient, k4q, k4w);
        for (int i = 0; i < 4; i++) p->q[i] += h / 6.0 * (k1q[i] + 2 * k2q[i] + 2 * k3q[i] + k4q[i]);
        for (int i = 0; i < 3; i++) p->w[i] += h / 6.0 * (k1w[i] + 2 * k2w[i] + 2 * k3w[i] + k4w[i]);
        normalize_q(p->q);
    }

    /* record the total applied torque (incl. gravity gradient at the start of the tick) */
    if (p->cfg.gravity_gradient) {
        double rb[3], Jr[3], tg[3];
        eci_to_body(p->q, r_hat, rb);
        matvec(J_SC, rb, Jr);
        crossd(rb, Jr, tg);
        for (int i = 0; i < 3; i++) tq[i] += 3.0 * n2 * tg[i];
    }
    memcpy(p->torque, tq, sizeof(tq));

    /* ---- orbit ---- */
    orbit_rk4(p->r, p->v, PLANT_DT);
    p->t += PLANT_DT;
}

void plant_get_truth(const plant_t* p, plant_truth_t* truth) {
    memset(truth, 0, sizeof(*truth));
    truth->t_s = p->t;
    truth->jd = jd_now(p);
    memcpy(truth->r_eci_m, p->r, sizeof(p->r));
    memcpy(truth->v_eci_mps, p->v, sizeof(p->v));
    memcpy(truth->q_body_to_eci, p->q, sizeof(p->q));
    memcpy(truth->w_body_rad_s, p->w, sizeof(p->w));
    field_eci(p, p->r, truth->b_eci_T);
    sun_eci(p, truth->sun_eci);
    truth->sunlit = sunlit(p->r, truth->sun_eci);
    memcpy(truth->torque_body_Nm, p->torque, sizeof(p->torque));
    memcpy(truth->gyro_bias_rad_s, p->gyro_bias, sizeof(p->gyro_bias));
    truth->tle_fit_resid_m = p->tle_fit_resid_m;
}

void plant_call_body(const plant_body_inputs_t* in, float* currents) {
    if (in->has_posn_update) body_set_tle_metadata(in->bstar, in->tle_epoch_jd);
    body(in->last_magnetometer, in->magnetometer, in->gyro, in->photodiodes,
         in->has_posn_update ? in->posn_update : NULLPTR, in->unix_time, in->jd_scalar, in->jd_frac,
         in->dt, currents);
}

int plant_fsw_step_host(const plant_body_inputs_t* in, float* currents, void* ctx) {
    (void)ctx;
    plant_call_body(in, currents);
    return 0;
}

long plant_run(plant_t* p, const plant_config_t* cfg, long n_ticks, plant_fsw_step_fn fsw_step,
               void* fsw_ctx, plant_tick_fn on_tick, void* tick_ctx) {
    if (plant_init(p, cfg) != 0) return -1;
    if (fsw_step == NULL) fsw_step = plant_fsw_step_host;
    for (long k = 0; k < n_ticks; k++) {
        plant_body_inputs_t in;
        float currents[3] = {0.0f, 0.0f, 0.0f};
        plant_sense(p, &in);                              /* 1. truth -> sensors / clock / TLE */
        if (fsw_step(&in, currents, fsw_ctx) != 0) return k; /* 2. flight software -> currents */
        plant_step(p, currents);                          /* 3. currents -> torque -> next state */
        if (on_tick != NULL) {
            plant_truth_t truth;
            plant_get_truth(p, &truth);
            on_tick(k, &in, currents, &truth, tick_ctx);
        }
    }
    return n_ticks;
}
