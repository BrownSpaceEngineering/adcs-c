#ifndef PLANT_H
#define PLANT_H
/*
 * Truth plant for PIL / host closed-loop testing of body() (NOT flight code).
 *
 * C port of the Simulink C-ADCS harness (asbCubeSat/Vehicle Model/C ADCS):
 *   - orbit:      J2 two-body integrated in double, initialised from the given TLE with the flight
 *                 SGP4; ground "uplinks" refit SGP4 mean elements to the truth state
 *                 (port of the "Orbit + Time (ground TLE)" block, fit = tle_fit.c)
 *   - sensors:    port of the "Sensors (noise)" block (magnetometer, gyro, 18 photodiodes,
 *                 eclipse, albedo, NaN dropouts)
 *   - actuators:  port of the "ADCS_Actuators" block (coil currents -> dipole with errors)
 *   - dynamics:   rigid body, Spacecraft Dynamics inertia, gravity gradient, RK4
 *   - field/sun:  the same WMM / sun models the flight code uses
 *
 * One control tick (PLANT_DT = 0.1 s):
 *     plant_sense(&p, &in);                  // sensors at the current truth state
 *     if (in.has_posn_update) body_set_tle_metadata(in.bstar, in.tle_epoch_jd);
 *     body(in.last_magnetometer, in.magnetometer, in.gyro, in.photodiodes,
 *          in.has_posn_update ? in.posn_update : NULLPTR,
 *          in.unix_time, in.jd_scalar, in.jd_frac, in.dt, currents);
 *     plant_step(&p, currents);              // apply currents, integrate to the next tick
 * plant_call_body() does the middle two lines for a host run. In PIL, send `in` to the target
 * and the 3 currents back.
 *
 * Frames / units match body.h: magnetometer body frame in TESLA, gyro rad/s, photodiodes raw,
 * posn_update = TLE mean [n (rev/day), e, i, RAAN, argp, M] with angles in radians.
 * Quaternions are Hamilton scalar-first, ACTIVE body->ECI (same as the filter estimate).
 */
#include <stdbool.h>
#include <stdint.h>
#include "include/photodiode_determination.h"

#define PLANT_DT 0.1 /* s, FSW tick (body() rate) */

typedef struct {
    /* Initial orbit: TLE mean elements in body() posn_update order, B* and epoch (UTC JD) */
    float tle_mean[6];
    float bstar;
    double tle_epoch_jd;
    double start_jd;                /* UTC JD at t = 0; may be after the TLE epoch */

    /* Initial truth attitude / rate */
    double q0_body_to_eci[4];       /* normalised internally */
    double w0_body_rad_s[3];        /* body rate w.r.t. ECI, body axes */

    int worst_case;                 /* 0 nominal, 1 worst case (adcs_worst_case in Simulink) */
    uint64_t seed;                  /* RNG seed (runs are reproducible) */
    bool nan_dropouts;              /* sensor NaN fault injection (Simulink default: on) */
    bool gravity_gradient;          /* Spacecraft Dynamics useGravGrad (Simulink: on) */

    /* Tumble kick (ADCS_Actuators): adds J * kick_dw / kick_len torque for kick_len s */
    bool kick_on;
    double kick_t_s, kick_len_s, kick_dw[3];
} plant_config_t;

/* Exactly the arguments body() takes for one tick (plus the TLE metadata call). */
typedef struct {
    float last_magnetometer[3];     /* T, body (previous tick's sample; NaN during dropouts) */
    float magnetometer[3];          /* T, body */
    float gyro[3];                  /* rad/s, body */
    float photodiodes[NUM_DIODES];  /* raw, PHOTODIODES order */
    bool has_posn_update;           /* true on an uplink tick: pass posn_update, else NULLPTR */
    float posn_update[6];
    float bstar;                    /* call body_set_tle_metadata(bstar, tle_epoch_jd) first */
    double tle_epoch_jd;
    int unix_time;
    int jd_scalar;
    float jd_frac;
    float dt;
} plant_body_inputs_t;

/* Truth, for logging / scoring */
typedef struct {
    double t_s;                     /* s since start */
    double jd;                      /* UTC */
    double r_eci_m[3], v_eci_mps[3];
    double q_body_to_eci[4];
    double w_body_rad_s[3];
    double b_eci_T[3];
    double sun_eci[3];
    bool sunlit;
    double torque_body_Nm[3];       /* last applied (magnetorquer + gravity gradient + kick) */
    double gyro_bias_rad_s[3];      /* true gyro bias */
    double tle_fit_resid_m;         /* residual of the last uplinked TLE fit at its epoch */
} plant_truth_t;

typedef struct {
    plant_config_t cfg;
    double t;                       /* s since start */
    double r[3], v[3];              /* m, m/s, GCRF */
    double q[4], w[3];              /* body->ECI, rad/s */
    double torque[3];
    /* sensor state */
    double gyro_bias[3];
    float prev_mag[3];
    bool have_prev_mag;
    int out_mag, out_gyro, out_pd[NUM_DIODES];
    double pd_gain[NUM_DIODES];
    /* uplink state */
    double next_uplink_jd;
    int n_uplinks;
    double tle_fit_resid_m;
    /* RNG */
    uint64_t rng[4];
    bool have_spare_normal;
    double spare_normal;
} plant_t;

/* Fills cfg with the Simulink defaults (nominal, dropouts on, GG on, kick off, rest at identity) */
void plant_default_config(plant_config_t* cfg);
/* Returns 0 on success, -1 if the TLE can't be initialised / propagated to start_jd */
int plant_init(plant_t* p, const plant_config_t* cfg);
/* Sample the sensors / clock / uplink at the current truth state -> body() inputs */
void plant_sense(plant_t* p, plant_body_inputs_t* in);
/* Apply coil currents (A) for one tick and integrate the truth state forward PLANT_DT */
void plant_step(plant_t* p, const float* currents);
void plant_get_truth(const plant_t* p, plant_truth_t* truth);
/* Host convenience: metadata + body() call for one tick */
void plant_call_body(const plant_body_inputs_t* in, float* currents);

/*
 * One flight-software step: take this tick's body() inputs, write the 3 coil currents (A).
 * Return 0 on success, nonzero to abort the run (e.g. a PIL link timeout).
 *   host:  fsw_step = plant_fsw_step_host (calls body() directly)
 *   PIL:   your function that sends `in` to the target and reads the currents back
 */
typedef int (*plant_fsw_step_fn)(const plant_body_inputs_t* in, float* currents, void* ctx);
/* Called after every tick with the truth state and the currents that were applied (may be NULL) */
typedef void (*plant_tick_fn)(long tick, const plant_body_inputs_t* in, const float* currents,
                              const plant_truth_t* truth, void* ctx);

/* plant_fsw_step_fn wrapper around plant_call_body (ctx unused) */
int plant_fsw_step_host(const plant_body_inputs_t* in, float* currents, void* ctx);

/*
 * Closed loop: plant_init, then n_ticks x (plant_sense -> fsw_step -> plant_step -> on_tick).
 * Does NOT call body_reset(); do that first for a host run (on PIL, reset the target).
 * Returns the number of ticks completed, or -1 if plant_init failed.
 */
long plant_run(plant_t* p, const plant_config_t* cfg, long n_ticks, plant_fsw_step_fn fsw_step,
               void* fsw_ctx, plant_tick_fn on_tick, void* tick_ctx);
#endif
