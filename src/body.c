#include "include/body.h"
#include "arm_math.h"
#include "include/bdot.h"
#include "include/pointing_error.h"
#include "include/ecef2eci.h"
#include "include/iterate.h"
#include "include/laextension.h"
#include "include/magnetosphere.h"
#include "include/moments2currents.h"
#include "include/pd.h"
#include "include/newpointing.h"
#include "include/photodiode_determination.h"
#include "include/quat.h"
#include "include/quest.h"
#include "include/sensor_frames.h"
#include "include/sgp4_step.h"
#include "include/sunvec.h"
#include "include/torque2moments.h"
#include "math.h"
#include "string.h"
#include <limits.h>
#include <stdbool.h>




// Magnetometer [U8, U9, U10] and gyro readings arrive in sensor axes; body() maps them to the body
// frame (sensor_frames.h) before anything else sees them.

// Providence, lla2ecef(41.825226, -71.418884, 0) (see eceftoeci.m)
float PVD_ECEF[3] = {6.3761e6f, -0.1387e6f, 0.0807e6f};
const float DETUMBLING_GYRO_THRESHOLD = 0.005; // rad/s (~4-5x orbit rate; pointing failed to converge from 0.014)
const float DETUMBLING_MAG_THRESHOLD = 4e-7f;  // |dB/dt|, T/s. TODO: tune (anything in T/s passes)
// B-dot gain for m = -k dB/dt (A m^2 per T/s). Sized from the optimal-gain result
// k* = 2 w_orb (1 + sin xi_m) J_min / |B|^2 (Avanzini & Giulietti 2012) with w_orb = 1.07e-3 rad/s,
// xi_m ~ 60 deg, J_min = 0.0124 kg m^2, |B| ~ 30 uT -> ~5e4. That optimum assumes continuous,
// unsaturated control; with the 1 s clean-sample B-dot cycle (80% duty) and the 0.156 A coil limit,
// a host Monte Carlo (plant.c, 8-10 seeds per case) handed off much sooner with 4x that:
//   time to hand-off, median (5e4 -> 2e5): 3 deg/s 4223 -> 1198 s, 17 deg/s 7056 -> 4927 s,
//   30 deg/s 8671 -> 3998 s; no resets. 2e6 was much worse (saturated chatter, ~19000 s).
// Saturates the coils above ~0.075 rad/s.
#ifndef BODY_BDOT_GAIN
#define BODY_BDOT_GAIN 2e5f
#endif
const float BDOT_GAIN = BODY_BDOT_GAIN;
const float* Imax = NULLPTR;                     // NULLPTR uses hardware defaults (see moments2currents.c)

// Filter tuning. Measurement vectors are unit vectors, so R is in (unit vector)^2.
// Initial bias variance and process noise follow simulink.m (Q there is per 0.1 s step).
// clang-format off
const float INIT_P[6 * 6] = {
    0.01f, 0,     0,     0,       0,       0,
    0,     0.01f, 0,     0,       0,       0,
    0,     0,     0.01f, 0,       0,       0,
    0,     0,     0,     3.05e-6f, 0,       0,
    0,     0,     0,     0,       3.05e-6f, 0,
    0,     0,     0,     0,       0,       3.05e-6f};
// Process noise per second (multiplied by dt each step). Spectral densities for the
// Murata SCH16T-K01, inflated to cover unmodelled errors (field-model mismatch, albedo,
// scale factor/misalignment, thermal bias drift) and to keep P well above the float32 /
// SIGMA_JITTER floor:
//   attitude: ARW 0.4 mdps/rtHz -> sigma_v^2 = 5e-11 rad^2/s; use 1e-8 (tune 1e-8..1e-7)
//   bias:     BI 0.5 deg/h -> sigma_u^2 ~ 2e-14 rad^2/s^3;      use 1e-12 (tune 1e-13..1e-11)
const float Q_RATE[6 * 6] = {
    1e-8f, 0,     0,     0,      0,      0,
    0,     1e-8f, 0,     0,      0,      0,
    0,     0,     1e-8f, 0,      0,      0,
    0,     0,     0,     1e-12f, 0,      0,
    0,     0,     0,     0,      1e-12f, 0,
    0,     0,     0,     0,      0,      1e-12f};
// [magnetometer (3), photodiode sun vector (3)]
const float R_IN_SUN[6 * 6] = {
    4e-4f, 0,     0,     0,       0,       0,
    0,     4e-4f, 0,     0,       0,       0,
    0,     0,     4e-4f, 0,       0,       0,
    0,     0,     0,     2.5e-3f, 0,       0,
    0,     0,     0,     0,       2.5e-3f, 0,
    0,     0,     0,     0,       0,       2.5e-3f};
// Magnetometer only (eclipse), like mode 1 of simulink.m
const float R_IN_SHADOW[3 * 3] = {
    4e-4f, 0,     0,
    0,     4e-4f, 0,
    0,     0,     4e-4f};
// clang-format on

// Reject filter output that jumps further than this from the gyro-propagated estimate
#define MAX_FILTER_JUMP_RAD 0.35f
#define MAX_ATTITUDE_VARIANCE 10.0f
// Declare divergence if a measured vector disagrees with the estimate by more than this for
// this many consecutive steps
#define MAX_INNOVATION_RAD 0.45f
#define MAX_INNOVATION_STRIKES 5
// A rate above this while pointing indicates a new tumble or a gyro fault.
#define MAX_POINTING_RATE_RAD_S 0.3f
#define ORBIT_MAX_FAIL_TICK 100u

// ---- Module state (previously `static` in body.h, which gave every includer its own copy) ----
static float tle_mean[6];
static float tle_bstar;
static double tle_epoch_jd;
static float pending_bstar;
static float last_r_eci[3];
static bool position_valid = false;
static double pending_epoch_jd;
static bool tle_metadata_pending = false;
static float estimated_quat[4];
static float estimated_gyro_bias[3];
static float error_quat_state[6];
static float error_quat_cov[6 * 6];
static float q_want_prev[4] = {1, 0, 0, 0};

static bool pointing = false;
static bool tle_initialized = false;
static bool attitude_initialized = false;
static int filter_resets = 0;
static int innovation_strikes = 0;
static unsigned int pointing_tick = 0;
static unsigned int mag_fault_ticks = 0; // consecutive ticks with an unusable magnetometer sample
static unsigned int ready_ticks = 0;      // consecutive ticks the detumble -> pointing check passed
static float dbdt_filt[3];                // low-passed dB/dt (T/s) for the hand-off check
static bool dbdt_filt_valid = false;
// B-dot cycle (see BODY_BDOT_PERIOD_TICKS): one clean magnetometer sample per cycle
static unsigned int bdot_tick = 0;        // 1..BODY_BDOT_PERIOD_TICKS, sample on the last tick
static float bdot_prev_mag[3];            // last clean sample (T, body)
static bool bdot_prev_valid = false;
static float bdot_elapsed_s = 0.0f;       // time since bdot_prev_mag was taken
static float bdot_moment[3];              // dipole held between samples (A m^2)
// Sensor sanity state (stuck = bit-identical repeats; see BODY_SENSOR_STUCK_TICKS)
static float gyro_last_raw[3], mag_last_raw[3], pd_last_raw[NUM_DIODES];
static unsigned int gyro_repeat_ticks = 0, mag_repeat_ticks = 0, pd_repeat_ticks[NUM_DIODES];
static float gyro_ref[3];                 // last accepted gyro sample (spike check)
static bool gyro_ref_valid = false;
static unsigned int gyro_step_ticks = 0;  // consecutive samples far from gyro_ref
static body_sensor_rejects_t rejects;
static unsigned int gate_streak = 0;     // consecutive measurement updates with a gated-out vector
// Magnetometer hard-iron calibration (RLS, in uT)
static float magcal_b[3];                 // uT
static float magcal_P[9];                 // uT^2
static unsigned int magcal_n = 0, magcal_ticks = 0;
// Magnetorquer configuration, applied to every command (set from the ground via body_set_coil_config)
static float coil_sign[3] = {1.0f, 1.0f, 1.0f};
static bool coil_enabled[3] = {true, true, true};
// Time (s) since the filter was last propagated. Accumulates every pointing tick and is consumed by
// the next UKF step, so ticks skipped by an early return (bad gyro, etc.) aren't lost.
static float unpropagated_s = 0.0f;

static const float ZERO3[3] = {0, 0, 0};

static void reset_bdot_cycle(void) {
    bdot_tick = 0;
    bdot_prev_valid = false;
    bdot_elapsed_s = 0.0f;
    memset(bdot_prev_mag, 0, sizeof(bdot_prev_mag));
    memset(bdot_moment, 0, sizeof(bdot_moment));
}

void body_reset(void) {
    memset(tle_mean, 0, sizeof(tle_mean));
    tle_bstar = 0.0;
    tle_epoch_jd = 0.0;
    pending_bstar = 0.0;
    pending_epoch_jd = 0.0;
    tle_metadata_pending = false;
    memset(estimated_quat, 0, sizeof(estimated_quat));
    estimated_quat[0] = 1.0f;
    memset(estimated_gyro_bias, 0, sizeof(estimated_gyro_bias));
    memset(error_quat_state, 0, sizeof(error_quat_state));
    memset(error_quat_cov, 0, sizeof(error_quat_cov));
    q_want_prev[0] = 1.0f;
    q_want_prev[1] = 0.0f;
    q_want_prev[2] = 0.0f;
    q_want_prev[3] = 0.0f;
    pointing = false;
    tle_initialized = false;
    position_valid = false;
    memset(last_r_eci, 0, sizeof(last_r_eci));
    attitude_initialized = false;
    filter_resets = 0;
    innovation_strikes = 0;
    pointing_tick = 0;
    mag_fault_ticks = 0;
    unpropagated_s = 0.0f;
    ready_ticks = 0;
    dbdt_filt_valid = false;
    reset_bdot_cycle();
    // Sensor sanity state belongs to the sensors, so only a full reset clears it
    gyro_repeat_ticks = 0;
    mag_repeat_ticks = 0;
    memset(pd_repeat_ticks, 0, sizeof(pd_repeat_ticks));
    memset(gyro_last_raw, 0, sizeof(gyro_last_raw));
    memset(mag_last_raw, 0, sizeof(mag_last_raw));
    memset(pd_last_raw, 0, sizeof(pd_last_raw));
    gyro_ref_valid = false;
    gyro_step_ticks = 0;
    gate_streak = 0;
    memset(&rejects, 0, sizeof(rejects));
    memset(magcal_b, 0, sizeof(magcal_b));
    memset(magcal_P, 0, sizeof(magcal_P));
    for (int i = 0; i < 3; i++) magcal_P[4 * i] = BODY_MAGCAL_P0_UT2;
    magcal_n = 0;
    magcal_ticks = BODY_MAGCAL_MIN_TICKS;
    for (int i = 0; i < 3; i++) {
        coil_sign[i] = 1.0f;
        coil_enabled[i] = true;
    }
}

void body_get_sensor_rejects(body_sensor_rejects_t* out) { *out = rejects; }

bool body_set_tle_metadata(float bstar, double epoch_jd) {
    if (!isfinite(bstar) || !isfinite(epoch_jd)) {
        return false;
    }
    pending_bstar = bstar;
    pending_epoch_jd = epoch_jd;
    tle_metadata_pending = true;
    return true;
}

bool body_get_tle(float* tle_mean_out, float* bstar, double* epoch_jd) {
    if (!tle_initialized) {
        return false;
    }
    memcpy(tle_mean_out, tle_mean, sizeof(tle_mean));
    *bstar = tle_bstar;
    *epoch_jd = tle_epoch_jd;
    return true;
}

bool body_get_position(float* r_eci_m) {
    if (!position_valid) {
        return false;
    }
    memcpy(r_eci_m, last_r_eci, sizeof(last_r_eci));
    return true;
}

bool body_is_pointing(void) { return pointing; }

bool body_get_attitude(float* q_body_to_eci) {
    memcpy(q_body_to_eci, estimated_quat, sizeof(estimated_quat));
    return attitude_initialized;
}

void body_get_gyro_bias(float* bias) { memcpy(bias, estimated_gyro_bias, sizeof(estimated_gyro_bias)); }

int body_get_filter_resets(void) { return filter_resets; }

static void zero_currents(float* output_currents) { memcpy(output_currents, ZERO3, sizeof(ZERO3)); }

static void command_zero_moment(float* output_currents) {
    moment2current3axis(ZERO3, Imax, output_currents);
}

static void reset_filter_to_detumble(void) {
    memset(estimated_quat, 0, sizeof(estimated_quat));
    estimated_quat[0] = 1.0f;
    memset(estimated_gyro_bias, 0, sizeof(estimated_gyro_bias));
    memset(error_quat_state, 0, sizeof(error_quat_state));
    memset(error_quat_cov, 0, sizeof(error_quat_cov));
    attitude_initialized = false;
    pointing = false;
    innovation_strikes = 0;
    gate_streak = 0;
    pointing_tick = 0;
    mag_fault_ticks = 0;
    unpropagated_s = 0.0f;
    ready_ticks = 0;
    dbdt_filt_valid = false;
    reset_bdot_cycle();
    filter_resets++;
}

// Largest angle between a measured body vector and the one predicted by q_body_to_ref
static float worst_innovation(const float* q_body_to_ref, const float* body_vecs, const float* ref_vecs,
                              int num_vecs) {
    float q_ref_to_body[4];
    quat_conj(q_body_to_ref, q_ref_to_body);
    float worst = 0.0f;
    for (int v = 0; v < num_vecs; v++) {
        float predicted[3], c[3];
        quat_apply(q_ref_to_body, &ref_vecs[3 * v], predicted);
        cross(&body_vecs[3 * v], predicted, c);
        float angle = atan2f(l2_norm(c, 3), dot3(&body_vecs[3 * v], predicted));
        worst = fmaxf(worst, angle);
    }
    return worst;
}

/**
 * Declares the filter failed if its output is non-finite, its covariance is unhealthy, or its
 * estimate jumped far from where the (bias-corrected) gyro says it should be.
 *
 * The previous version checked the OLD covariance for NaNs, and compared |dq/dt| / |gyro| > 10,
 * which divides by ~0 (and always "fails") whenever the satellite is nearly still.
 */
static bool filter_failure(const float* new_quat, const float* gyro, const float* new_P, float dt) {
    if (!all_finite(new_quat, 4) || !all_finite(new_P, 36)) {
        return true;
    }
    for (int i = 0; i < 6; i++) {
        if (!(new_P[i * 6 + i] > 0.0f)) {
            return true;
        }
    }
    for (int i = 0; i < 3; i++) {
        if (new_P[i * 6 + i] > MAX_ATTITUDE_VARIANCE) {
            return true;
        }
    }

    float omega_dt[3];
    for (int i = 0; i < 3; i++) {
        omega_dt[i] = (gyro[i] - estimated_gyro_bias[i]) * dt;
    }
    float dq[4], predicted[4], diff[4], diff_rot[3];
    rotationvec2quat(omega_dt, dq);
    quat_multiply(estimated_quat, dq, predicted);
    quat_diff(predicted, new_quat, diff);
    quat2rotationvec(diff, diff_rot);
    return l2_norm(diff_rot, 3) > MAX_FILTER_JUMP_RAD;
}

static bool mag_ok(const float* measurement) {
    if (measurement == NULLPTR) {
        rejects.null_ptr++;
        return false;
    }
    if (!all_finite(measurement, 3)) {
        rejects.mag_nan++;
        return false;
    }
    float norm = l2_norm(measurement, 3);
    if (!(norm > 1e-30f) || !isfinite(norm)) {
        rejects.mag_norm++; // all-zero (or overflowing) vector: an impossible field magnitude
        return false;
    }
    return true;
}

/**
 * True when it's safe to stop detumbling: the gyro rate AND the low-passed |dB/dt| are both under
 * their thresholds, the sun is visible (QUEST needs it), and all of that has held for
 * BODY_READY_TO_POINT_TICKS consecutive ticks. Any failing tick resets the count.
 *
 * Why both, and why debounced: with an OR, a single noisy dB/dt sample under the threshold (or a
 * spin about the field line, where dB/dt ~ 0 at any rate) handed over to pointing at 0.25 rad/s.
 * The low-pass keeps magnetometer noise in dB/dt well under the threshold.
 *
 * dbdt_filt is updated once per B-dot cycle from clean (coils-off) samples by update_dbdt_filt();
 * this check runs every tick on the latest value.
 *
 * Fixes vs the original version: dB used last_mag[2] for the y and z components (and mag[1] for
 * z), and was multiplied by dt instead of divided.
 */
static void update_dbdt_filt(const float* dbdt) {
    if (!dbdt_filt_valid) {
        memcpy(dbdt_filt, dbdt, sizeof(dbdt_filt));
        dbdt_filt_valid = true;
    } else {
        for (int i = 0; i < 3; i++) {
            dbdt_filt[i] += BODY_DBDT_FILTER_ALPHA * (dbdt[i] - dbdt_filt[i]);
        }
    }
}

static bool ready_to_point(const float* photodiode_measurements, const float* gyro, bool gyro_valid) {
    if (gyro == NULLPTR || photodiode_measurements == NULLPTR || !dbdt_filt_valid) {
        ready_ticks = 0;
        return false;
    }
    float photodiode_sun_vec[3];
    bool in_sun = get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);
   

    bool calm = l2_norm(dbdt_filt, 3) < DETUMBLING_MAG_THRESHOLD
                && l2_norm(gyro, 3) < DETUMBLING_GYRO_THRESHOLD;

    if(!gyro_valid) {
        calm = false;
    }
    if (calm && in_sun) {
        if (ready_ticks < UINT_MAX) ready_ticks++;
    } else {
        ready_ticks = 0;
    }
    return ready_ticks >= BODY_READY_TO_POINT_TICKS;
}

static bool valid_tle_mean(const float* elements) {
    if (!isfinite(elements[0]) || elements[0] <= 0.0 || !isfinite(elements[1])
        || elements[1] < 0.0 || elements[1] >= 1.0) return false;
    for (int i = 2; i < 6; i++) {
        if (!isfinite(elements[i])) return false;
    }
    return true;
}

/**
 * Counts consecutive bit-identical samples. A working sensor always has noise, so N exact repeats
 * mean a stale driver buffer / frozen part, not a steady field (see BODY_SENSOR_STUCK_TICKS).
 * Returns true once the sample has repeated BODY_SENSOR_STUCK_TICKS times in a row.
 */
static bool repeat_check(const float* sample, float* last, int n, unsigned int* repeats) {
    if (memcmp(sample, last, sizeof(float) * (size_t)n) == 0) {
        if (*repeats < UINT_MAX) (*repeats)++;
    } else {
        *repeats = 0;
        memcpy(last, sample, sizeof(float) * (size_t)n);
    }
    return *repeats >= BODY_SENSOR_STUCK_TICKS;
}

/**
 * Gyro sample is usable: finite, not stuck, and not a one-tick jump. No torque on board can change
 * the rate by BODY_GYRO_MAX_STEP_RAD_S in one tick, so such a sample is a spike and is dropped
 * like a NaN; if the new value persists for BODY_GYRO_STEP_ACCEPT_TICKS it is real (a collision /
 * deployment kick) and is accepted, and the rate checks then act on it.
 */

static bool gyro_sample_ok(const float* gyro) {
    if (!all_finite(gyro, 3)) return false;
    if (repeat_check(gyro, gyro_last_raw, 3, &gyro_repeat_ticks)) {
        rejects.gyro_stuck++;
        return false;
    }
    if (gyro_ref_valid) {
        float d[3] = {gyro[0] - gyro_ref[0], gyro[1] - gyro_ref[1], gyro[2] - gyro_ref[2]};
        if (l2_norm(d, 3) > BODY_GYRO_MAX_STEP_RAD_S) {
            if (++gyro_step_ticks < BODY_GYRO_STEP_ACCEPT_TICKS) {
                rejects.gyro_spike++;
                return false;
            }
        }
    }
    memcpy(gyro_ref, gyro, sizeof(gyro_ref));
    gyro_ref_valid = true;
    gyro_step_ticks = 0;
    return true;
}

static bool magcal_ready(void) {
    if (magcal_n < 100u) return false;
    for (int i = 0; i < 3; i++) {
        if (!(magcal_P[4 * i] < BODY_MAGCAL_READY_UT * BODY_MAGCAL_READY_UT)) return false;
    }
    return l2_norm(magcal_b, 3) < BODY_MAGCAL_MAX_UT;
}

bool body_get_mag_bias(float* bias_T) {
    for (int i = 0; i < 3; i++) bias_T[i] = magcal_b[i] * 1e-6f;
    return magcal_ready();
}

/**
 * One RLS step of the attitude-independent hard-iron fit. With B_m = B + b:
 *   |B_m|^2 - |B_ref|^2 = 2 B_m.b - |b|^2   ->   y = |B_m|^2 - |B_ref|^2 + |b_hat|^2,  h = 2 B_m
 * Only call with a clean (coils-off) sample. Units: uT.
 */
static void magcal_update(const float* mag_T, float b_ref_norm_T) {
    float Bm[3] = {mag_T[0] * 1e6f, mag_T[1] * 1e6f, mag_T[2] * 1e6f};
    float bref = b_ref_norm_T * 1e6f;
    float y = dot3(Bm, Bm) - bref * bref + dot3(magcal_b, magcal_b);
    float h[3] = {2.0f * Bm[0], 2.0f * Bm[1], 2.0f * Bm[2]};
    float Ph[3];
    for (int i = 0; i < 3; i++) Ph[i] = magcal_P[3 * i] * h[0] + magcal_P[3 * i + 1] * h[1] + magcal_P[3 * i + 2] * h[2];
    float s = BODY_MAGCAL_FORGET * BODY_MAGCAL_R_UT4 + dot3(h, Ph);
    if (!(s > 0.0f) || !isfinite(s)) return;
    float innov = y - dot3(h, magcal_b);
    float b_new[3], P_new[9];
    for (int i = 0; i < 3; i++) b_new[i] = magcal_b[i] + Ph[i] / s * innov;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            P_new[3 * i + j] = (magcal_P[3 * i + j] - Ph[i] * Ph[j] / s) / BODY_MAGCAL_FORGET;
        }
    }
    if (!all_finite(b_new, 3) || !all_finite(P_new, 9)) return;
    memcpy(magcal_b, b_new, sizeof(magcal_b));
    memcpy(magcal_P, P_new, sizeof(magcal_P));
    // keep the prior as a ceiling so a long unobservable stretch can't blow P up
    for (int i = 0; i < 3; i++) {
        if (magcal_P[4 * i] > BODY_MAGCAL_P0_UT2) magcal_P[4 * i] = BODY_MAGCAL_P0_UT2;
    }
    if (magcal_n < UINT_MAX) magcal_n++;
    magcal_ticks = 0;
}

/* ---------------- magnetorquer configuration ---------------- */
void body_set_coil_config(const float* sign, const bool* enabled) {
    for (int i = 0; i < 3; i++) {
        if (sign != NULLPTR) coil_sign[i] = sign[i] < 0.0f ? -1.0f : 1.0f;
        if (enabled != NULLPTR) coil_enabled[i] = enabled[i];
    }
}

void body_get_coil_config(float* sign, bool* enabled) {
    for (int i = 0; i < 3; i++) {
        if (sign != NULLPTR) sign[i] = coil_sign[i];
        if (enabled != NULLPTR) enabled[i] = coil_enabled[i];
    }
}

/* Desired coil currents -> what is sent: polarity / enable from the coil configuration */
static void apply_coil_config(float* currents) {
    for (int i = 0; i < 3; i++) currents[i] = coil_enabled[i] ? coil_sign[i] * currents[i] : 0.0f;
}

/* With a coil disabled, re-solve the dipole on the remaining coils so they produce the closest
 * torque to what the 3-axis dipole m would have: min |sum_j u_j (e_j x b) - m x b| over enabled j. */
static void project_dipole_to_enabled(float* m, const float* b_body) {
    int idx[3], ne = 0;
    for (int i = 0; i < 3; i++) {
        if (coil_enabled[i]) idx[ne++] = i;
    }
    if (ne == 3) return;
    float tau[3];
    cross(m, b_body, tau);
    float u[3] = {0.0f, 0.0f, 0.0f};
    if (ne > 0) {
        float C[3][3]; // columns e_j x b
        for (int k = 0; k < ne; k++) {
            float e[3] = {0.0f, 0.0f, 0.0f};
            e[idx[k]] = 1.0f;
            float col[3];
            cross(e, b_body, col);
            for (int r = 0; r < 3; r++) C[r][k] = col[r];
        }
        // normal equations (ne <= 2)
        float N[4] = {0}, g[2] = {0};
        for (int a = 0; a < ne; a++) {
            for (int r = 0; r < 3; r++) g[a] += C[r][a] * tau[r];
            for (int b2 = 0; b2 < ne; b2++) {
                for (int r = 0; r < 3; r++) N[a * 2 + b2] += C[r][a] * C[r][b2];
            }
        }
        if (ne == 1) {
            if (N[0] > 0.0f) u[0] = g[0] / N[0];
        } else {
            float det = N[0] * N[3] - N[1] * N[2];
            if (fabsf(det) > 1e-30f) {
                u[0] = (N[3] * g[0] - N[1] * g[1]) / det;
                u[1] = (N[0] * g[1] - N[2] * g[0]) / det;
            }
        }
    }
    float out[3] = {0.0f, 0.0f, 0.0f};
    for (int k = 0; k < ne; k++) {
        int i = idx[k];
        out[i] = fmaxf(-NEW_PD_MAX_DIPOLE[i], fminf(u[k], NEW_PD_MAX_DIPOLE[i]));
    }
    memcpy(m, out, sizeof(out));
}

/* Magnetometer frozen (only checked on finite samples, so NaN dropouts don't count) */
static bool mag_stuck(const float* mag) {
    bool stuck = repeat_check(mag, mag_last_raw, 3, &mag_repeat_ticks);
    if (stuck) rejects.mag_stuck++;
    return stuck;
}

/**
 * Copies the photodiode readings, replacing implausible ones with NaN (which the sun-vector code
 * already treats as a dropout):
 *   - above BODY_PD_MAX_PLAUSIBLE_FRAC * MAX_READING or below -BODY_PD_NEG_TOL: railed / shorted
 *   - stuck above the dark level: exact repeats (a dark diode legitimately reads exactly 0 again
 *     and again, so the floor is excluded from the stuck check)
 */
static void sanitize_photodiodes(const float* in, float* out) {
    for (int k = 0; k < NUM_DIODES; k++) {
        float x = in[k];
        out[k] = x;
        if (!isfinite(x)){
            rejects.pd_rejected[k]++;
            rejects.pd_rejects_streak[k]++;
            continue;
        } 
        bool stuck = false;
        if (x > MIN_READING) {
            stuck = repeat_check(&x, &pd_last_raw[k], 1, &pd_repeat_ticks[k]);
        } else {
            pd_repeat_ticks[k] = 0;
            pd_last_raw[k] = x;
        }
        if (stuck || x > BODY_PD_MAX_PLAUSIBLE_FRAC * MAX_READING || x < -BODY_PD_NEG_TOL) {
            out[k] = NAN;
            rejects.pd_rejected[k]++;
            rejects.pd_rejects_streak[k]++;
        }
        else {
            rejects.pd_rejects_streak[k] = 0;
        }
    }
}

/**
 * Position for this tick from the saved TLE, or from posn_update if it is valid (finite, sane
 * elements, metadata from body_set_tle_metadata). An invalid update is ignored rather than
 * skipping the tick. Returns false if there is no orbit yet or propagation fails.
 */
static unsigned int orbit_err = BODY_ORBIT_OK; // why update_orbit last returned false

/* A new TLE failed: count it, record why, and use up its metadata (a retry must resend both) */
static void reject_new_tle(unsigned int reason) {
    rejects.tle_rejected++;
    rejects.tle_last_reject = reason;
    tle_metadata_pending = false;
}

static bool update_orbit(const float* posn_update, int jd_scalar, float jd_frac, float* r_eci, float* v_eci) {
    double current_jd = (double)jd_scalar + (double)jd_frac;
    if (!isfinite(current_jd)) {
        orbit_err = BODY_ORBIT_BAD_TIME;
        return false;
    }
    float new_tle[6];
    bool adopt = false;
    // Without its epoch the elements can't be placed in time, so a TLE with no metadata is rejected.
    // (No metadata does NOT use up anything: there is none pending.)
    if (posn_update != NULLPTR) {
        if (!tle_metadata_pending) {
            rejects.tle_rejected++;
            rejects.tle_last_reject = BODY_ORBIT_NO_METADATA;
        } else if (!all_finite(posn_update, 6) || !valid_tle_mean(posn_update)) {
            reject_new_tle(BODY_ORBIT_BAD_ELEMENTS);
        } else {
            memcpy(new_tle, posn_update, sizeof(new_tle));
            adopt = true;
        }
    }
    if (!adopt && !tle_initialized) { // No orbit to propagate yet
        orbit_err = BODY_ORBIT_NO_TLE;
        return false;
    }
    const float* active_tle = adopt ? new_tle : tle_mean;
    float active_bstar = adopt ? pending_bstar : tle_bstar;
    double active_epoch_jd = adopt ? pending_epoch_jd : tle_epoch_jd;

    // Propagate the orbit with SGP4 to the current time
    double seconds_since_epoch = (current_jd - active_epoch_jd) * 86400.0;
    if (!isfinite(seconds_since_epoch)) {
        if (adopt) reject_new_tle(BODY_ORBIT_BAD_TIME);
        orbit_err = BODY_ORBIT_BAD_TIME;
        return false;
    }
    float r_km[3], v_km_s[3], oe_osc[7];
    int sgp4_err = sgp4_step(active_tle, active_bstar, active_epoch_jd, (float)seconds_since_epoch, r_km, v_km_s,
                             oe_osc);
    if (sgp4_err != 0) {
        if (adopt) reject_new_tle(BODY_ORBIT_SGP4_BASE + (unsigned int)sgp4_err);
        orbit_err = BODY_ORBIT_SGP4_BASE + (unsigned int)sgp4_err;
        return false;
    }
    // WMM and pointing use metres; SGP4 returns kilometres and km/s.
    for (int i = 0; i < 3; i++) {
        r_eci[i] = (float)(r_km[i] * 1000.0);
        v_eci[i] = (float)(v_km_s[i] * 1000.0);
    }
    // Reject a degenerate propagated state
    if (!all_finite(r_eci, 3) || !all_finite(v_eci, 3) || !(l2_norm(r_eci, 3) > 0.0f)
        || !(l2_norm(v_eci, 3) > 0.0f)) {
        if (adopt) reject_new_tle(BODY_ORBIT_DEGENERATE);
        orbit_err = BODY_ORBIT_DEGENERATE;
        return false;
    }
    // Propagation succeeded, so commit the new TLE as the saved one
    if (adopt) {
        rejects.tle_accepted++;
        memcpy(tle_mean, new_tle, sizeof(tle_mean));
        tle_bstar = active_bstar;
        tle_epoch_jd = active_epoch_jd;
        tle_initialized = true;
        tle_metadata_pending = false;
    }
    memcpy(last_r_eci, r_eci, sizeof(last_r_eci));
    position_valid = true;
    orbit_err = BODY_ORBIT_OK;
    return true;
}

// The ADCS step proper. Magnetometer and gyro are already in the BODY frame here.
static void body_core(const float* last_magnetometer_measurements, // 1x3, body
                      const float* magnetometer_measurements,      // 1x3, body
                      const float* gyro_measurements,              // 1x3, body, rad/s
                      const float* photodiode_measurements,        // 1xNUM_DIODES, raw readings
                      const float* posn_update,                    // TLE mean elements, or NULLPTR
                      int unix_time,                               // unix time
                      int jd_scalar,                               // scalar JD
                      float jd_frac,                               // fractional JD
                      float dt,                                    // time since last update
                      float* output_currents) {
    // memzero output currents to start. ensures that any "bail" automatically sets output currents to 0
    if (output_currents != NULLPTR) zero_currents(output_currents);
    else return;

    // Bail on null pointers or bad dt
    if (gyro_measurements == NULLPTR || photodiode_measurements == NULLPTR) {
        rejects.null_ptr++;
        return;
    }

    if (!(dt > 0.0f) || !isfinite(dt)) {
        rejects.bad_dt++;
        return;
    }

    // While pointing, this tick's dt is owed to the filter even if we return before propagating;
    // while detumbling it counts toward the time between clean B-dot samples.
    if (pointing) unpropagated_s += dt;
    else bdot_elapsed_s += dt;

    // ---- Sensor sanity: non-finite, stuck (bit-identical repeats), out of range, one-tick spikes ----
    // Gyro: a bad sample skips the tick. A short gap is bridged by the next good tick (which
    // propagates over unpropagated_s); a long one means the gyro is gone, so stop pointing.

    bool gyro_valid = gyro_sample_ok(gyro_measurements);
    if (!gyro_valid) {
        ++rejects.gyro_fail_streak;
        if (pointing && rejects.gyro_fail_streak >= BODY_GYRO_FAULT_RESET_TICKS) {
            reset_filter_to_detumble();
        }
        //Only skip tick if pointing, no need to skip if BDOT
        if(pointing) {
            return;
        }
    }
    else {
        rejects.gyro_fail_streak = 0;
    }
    // Photodiodes: implausible or stuck readings become NaN (handled like dropouts downstream)
    float photodiodes[NUM_DIODES];
    sanitize_photodiodes(photodiode_measurements, photodiodes);

    // Magnetometer: finite and not stuck here; |B| vs the WMM is checked once the orbit is known.
    // While pointing only the current sample matters, and a bad sample is tolerated: the filter
    // coasts on the gyro and control uses the WMM field. Only a persistent fault resets the filter.
    float mag_body_unit[3] = {0.0f, 0.0f, 0.0f};
    bool mag_valid = mag_ok(magnetometer_measurements) && !mag_stuck(magnetometer_measurements);
    if (mag_valid) {
        memcpy(mag_body_unit, magnetometer_measurements, sizeof(mag_body_unit));
        mag_valid = normalize_vec(mag_body_unit, 3);
    }

    // Orbit and field model. B-dot doesn't need them, so detumbling runs before the first TLE;
    // pointing (and the hand-off to it) does.
    float r_eci[3], v_eci[3], expected_mag[3] = {0.0f, 0.0f, 0.0f}, mag_ref_unit[3] = {0.0f, 0.0f, 0.0f};
    bool have_orbit = update_orbit(posn_update, jd_scalar, jd_frac, r_eci, v_eci);
    if (have_orbit) {
        rejects.orbit_fail_streak = 0;
        rejects.orbit_fail_reason = BODY_ORBIT_OK;
    } else {
        if (rejects.orbit_fail_streak < UINT_MAX) rejects.orbit_fail_streak++;
        rejects.orbit_fail_reason = orbit_err;
    }
    float mag_corrected[3] = {0.0f, 0.0f, 0.0f};
    if (magcal_ticks < UINT_MAX) magcal_ticks++;
    if (have_orbit) {
        // Expected magnetometer reading (T, ECI)
        wmm_eci_embedded_v2(r_eci, jd_scalar, jd_frac, expected_mag);
        memcpy(mag_ref_unit, expected_mag, sizeof(mag_ref_unit));
        // Degenerate model field (exactly zero): can't use the magnetometer at all
        if (!normalize_vec(mag_ref_unit, 3)) {
            if (pointing) reset_filter_to_detumble();
            return;
        }
        // Hard-iron correction, once the on-orbit calibration has converged
        if (mag_valid) {
            for (int i = 0; i < 3; i++) mag_corrected[i] = magnetometer_measurements[i];
            if (magcal_ready()) {
                for (int i = 0; i < 3; i++) mag_corrected[i] -= magcal_b[i] * 1e-6f;
                memcpy(mag_body_unit, mag_corrected, sizeof(mag_body_unit));
                mag_valid = normalize_vec(mag_body_unit, 3);
            }
        }
        // |B| is attitude-independent: a railed or wildly wrong magnetometer disagrees with the WMM
        if (mag_valid) {
            float b_model = l2_norm(expected_mag, 3);
            float b_meas = l2_norm(mag_corrected, 3);
            if (fabsf(b_meas - b_model) > BODY_MAG_NORM_TOL_FRAC * b_model) {
                mag_valid = false;
                rejects.mag_norm++;
            }
        }
    }

    if (mag_valid) {
        mag_fault_ticks = 0;
    } else if (mag_fault_ticks < UINT_MAX) {
        mag_fault_ticks++;
    }
    if (pointing && mag_fault_ticks >= BODY_MAG_FAULT_RESET_TICKS) {
        reset_filter_to_detumble();
        return;
    }

    //If pointing and persistently loose orbit, then Bdot
    if (pointing && rejects.orbit_fail_streak >= ORBIT_MAX_FAIL_TICK) {
        reset_filter_to_detumble();
        return;
    }

    // Pointing needs a position this tick (a failed propagation skips the tick, as before)
    if (pointing && !have_orbit) return;

   

    // During detumbling...
    if (!pointing) {
        // B-dot cycle (BODY_BDOT_PERIOD_TICKS): the coils are off for BODY_BDOT_QUIET_TICKS before
        // the last tick, so only that tick's magnetometer sample is free of coil field. dB/dt comes
        // from consecutive clean samples, so last_magnetometer_measurements is no longer used.
        (void)last_magnetometer_measurements;
        bdot_tick = bdot_tick >= BODY_BDOT_PERIOD_TICKS ? 1u : bdot_tick + 1u;
        bool sample_tick = bdot_tick == BODY_BDOT_PERIOD_TICKS;
        if (sample_tick) {
            // A reference older than two cycles (ticks skipped by early returns) is too stale
            bool prev_fresh = bdot_prev_valid
                              && bdot_elapsed_s <= 2.0f * (float)BODY_BDOT_PERIOD_TICKS * dt;
            if (mag_valid) {
                if (prev_fresh) {
                    float dbdt[3];
                    for (int i = 0; i < 3; i++) {
                        dbdt[i] = (magnetometer_measurements[i] - bdot_prev_mag[i]) / bdot_elapsed_s;
                    }
                    update_dbdt_filt(dbdt);
                    Bdot(magnetometer_measurements, bdot_prev_mag, BDOT_GAIN, bdot_elapsed_s, bdot_moment);
                } else {
                    memset(bdot_moment, 0, sizeof(bdot_moment));
                }
                memcpy(bdot_prev_mag, magnetometer_measurements, sizeof(bdot_prev_mag));
                bdot_prev_valid = true;
            } else {
                // No clean sample this cycle: drop the reference and hold no dipole
                bdot_prev_valid = false;
                memset(bdot_moment, 0, sizeof(bdot_moment));
            }
            bdot_elapsed_s = 0.0f;
        }

        // Clean sample: feed the hard-iron calibration (at most every BODY_MAGCAL_MIN_TICKS)
        if (sample_tick && mag_valid && have_orbit && magcal_ticks >= BODY_MAGCAL_MIN_TICKS) {
            magcal_update(magnetometer_measurements, l2_norm(expected_mag, 3));
        }

        // Case 1: Measurements tell us to stop detumbling. Hand off only on a sample tick, so QUEST
        // gets a clean magnetometer vector.
        if (have_orbit && ready_to_point(photodiodes, gyro_measurements, gyro_valid) && sample_tick && mag_valid) {
            // We only get here right after deployment or after a filter reset, so start the
            // bias estimate from zero

            // Measured sun direction (body frame) and modeled sun direction (ECI)
            float photodiode_sun_vec[3];
            get_vec_from_photodiode_readings(photodiodes, photodiode_sun_vec);
            float expected_sun_vec[3];
            sun_vec(unix_time, expected_sun_vec);

            // Pair up body-frame and reference-frame vectors (mag, sun) for QUEST
            float body_vecs[6] = {mag_body_unit[0],      mag_body_unit[1],      mag_body_unit[2],
                                  photodiode_sun_vec[0], photodiode_sun_vec[1], photodiode_sun_vec[2]};
            float ref_vecs[6] = {mag_ref_unit[0],    mag_ref_unit[1],     mag_ref_unit[2],
                                 expected_sun_vec[0], expected_sun_vec[1], expected_sun_vec[2]};

            // Uses QUEST to get initial attitude estimate
            float estimated_q[4];
            if (quest(body_vecs, ref_vecs, 2, estimated_q)) {
                // Seed the filter with the QUEST attitude and switch to pointing
                memcpy(estimated_quat, estimated_q, sizeof(float) * 4);
                memset(estimated_gyro_bias, 0, sizeof(estimated_gyro_bias));
                memset(error_quat_state, 0, sizeof(error_quat_state));
                memcpy(error_quat_cov, INIT_P, sizeof(error_quat_cov));
                attitude_initialized = true;
                pointing = true;
                pointing_tick = 0;
                command_zero_moment(output_currents);
                return; // Start the first filter cycle on the next tick.
            }
            // Parallel sun and field vectors cannot initialize QUEST; keep detumbling.
        }
        // Otherwise rely on B-dot to keep detumbling or to hold our attitude steady.
        // Quiet period: coils off on the ticks just before the next sample tick.
        if (!sample_tick && bdot_tick >= BODY_BDOT_PERIOD_TICKS - BODY_BDOT_QUIET_TICKS) {
            command_zero_moment(output_currents);
            return;
        }
        moment2current3axis(bdot_moment, Imax, output_currents);
        apply_coil_config(output_currents);
        return;
    }

    // ---- During pointing ----
    // Spinning too fast to point: go back to detumbling
    if (l2_norm(gyro_measurements, 3) > MAX_POINTING_RATE_RAD_S) {
        reset_filter_to_detumble();
        return;
    }
    // Advance the tick counter (wraps at the update period); only the last tick takes a measurement
    pointing_tick = pointing_tick == BODY_MEASUREMENT_UPDATE_PERIOD_TICKS
                    ? 1u : pointing_tick + 1u;
    bool do_update = pointing_tick == BODY_MEASUREMENT_UPDATE_PERIOD_TICKS;
    // The update tick's sample is coil-quiet: feed the hard-iron calibration
    if (do_update && mag_valid && magcal_ticks >= BODY_MAGCAL_MIN_TICKS) {
        magcal_update(magnetometer_measurements, l2_norm(expected_mag, 3));
    }
    // On ticks where filter updates with measurement, check whether the sun is visible
    float photodiode_sun_vec[3];
    bool in_sun = do_update
                  && get_vec_from_photodiode_readings(photodiodes, photodiode_sun_vec);

    // Choose which vectors and noise matrix the UKF gets during an update tick
    float body_vecs[6] = {mag_body_unit[0], mag_body_unit[1], mag_body_unit[2], 0, 0, 0};
    float ref_vecs[6] = {mag_ref_unit[0], mag_ref_unit[1], mag_ref_unit[2], 0, 0, 0};
    int num_vecs;
    const float* R;
    if (!do_update || !mag_valid) {
        num_vecs = 0; // gyro propagation only (between measurement ticks, or a bad magnetometer)
        R = NULLPTR;
    } else if (in_sun) {
        // Sunlit: use both magnetometer and sun vectors
        float expected_sun_vec[3];
        sun_vec(unix_time, expected_sun_vec);
        memcpy(&body_vecs[3], photodiode_sun_vec, sizeof(float) * 3);
        memcpy(&ref_vecs[3], expected_sun_vec, sizeof(float) * 3);
        num_vecs = 2;
        R = R_IN_SUN;
    } else {
        // Eclipse: magnetometer only; the bias estimate keeps coasting (mode 1 in simulink.m).
        // This replaces the finite-difference "dB/dt" second vector, which divided/multiplied by
        // dt the wrong way round and used the in-sun noise matrix.
        num_vecs = 1;
        R = R_IN_SHADOW;
    }

    // Propagate over everything since the last UKF step (normally dt; longer after skipped ticks)
    float dt_prop = unpropagated_s;
    unpropagated_s = 0.0f;

    // Scale process noise by the propagation time
    float Q_step[6 * 6];
    for (int i = 0; i < 36; i++) {
        Q_step[i] = Q_RATE[i] * dt_prop;
    }

    // Predict every tick, with only every K-th tick being a predict + correction tick.
    float new_error_state[6];
    float new_estimated_quat[4];
    float new_P[6 * 6];
    // Innovation gate (drops a spiking / inconsistent vector from this update). After
    // BODY_GATE_MAX_STREAK gated updates in a row the filter, not the sensor, is suspect, so one
    // update runs ungated (it then converges, or the innovation-strike check below resets it).
    float gate = gate_streak >= BODY_GATE_MAX_STREAK ? INFINITY : BODY_INNOVATION_GATE;
    unsigned int gated_mask = 0u;
    ukf_status_t status = iterate_gated(error_quat_state, estimated_quat, error_quat_cov, body_vecs, ref_vecs,
                                        num_vecs, gyro_measurements, Q_step, R, dt_prop, gate, new_error_state,
                                        new_estimated_quat, new_P, &gated_mask);
    if (num_vecs > 0) {
        if (gated_mask != 0u) {
            gate_streak++;
            for (int v = 0; v < num_vecs; v++) {
                if (gated_mask & (1u << v)) rejects.innovation_gated++;
            }
        } else {
            gate_streak = 0u;
        }
    }

    // Filter errored or produced a bad estimate
    if (status != UKF_OK || filter_failure(new_estimated_quat, gyro_measurements, new_P, dt_prop)) {
        // Reset everything and return to detumbling just in case
        reset_filter_to_detumble();
        return;
    }

    // Persistent disagreement with the sensors means we've diverged (e.g. after a gyro glitch)
    if (num_vecs > 0) {
        if (worst_innovation(new_estimated_quat, body_vecs, ref_vecs, num_vecs) > MAX_INNOVATION_RAD) {
            if (++innovation_strikes >= MAX_INNOVATION_STRIKES) {
                reset_filter_to_detumble();
                return;
            }
        } else {
            innovation_strikes = 0;
        }
    }

    // Prepare for next run (iterate already zeroed the attitude part)
    memcpy(error_quat_state, new_error_state, sizeof(float) * 6);
    memcpy(error_quat_cov, new_P, sizeof(float) * 6 * 6);
    memcpy(estimated_quat, new_estimated_quat, sizeof(float) * 4);
    memcpy(estimated_gyro_bias, &new_error_state[3], sizeof(float) * 3);

    // Model field in body axes (for the reduced-coil re-solve)
    float q_eci_to_body[4], b_body[3];
    quat_conj(estimated_quat, q_eci_to_body);
    quat_apply(q_eci_to_body, expected_mag, b_body);

    // Coils settle before the update; convert a zero moment through the current model.
    if (pointing_tick >= BODY_MAGNETORQUER_QUIET_START_TICK) {
        command_zero_moment(output_currents);
        return;
    }

    // ---- Control: pointing_error.m -> PD_loop.m -> torque2moment3axis.m -> moment2current3axis.m
    // Bias-corrected angular rate
    float omega[3];
    for (int i = 0; i < 3; i++) {
        omega[i] = gyro_measurements[i] - estimated_gyro_bias[i];
    }

    // PD torque command
    float m_out[3];
    new_pd(r_eci, v_eci, estimated_quat, expected_mag, omega, m_out);
    // A disabled coil: put its share on the coils that still work
    project_dipole_to_enabled(m_out, b_body);

    // Moment -> coil currents (desired), then the coil configuration (polarity / enable)
    moment2current3axis(m_out, Imax, output_currents);
    apply_coil_config(output_currents);
}

// Public entry point: takes raw sensor-axis readings, maps them to the body frame, runs the step.
// NULL inputs stay NULL so body_core's bail-outs (and the null_ptr telemetry) behave as before. The
// magnetorquers are wired as body X/Y/Z, so the output currents need no mapping.
void body(const float* last_magnetometer_measurements, // 1x3, [U8, U9, U10]
          const float* magnetometer_measurements,      // 1x3, [U8, U9, U10]
          const float* gyro_measurements,              // 1x3, gyro axes, rad/s
          const float* photodiode_measurements,        // 1xNUM_DIODES, raw readings
          const float* posn_update,                    // TLE mean elements, or NULLPTR
          int unix_time,                               // unix time
          int jd_scalar,                               // scalar JD
          float jd_frac,                               // fractional JD
          float dt,                                    // time since last update
          float* output_currents) {
    float last_mag_body[3], mag_body[3], gyro_body[3];
    const float* last_mag_in = NULLPTR;
    const float* mag_in = NULLPTR;
    const float* gyro_in = NULLPTR;
    if (last_magnetometer_measurements != NULLPTR) {
        mag_sensor_to_body(last_magnetometer_measurements, last_mag_body);
        last_mag_in = last_mag_body;
    }
    if (magnetometer_measurements != NULLPTR) {
        mag_sensor_to_body(magnetometer_measurements, mag_body);
        mag_in = mag_body;
    }
    if (gyro_measurements != NULLPTR) {
        gyro_sensor_to_body(gyro_measurements, gyro_body);
        gyro_in = gyro_body;
    }
    body_core(last_mag_in, mag_in, gyro_in, photodiode_measurements, posn_update, unix_time, jd_scalar,
              jd_frac, dt, output_currents);
}
