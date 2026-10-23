#include "include/body.h"
#include "arm_math.h"
#include "include/bdot.h"
#include "include/down_quat.h"
#include "include/ecef2eci.h"
#include "include/iterate.h"
#include "include/laextension.h"
#include "include/magnetosphere.h"
#include "include/moments2currents.h"
#include "include/pd.h"
#include "include/photodiode_determination.h"
#include "include/quat.h"
#include "include/quest.h"
#include "include/sgp4_step.h"
#include "include/sunvec.h"
#include "include/torque2moments.h"
#include "math.h"
#include "string.h"

// Providence, lla2ecef(41.825226, -71.418884, 0) (see eceftoeci.m)
float PVD_ECEF[3] = {6.3761e6f, -0.1387e6f, 0.0807e6f};
const float DETUMBLING_GYRO_THRESHOLD = 0.17f; // rad/s
const float DETUMBLING_MAG_THRESHOLD = 0.01f;  // |dB/dt|, T/s. TODO: tune (anything in T/s passes)
const float BDOT_GAIN = 1.0f;                  // TODO: tune (placeholder like the MATLAB k)
const float* Imax = NULLPTR;                     // NULLPTR uses hardware defaults (see moments2currents.c)

// Filter tuning. Measurement vectors are unit vectors, so R is in (unit vector)^2.
// Initial bias variance and process noise follow simulink.m (Q there is per 0.1 s step).
// clang-format off
const float INIT_P[6 * 6] = {
    0.01f, 0,     0,     0,       0,       0,
    0,     0.01f, 0,     0,       0,       0,
    0,     0,     0.01f, 0,       0,       0,
    0,     0,     0,     2.5e-5f, 0,       0,
    0,     0,     0,     0,       2.5e-5f, 0,
    0,     0,     0,     0,       0,       2.5e-5f};
// Process noise per second (multiplied by dt each step)
const float Q_RATE[6 * 6] = {
    4e-5f, 0,     0,     0,     0,     0,
    0,     4e-5f, 0,     0,     0,     0,
    0,     0,     4e-5f, 0,     0,     0,
    0,     0,     0,     1e-9f, 0,     0,
    0,     0,     0,     0,     1e-9f, 0,
    0,     0,     0,     0,     0,     1e-9f};
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
#define MAX_POINTING_RATE_RAD_S 1.0f

// ---- Module state (previously `static` in body.h, which gave every includer its own copy) ----
static float tle_mean[6];
static float tle_bstar;
static double tle_epoch_jd;
static float pending_bstar;
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

static const float ZERO3[3] = {0, 0, 0};

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
    attitude_initialized = false;
    filter_resets = 0;
    innovation_strikes = 0;
    pointing_tick = 0;
}

bool body_set_tle_metadata(float bstar, double epoch_jd) {
    if (!isfinite(bstar) || !isfinite(epoch_jd)) {
        return false;
    }
    pending_bstar = bstar;
    pending_epoch_jd = epoch_jd;
    tle_metadata_pending = true;
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
    pointing_tick = 0;
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
    if (measurement == NULLPTR || !all_finite(measurement, 3)) return false;
    float norm = l2_norm(measurement, 3);
    return norm > 1e-30f && isfinite(norm);
}

/**
 * True when rates are low enough to stop detumbling and we can see the sun (needed for QUEST).
 *
 * Fixes vs the previous version: dB used last_mag[2] for the y and z components (and mag[1] for
 * z), and was multiplied by dt instead of divided.
 */
static bool ready_to_point(const float* photodiode_measurements, const float* gyro,
                           const float* last_magnetometer_measurements,
                           const float* magnetometer_measurements, float dt) {
    if (gyro == NULLPTR || photodiode_measurements == NULLPTR || !(dt > 0.0f)
        || !mag_ok(last_magnetometer_measurements) || !mag_ok(magnetometer_measurements)) {
        return false;
    }
    float dM[3];
    for (int i = 0; i < 3; i++) {
        dM[i] = (magnetometer_measurements[i] - last_magnetometer_measurements[i]) / dt;
    }
    float photodiode_sun_vec[3];
    bool in_sun = get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);

    return (l2_norm(dM, 3) < DETUMBLING_MAG_THRESHOLD || l2_norm(gyro, 3) < DETUMBLING_GYRO_THRESHOLD)
           && in_sun;
}

static bool valid_tle_mean(const float* elements) {
    if (!isfinite(elements[0]) || elements[0] <= 0.0 || !isfinite(elements[1])
        || elements[1] < 0.0 || elements[1] >= 1.0) return false;
    for (int i = 2; i < 6; i++) {
        if (!isfinite(elements[i])) return false;
    }
    return true;
}

void body(const float* last_magnetometer_measurements, // 1x3
          const float* magnetometer_measurements,      // 1x3
          const float* gyro_measurements,              // 1x3, radians
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

    // Bad magnetometer data: drop back to detumbling and bail
    if (!mag_ok(last_magnetometer_measurements) || !mag_ok(magnetometer_measurements)) {
        if (pointing) reset_filter_to_detumble();
        return;
    }
    // Bail on null pointers, bad dt, or non-finite gyro
    if (output_currents == NULLPTR || gyro_measurements == NULLPTR
        || photodiode_measurements == NULLPTR || !(dt > 0.0f)
        || !isfinite(dt) || !all_finite(gyro_measurements, 3)) {
        return;
    }

    // Unit-vector the measured field; a zero-length reading is unusable
    float mag_body_unit[3] = {magnetometer_measurements[0], magnetometer_measurements[1],
                              magnetometer_measurements[2]};
    if (!normalize_vec(mag_body_unit, 3)) {
        if (pointing) reset_filter_to_detumble();
        return;
    }

    // A new TLE replaces the saved epoch elements; otherwise propagate the saved TLE.
    double current_jd = (double)jd_scalar + (double)jd_frac;
    if (!isfinite(current_jd)) return;
    float new_tle[6];
    const float* active_tle = tle_mean;
    float active_bstar = tle_bstar;
    double active_epoch_jd = tle_epoch_jd;
    if (posn_update != NULLPTR) {
        // Validate the incoming TLE before adopting it
        if (!all_finite(posn_update, 6)) return;
        memcpy(new_tle, posn_update, sizeof(new_tle));
        if (!valid_tle_mean(new_tle)) return;
        active_tle = new_tle;
        active_bstar = tle_metadata_pending ? pending_bstar : 0.0;
        active_epoch_jd = tle_metadata_pending ? pending_epoch_jd : current_jd;
    } else if (!tle_initialized) {
        return; // No orbit to propagate yet.
    }

    // Propagate the orbit with SGP4 to the current time
    double seconds_since_epoch = (current_jd - active_epoch_jd) * 86400.0;
    if (!isfinite(seconds_since_epoch)) return;
    float r_km[3], v_km_s[3], oe_osc[7];
    if (sgp4_step(active_tle, active_bstar, active_epoch_jd, (float)seconds_since_epoch,
                  r_km, v_km_s, oe_osc) != 0) return;

    // WMM and pointing use metres; SGP4 returns kilometres and km/s.
    float r_eci[3], v_eci[3];
    for (int i = 0; i < 3; i++) {
        r_eci[i] = (float)(r_km[i] * 1000.0);
        v_eci[i] = (float)(v_km_s[i] * 1000.0);
    }
    // Reject a degenerate propagated state
    if (!all_finite(r_eci, 3) || !all_finite(v_eci, 3)
        || !(l2_norm(r_eci, 3) > 0.0f) || !(l2_norm(v_eci, 3) > 0.0f)) return;
    // Propagation succeeded, so commit the new TLE as the saved one
    if (posn_update != NULLPTR) {
        memcpy(tle_mean, new_tle, sizeof(tle_mean));
        tle_bstar = active_bstar;
        tle_epoch_jd = active_epoch_jd;
        tle_initialized = true;
        tle_metadata_pending = false;
    }

    // Expected magnetometer reading (T, ECI)
    float expected_mag[3];
    wmm_eci_embedded_v2(r_eci, jd_scalar, jd_frac, expected_mag);

    // Measured / expected normalized unit vectors.
    float mag_ref_unit[3] = {expected_mag[0], expected_mag[1], expected_mag[2]};
    // If degenerate model field, bail
    if (!normalize_vec(mag_ref_unit, 3)) {
        if (pointing) reset_filter_to_detumble();
        return;
    }

    // During detumbling...
    if (!pointing) {
        // Case 1: Measurements tell us to stop detumbling
        if (ready_to_point(photodiode_measurements, gyro_measurements,
                           last_magnetometer_measurements, magnetometer_measurements, dt)) {
            // We only get here right after deployment or after a filter reset, so start the
            // bias estimate from zero

            // Measured sun direction (body frame) and modeled sun direction (ECI)
            float photodiode_sun_vec[3];
            get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);
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
        // Otherwise rely on B-dot to keep detumbling or to hold our attitude steady
        float moments[3];
        Bdot(magnetometer_measurements, last_magnetometer_measurements, BDOT_GAIN, dt, moments);
        moment2current3axis(moments, Imax, output_currents);
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
    // On ticks where filter updates with measurement, check whether the sun is visible
    float photodiode_sun_vec[3];
    bool in_sun = do_update
                  && get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);

    // Choose which vectors and noise matrix the UKF gets during an update tick
    float body_vecs[6] = {mag_body_unit[0], mag_body_unit[1], mag_body_unit[2], 0, 0, 0};
    float ref_vecs[6] = {mag_ref_unit[0], mag_ref_unit[1], mag_ref_unit[2], 0, 0, 0};
    int num_vecs;
    const float* R;
    if (!do_update) {
        num_vecs = 0; // scheduled gyro propagation between measurement ticks
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

    // Scale process noise by the timestep
    float Q_step[6 * 6];
    for (int i = 0; i < 36; i++) {
        Q_step[i] = Q_RATE[i] * dt;
    }

    // Predict every tick, with only every K-th tick being a predict + correction tick.
    float new_error_state[6];
    float new_estimated_quat[4];
    float new_P[6 * 6];
    ukf_status_t status = iterate(error_quat_state, estimated_quat, error_quat_cov, body_vecs, ref_vecs,
                                  num_vecs, gyro_measurements, Q_step, R, dt, new_error_state,
                                  new_estimated_quat, new_P);

    // Filter errored or produced a bad estimate
    if (status != UKF_OK || filter_failure(new_estimated_quat, gyro_measurements, new_P, dt)) {
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

    // Coils settle before the update; convert a zero moment through the current model.
    if (pointing_tick >= BODY_MAGNETORQUER_QUIET_START_TICK) {
        command_zero_moment(output_currents);
        return;
    }

    // ---- Control: pointing_error.m -> PD_loop.m -> torque2moment3axis.m -> moment2current3axis.m
    // Target ground location in ECI
    float pvd_eci[3];
    ecef_2_eci(PVD_ECEF, pvd_eci, unix_time);

    // Body-frame error quaternion. (The old code used goal * est^-1, an ECI-frame error, and
    // fed it to a controller that works in body coordinates.)
    float q_err[4];
    // No valid pointing solution: coast with the coils off
    if (!pointing_error(r_eci, v_eci, estimated_quat, pvd_eci, q_want_prev, q_err, NULLPTR)) {
        command_zero_moment(output_currents);
        return;
    }

    // Bias-corrected angular rate
    float omega[3];
    for (int i = 0; i < 3; i++) {
        omega[i] = gyro_measurements[i] - estimated_gyro_bias[i];
    }

    // PD torque command
    float t[3];
    pd_loop(q_err, omega, t);

    // Torque -> magnetic moment (perpendicular to the measured field)
    float m[3];
    torque_2_moments(magnetometer_measurements, t, m);

    // Moment -> coil currents
    moment2current3axis(m, Imax, output_currents);
}
