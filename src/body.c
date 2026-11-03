#include "include/body.h"
#include "arm_math.h"
#include "include/bdot.h"
#include "include/pointing_error.h"
#include "include/ecef2eci.h"
#include "include/iterate.h"
#include "include/kepler.h"
#include "include/laextension.h"
#include "include/magnetosphere.h"
#include "include/moments2currents.h"
#include "include/orbital2eci.h"
#include "include/pd.h"
#include "include/photodiode_determination.h"
#include "include/quat.h"
#include "include/quest.h"
#include "include/sunvec.h"
#include "include/torque2moments.h"
#include "math.h"
#include "string.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define DEG2RAD ((float)M_PI / 180.0f)

// Providence, lla2ecef(41.825226, -71.418884, 0) (see eceftoeci.m)
float PVD_ECEF[3] = {6.3761e6f, -0.1387e6f, 0.0807e6f};
const float DETUMBLING_GYRO_THRESHOLD = 0.17f; // rad/s
const float DETUMBLING_MAG_THRESHOLD = 0.01f;  // |dB/dt|, T/s. TODO: tune (anything in T/s passes)
const float BDOT_GAIN = 1.0f;                  // TODO: tune (placeholder like the MATLAB k)
float Imax[3] = {1, 1, 1};                     // TEMP

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

// ---- Module state (previously `static` in body.h, which gave every includer its own copy) ----
static float kepler_posn[6];
static float estimated_quat[4];
static float estimated_gyro_bias[3];
static float error_quat_state[6];
static float error_quat_cov[6 * 6];
static float q_want_prev[4] = {1, 0, 0, 0};

static bool pointing = false;
static bool posn_initialized = false;
static bool attitude_initialized = false;
static int filter_resets = 0;
static int innovation_strikes = 0;

static const float ZERO3[3] = {0, 0, 0};

void body_reset(void) {
    memset(kepler_posn, 0, sizeof(kepler_posn));
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
    posn_initialized = false;
    attitude_initialized = false;
    filter_resets = 0;
    innovation_strikes = 0;
}

bool body_is_pointing(void) { return pointing; }

bool body_get_attitude(float* q_body_to_eci) {
    memcpy(q_body_to_eci, estimated_quat, sizeof(estimated_quat));
    return attitude_initialized;
}

void body_get_gyro_bias(float* bias) { memcpy(bias, estimated_gyro_bias, sizeof(estimated_gyro_bias)); }

int body_get_filter_resets(void) { return filter_resets; }

static void zero_currents(float* output_currents) { memcpy(output_currents, ZERO3, sizeof(ZERO3)); }

static void reset_filter_to_detumble(void) {
    memset(estimated_quat, 0, sizeof(estimated_quat));
    estimated_quat[0] = 1.0f;
    memset(estimated_gyro_bias, 0, sizeof(estimated_gyro_bias));
    memset(error_quat_state, 0, sizeof(error_quat_state));
    memset(error_quat_cov, 0, sizeof(error_quat_cov));
    attitude_initialized = false;
    pointing = false;
    innovation_strikes = 0;
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

/**
 * True when rates are low enough to stop detumbling and we can see the sun (needed for QUEST).
 *
 * Fixes vs the previous version: dB used last_mag[2] for the y and z components (and mag[1] for
 * z), and was multiplied by dt instead of divided.
 */
static bool ready_to_point(const float* photodiode_measurements, const float* gyro,
                           const float* last_magnetometer_measurements,
                           const float* magnetometer_measurements, float dt) {
    if (gyro == NULLPTR || last_magnetometer_measurements == NULLPTR
        || magnetometer_measurements == NULLPTR || photodiode_measurements == NULLPTR || !(dt > 0.0f)) {
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

// Satellite ECI position (m) and velocity (m/s) from the stored elements (angles in degrees)
static void current_eci(float* r_eci, float* v_eci) {
    float kepler_rad[6] = {kepler_posn[0],
                           kepler_posn[1],
                           kepler_posn[2] * DEG2RAD,
                           kepler_posn[3] * DEG2RAD,
                           kepler_posn[4] * DEG2RAD,
                           kepler_posn[5] * DEG2RAD};
    orbital_to_eci_posvel(kepler_rad, MU_EARTH_M3_S2, r_eci, v_eci);
}

void body(const float* last_magnetometer_measurements, // 1x3
          const float* magnetometer_measurements,      // 1x3
          const float* gyro_measurements,              // 1x3, radians
          const float* photodiode_measurements,        // 1xNUM_DIODES, raw readings
          const float* posn_update,                    // 1x6, may be NULLPTR
          int unix_time,                               // unix time
          int jd_scalar,                               // scalar JD
          float jd_frac,                               // fractional JD
          float dt,                                    // time since last update
          float* output_currents) {
    if (output_currents == NULLPTR) {
        return;
    }
    zero_currents(output_currents);

    if (magnetometer_measurements == NULLPTR || last_magnetometer_measurements == NULLPTR
        || gyro_measurements == NULLPTR || photodiode_measurements == NULLPTR || !(dt > 0.0f)
        || !isfinite(dt) || !all_finite(magnetometer_measurements, 3)
        || !all_finite(last_magnetometer_measurements, 3) || !all_finite(gyro_measurements, 3)) {
        return;
    }

    // Update position
    if (posn_update != NULLPTR && all_finite(posn_update, 6)) {
        memcpy(kepler_posn, posn_update, sizeof(float) * 6);
        posn_initialized = true;
    } else if (posn_initialized) {
        float new_posn[6];
        propogateOrbitalElements(kepler_posn[0], kepler_posn[1], kepler_posn[2], kepler_posn[3],
                                 kepler_posn[4], kepler_posn[5], dt, new_posn);
        memcpy(kepler_posn, new_posn, sizeof(float) * 6);
    }

    // Make sure position exists before doing anything else
    if (!posn_initialized) {
        return;
    }

    float r_eci[3], v_eci[3];
    current_eci(r_eci, v_eci);

    // Expected magnetometer reading (T, ECI)
    float expected_mag[3];
    wmm_eci_embedded_v2(r_eci, jd_scalar, jd_frac, expected_mag);

    // Measured / expected unit vectors. Only directions are used by QUEST and the UKF, so the
    // magnetometer's scale factor doesn't matter here.
    float mag_body_unit[3] = {magnetometer_measurements[0], magnetometer_measurements[1],
                              magnetometer_measurements[2]};
    float mag_ref_unit[3] = {expected_mag[0], expected_mag[1], expected_mag[2]};
    bool mag_ok = normalize_vec(mag_body_unit, 3) && normalize_vec(mag_ref_unit, 3);

    // During detumbling...
    if (!pointing) {
        // Case 1: Measurements tell us to stop detumbling
        if (mag_ok
            && ready_to_point(photodiode_measurements, gyro_measurements,
                              last_magnetometer_measurements, magnetometer_measurements, dt)) {
            // We only get here right after deployment or after a filter reset, so start the
            // bias estimate from zero
            float photodiode_sun_vec[3];
            get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);
            float expected_sun_vec[3];
            sun_vec(unix_time, expected_sun_vec);

            float body_vecs[6] = {mag_body_unit[0],      mag_body_unit[1],      mag_body_unit[2],
                                  photodiode_sun_vec[0], photodiode_sun_vec[1], photodiode_sun_vec[2]};
            float ref_vecs[6] = {mag_ref_unit[0],     mag_ref_unit[1],     mag_ref_unit[2],
                                 expected_sun_vec[0], expected_sun_vec[1], expected_sun_vec[2]};

            // Uses QUEST to get initial attitude estimate
            float estimated_q[4];
            if (!quest(body_vecs, ref_vecs, 2, estimated_q)) {
                return; // sun and field (nearly) parallel: try again next step
            }

            memcpy(estimated_quat, estimated_q, sizeof(float) * 4);
            memset(estimated_gyro_bias, 0, sizeof(float) * 3);
            memset(error_quat_state, 0, sizeof(float) * 6);
            memcpy(error_quat_cov, INIT_P, sizeof(float) * 6 * 6);
            attitude_initialized = true;
            pointing = true;
            return;
        }
        // Otherwise rely on B-dot to keep detumbling or to hold our attitude steady
        float moments[3];
        Bdot(magnetometer_measurements, last_magnetometer_measurements, BDOT_GAIN, dt, moments);
        moment2current3axis(moments, Imax, output_currents);
        return;
    }

    // ---- During pointing ----
    float photodiode_sun_vec[3];
    bool in_sun = get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);

    float body_vecs[6] = {mag_body_unit[0], mag_body_unit[1], mag_body_unit[2], 0, 0, 0};
    float ref_vecs[6] = {mag_ref_unit[0], mag_ref_unit[1], mag_ref_unit[2], 0, 0, 0};
    int num_vecs;
    const float* R;
    if (!mag_ok) {
        num_vecs = 0; // no usable field: gyro propagation only
        R = R_IN_SHADOW;
    } else if (in_sun) {
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

    float Q_step[6 * 6];
    for (int i = 0; i < 36; i++) {
        Q_step[i] = Q_RATE[i] * dt;
    }

    // Run the filter
    float new_error_state[6];
    float new_estimated_quat[4];
    float new_P[6 * 6];
    ukf_status_t status = iterate(error_quat_state, estimated_quat, error_quat_cov, body_vecs, ref_vecs,
                                  num_vecs, gyro_measurements, Q_step, R, dt, new_error_state,
                                  new_estimated_quat, new_P);

    if (status != UKF_OK || filter_failure(new_estimated_quat, gyro_measurements, new_P, dt)) {
        // Reset everything and return to detumbling just in case
        reset_filter_to_detumble();
        return;
    }

    // Persistent disagreement with the sensors means we've diverged (e.g. after a gyro glitch)
    if (num_vecs > 0 && worst_innovation(new_estimated_quat, body_vecs, ref_vecs, num_vecs) > MAX_INNOVATION_RAD) {
        if (++innovation_strikes >= MAX_INNOVATION_STRIKES) {
            reset_filter_to_detumble();
            return;
        }
    } else {
        innovation_strikes = 0;
    }

    // Prepare for next run (iterate already zeroed the attitude part)
    memcpy(error_quat_state, new_error_state, sizeof(float) * 6);
    memcpy(error_quat_cov, new_P, sizeof(float) * 6 * 6);
    memcpy(estimated_quat, new_estimated_quat, sizeof(float) * 4);
    memcpy(estimated_gyro_bias, &new_error_state[3], sizeof(float) * 3);

    // ---- Control: pointing_error.m -> PD_loop.m -> torque2moment3axis.m -> moment2current3axis.m
    float pvd_eci[3];
    ecef_2_eci(PVD_ECEF, pvd_eci, unix_time);

    // Body-frame error quaternion. (The old code used goal * est^-1, an ECI-frame error, and
    // fed it to a controller that works in body coordinates.)
    float q_err[4];
    if (!pointing_error(r_eci, v_eci, estimated_quat, pvd_eci, q_want_prev, q_err)) {
        return;
    }

    float omega[3];
    for (int i = 0; i < 3; i++) {
        omega[i] = gyro_measurements[i] - estimated_gyro_bias[i];
    }

    float t[3];
    pd_loop(q_err, omega, t);

    float m[3];
    torque_2_moments(magnetometer_measurements, t, m);

    moment2current3axis(m, Imax, output_currents);
}
