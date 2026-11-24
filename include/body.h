#ifndef BODY
#define BODY
#include "include/photodiode_determination.h"
#include <stdbool.h>
#include <stdint.h>

#ifndef NULLPTR
#define NULLPTR 0x0
#endif

extern float PVD_ECEF[3];
extern const float DETUMBLING_MAG_THRESHOLD;
extern const float DETUMBLING_GYRO_THRESHOLD;
extern const float BDOT_GAIN;
extern const float* Imax;
extern const float INIT_P[6 * 6];
extern const float Q_RATE[6 * 6];
extern const float R_IN_SUN[6 * 6];
extern const float R_IN_SHADOW[3 * 3];

/* Update on tick 100; allow the magnetorquers to settle on ticks 95-100. */
#define BODY_MEASUREMENT_UPDATE_PERIOD_TICKS 100u
#define BODY_MAGNETORQUER_QUIET_START_TICK 95u
#if BODY_MAGNETORQUER_QUIET_START_TICK < 1u \
    || BODY_MAGNETORQUER_QUIET_START_TICK > BODY_MEASUREMENT_UPDATE_PERIOD_TICKS
#error "Magnetorquer quiet start must be within the measurement cycle"
#endif
/* Consecutive bad magnetometer ticks tolerated while pointing before the filter is reset.
 * Shorter runs are treated as glitches: the filter coasts on the gyro and skips the update. */
#define BODY_MAG_FAULT_RESET_TICKS (3u * BODY_MEASUREMENT_UPDATE_PERIOD_TICKS)
/* Consecutive non-finite gyro ticks tolerated while pointing before the filter is reset.
 * Shorter gaps are bridged: the next good tick propagates over the whole gap. */
#define BODY_GYRO_FAULT_RESET_TICKS 100u
/* ---- Sensor sanity checks (body.c: gyro_sample_ok, mag_stuck, sanitize_photodiodes) ---- */
/* Bit-identical consecutive repeats that mark a sensor as stuck. With RM3100 noise ~ 1 LSB, all 3
 * axes repeating once is ~1/70, so 10 in a row never happens to a live sensor. */
#define BODY_SENSOR_STUCK_TICKS 10u
/* Largest believable one-tick gyro change (rad/s). The coils give < 2e-4 rad/s per tick; noise is
 * ~1e-4 rad/s. A bigger step is a spike unless it persists for BODY_GYRO_STEP_ACCEPT_TICKS. */
#define BODY_GYRO_MAX_STEP_RAD_S 0.02f
#define BODY_GYRO_STEP_ACCEPT_TICKS 10u
/* |B| measured vs |B| WMM: reject the sample if they differ by more than this fraction of |B| WMM
 * (field-model error and calibrated biases are < ~5%; a railed axis is 800 uT). */
#define BODY_MAG_NORM_TOL_FRAC 0.25f
/* Photodiodes: readings above this fraction of MAX_READING are railed (sun + albedo stays below
 * ~1.2x). Light can't be negative: readings below -BODY_PD_NEG_TOL are a fault (small negatives are
 * noise / ADC offset around 0). Only diodes above MIN_READING get the stuck check (dark diodes
 * legitimately repeat 0). */
#define BODY_PD_MAX_PLAUSIBLE_FRAC 1.25f
#define BODY_PD_NEG_TOL 0.05f
/* UKF innovation gate: a measured vector with nu' S^-1 nu above this is left out of the update
 * (chi-square, 3 dof, p = 0.001: ~1 false rejection per 1000 good vectors). After this many gated
 * updates in a row, one update runs ungated so a drifted filter can't lock itself out. */
#ifndef BODY_INNOVATION_GATE
#define BODY_INNOVATION_GATE 16.27f
#endif
#define BODY_GATE_MAX_STREAK 5u

/* On-orbit magnetometer hard-iron calibration (attitude-independent: |B_meas - b| = |B_WMM|),
 * recursive least squares on coil-quiet samples, at most one per BODY_MAGCAL_MIN_TICKS. The
 * estimate is subtracted from the measurement once every axis' std is below BODY_MAGCAL_READY_UT. */
#define BODY_MAGCAL_MIN_TICKS 100u
#define BODY_MAGCAL_FORGET 0.998f       /* per sample: ~500 samples (~1.4 h) memory, tracks bias steps */
#define BODY_MAGCAL_P0_UT2 25.0f        /* prior variance per axis (5 uT std) */
#define BODY_MAGCAL_R_UT4 256.0f        /* |B|^2 noise: (2 |B| sigma_model)^2, |B| ~ 40 uT, ~200 nT model error */
#ifndef BODY_MAGCAL_READY_UT
#define BODY_MAGCAL_READY_UT 0.3f
#endif
#define BODY_MAGCAL_MAX_UT 10.0f        /* larger estimates are rejected (not a hard-iron residual) */
/* Current hard-iron estimate (T, body); returns true once it is being applied */
bool body_get_mag_bias(float* bias_T);

/* ---- Magnetorquer configuration ----
 * Per-axis polarity (+1 / -1) and enable, applied to every coil command. Set from the ground
 * (e.g. a wiring error found after launch, or a dead coil). With an axis disabled, the pointing
 * dipole is re-solved on the remaining coils. NOTE: body_reset() restores +1 / all enabled, so flight
 * software must re-apply a stored configuration after every reset. */
void body_set_coil_config(const float* sign, const bool* enabled);
void body_get_coil_config(float* sign, bool* enabled);

/* ---- Sensor -> body maps ----
 * Replace the magnetometer and/or gyro axis maps from the ground (a sensor found mounted or wired the
 * wrong way after launch). Codes are int8[3], one per body axis: +-1 / +-2 / +-3 = +- sensor x / y / z
 * (see sensor_frames.h); NULL leaves that sensor's map alone. Both are validated before either changes:
 * returns false, changing nothing, unless each names every sensor axis once and the gyro map is a
 * proper rotation (det +1). A real change drops the attitude estimate back to detumbling (not counted
 * in body_get_filter_resets()), and a magnetometer change also restarts the hard-iron calibration.
 * Re-applying the current maps changes nothing. NOTE: body_reset() restores the as-built maps, so
 * flight software must re-apply stored maps after every reset, as for the coil configuration. */
bool body_set_sensor_maps(const int8_t* mag_codes, const int8_t* gyro_codes);
void body_get_sensor_maps(int8_t* mag_codes, int8_t* gyro_codes);

/* Orbit / TLE status codes (tle_last_reject, orbit_fail_reason) */
#define BODY_ORBIT_OK 0u
#define BODY_ORBIT_NO_METADATA 1u    /* uplink: elements without fresh body_set_tle_metadata */
#define BODY_ORBIT_BAD_ELEMENTS 2u   /* uplink: non-finite, n <= 0, or e outside [0, 1) */
#define BODY_ORBIT_BAD_TIME 3u       /* non-finite current JD or time since epoch (defensive) */
#define BODY_ORBIT_DEGENERATE 4u     /* SGP4 returned a non-finite / zero position or velocity */
#define BODY_ORBIT_NO_TLE 5u         /* orbit: no TLE has ever been adopted (normal before 1st uplink) */
#define BODY_ORBIT_TIME_INVALID 6u   /* orbit: the caller passed time_valid = false (no orbit, B-dot only) */
#define BODY_ORBIT_TIME_EPOCH_MISMATCH 7u /* the time is implausibly far from the TLE epoch: the clock (or the
                                           * uplinked epoch) is wrong. See BODY_CLOCK_MAX_*_EPOCH_DAYS */
#define BODY_ORBIT_SGP4_BASE 10u     /* 10 + SGP4 error: 11 e out of range, 12 n <= 0, 14 p < 0,
                                      * 16 decayed (r < 1 Earth radius), 17 period >= 225 min */

/* Clock plausibility: the time passed to body() must be within this window around the TLE epoch, or the
 * orbit (and with it the WMM field and sun vector) is treated as unusable. A clock that reset to a default
 * date, or a wrong uplinked time, lands far outside it. The TLE side is generous so that a TLE left un-refreshed
 * for weeks keeps pointing (SGP4 error grows ~ km/day, which costs well under a degree at that age). */
#define BODY_CLOCK_MAX_BEFORE_EPOCH_DAYS 2.0
#define BODY_CLOCK_MAX_AFTER_EPOCH_DAYS 120.0

/* Sensor-check telemetry. Fields are running totals of samples rejected since body_reset() (they
 * only go up; the ground differences consecutive packets), except the *_streak fields, which count
 * consecutive bad ticks right now and return to 0 on the next good sample. */
typedef struct {
    unsigned int gyro_stuck, gyro_spike, mag_stuck, mag_norm, mag_nan, innovation_gated, null_ptr, bad_dt;
    uint32_t pd_rejected[NUM_DIODES];       /* per diode: NaN in, stuck, railed, or negative */
    uint32_t pd_rejects_streak[NUM_DIODES]; /* streak: consecutive rejected ticks, per diode */
    /* TLE uplinks (posn_update != NULLPTR): adopted, or rejected. A rejected uplink is ignored, the
     * previous TLE stays in use, and its metadata is used up: a retry needs body_set_tle_metadata
     * AND the elements again. tle_last_reject: BODY_ORBIT_* code of the most recent rejection.
     * TLE age is epoch_jd from body_get_tle. */
    unsigned int tle_accepted, tle_rejected;
    unsigned int tle_last_reject;
    /* Orbit state (not uplinks): consecutive ticks with no usable orbit (0 = orbit OK now) and why
     * (BODY_ORBIT_* code; 0 while OK). No orbit = no pointing; detumbling (B-dot) continues. */
    unsigned int orbit_fail_streak;
    unsigned int orbit_fail_reason;

    /* NOT a total: consecutive ticks the gyro has been unusable right now (any reason). 0 = healthy;
     * counts in every mode; returns to 0 on the next good sample. Pointing resets to detumble at
     * BODY_GYRO_FAULT_RESET_TICKS. */
    unsigned int gyro_fail_streak;
} body_sensor_rejects_t;
void body_get_sensor_rejects(body_sensor_rejects_t* out);

/* Detumbling B-dot cycle. The coils corrupt the magnetometer while they are driven, so (like the
 * pointing cycle) only one sample per cycle is used: the coils are commanded off for the
 * BODY_BDOT_QUIET_TICKS ticks before the sample tick (the last tick of the cycle), and dB/dt is the
 * difference of consecutive clean samples (one cycle apart). The new dipole is commanded on the
 * sample tick and held until the next quiet period. Duty cycle = 1 - QUIET/PERIOD.
 * One cycle (1 s) between samples aliases only above 180 deg/s and cuts dB/dt noise 10x vs 0.1 s. */
#define BODY_BDOT_PERIOD_TICKS 10u
#define BODY_BDOT_QUIET_TICKS 2u
#if BODY_BDOT_QUIET_TICKS < 1u || BODY_BDOT_QUIET_TICKS >= BODY_BDOT_PERIOD_TICKS
#error "B-dot quiet ticks must be at least 1 and shorter than the B-dot cycle"
#endif
/* Detumble -> pointing hand-off. Both |gyro| and the low-passed |dB/dt| must stay under their
 * thresholds (with the sun visible) for this many consecutive ticks; any tick that fails resets
 * the count. The hand-off itself happens on a B-dot sample tick (clean magnetometer for QUEST). */
#define BODY_READY_TO_POINT_TICKS 50u
/* First-order low-pass on the clean-sample dB/dt used by the hand-off check (alpha per B-dot
 * cycle; 0.5 per 1 s cycle is a ~1.4 s time constant). */
#define BODY_DBDT_FILTER_ALPHA 0.5f

/**
 * One ADCS step.
 *
 * Units / conventions:
 *   magnetometer: SENSOR axes [U8, U9, U10], TESLA (same units as the WMM model)
 *   gyro:         SENSOR (gyro chip) axes, rad/s
 *                 Both are mapped to the body frame internally; see sensor_frames.h.
 *   photodiodes:  NUM_DIODES raw readings in PHOTODIODES order
 *   posn_update:  TLE mean [n (rev/day), e, i, RAAN, argp, M (rad)] or NULLPTR to reuse the last TLE
 *                 Call body_set_tle_metadata() first; a new TLE without metadata is rejected.
 *   unix_time / jd_scalar + jd_frac: the same instant, UTC
 *   time_valid:   false if the caller's clock isn't trusted (never set since boot, restored after an outage of
 *                 unknown length, ...). body() then computes no orbit: it keeps detumbling with B-dot (which needs
 *                 only the magnetometer), doesn't hand off to pointing, falls back to B-dot after
 *                 ORBIT_MAX_FAIL_TICK ticks if it was pointing, and leaves a pending TLE untouched (not
 *                 counted as rejected) until the time is valid. orbit_fail_reason = BODY_ORBIT_TIME_INVALID.
 *   output_currents: magnetorquer currents (A), body X/Y/Z = C2/C1/Air; always written when non-NULL
 *                    Zero-moment commands also pass through the current converter.
 */
void body(const float* last_magnetometer_measurements, // 1x3, [U8, U9, U10]
          const float* magnetometer_measurements,      // 1x3, [U8, U9, U10]
          const float* gyro_measurements,              // 1x3, gyro axes, rad/s
          const float* photodiode_measurements,        // 1xNUM_DIODES, raw readings
          const float* posn_update,                    // 1x6, may be NULLPTR
          int unix_time,                               // unix time
          int jd_scalar,
          float jd_frac,
          bool time_valid,                             // false: clock not trusted, no orbit (B-dot only)
          float dt,                                    // time since last update
          float* output_currents);                     // 1x3

// Clears all filter / mode state (power-on state)
void body_reset(void);
bool body_is_pointing(void);
/* Detumble-only mode (ground command): body() keeps running B-dot and never hands off to pointing; if it was
 * pointing, it drops back to B-dot on the next step. All of the usual checks and telemetry keep running.
 * body_reset() turns it off; the caller re-applies its mode after a reset. */
void body_set_detumble_only(bool on);
bool body_is_detumble_only(void);
// Current body->ECI estimate; returns false if the attitude hasn't been initialized
bool body_get_attitude(float* q_body_to_eci);
void body_get_gyro_bias(float* bias);
// Number of times the filter has been declared failed and reset to detumbling
int body_get_filter_resets(void);
/* Apply B* and UTC epoch to the next non-NULL posn_update; consumed once. */
bool body_set_tle_metadata(float bstar, double epoch_jd);
// The TLE body() last committed (mean elements in posn_update order); returns false if there isn't one yet
bool body_get_tle(float* tle_mean_out, float* bstar, double* epoch_jd);
// ECI position (m) from the last successful propagation; returns false if there hasn't been one yet
bool body_get_position(float* r_eci_m);
#endif
