#ifndef BODY
#define BODY
#include <stdbool.h>

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

/**
 * One ADCS step.
 *
 * Units / conventions:
 *   magnetometer: body frame, TESLA (same units as the WMM model and torque_2_moments)
 *   gyro:         body frame, rad/s
 *   photodiodes:  NUM_DIODES raw readings in PHOTODIODES order
 *   posn_update:  TLE mean [n (rev/day), e, i, RAAN, argp, M (rad)] or NULLPTR to reuse the last TLE
 *                 Without metadata, a new TLE uses B* = 0 and the current Julian date as epoch.
 *   unix_time / jd_scalar + jd_frac: the same instant, UTC
 *   output_currents: magnetorquer currents (A); always written when non-NULL
 *                    Zero-moment commands also pass through the current converter.
 */
void body(const float* last_magnetometer_measurements, // 1x3
          const float* magnetometer_measurements,      // 1x3
          const float* gyro_measurements,              // 1x3, radians
          const float* photodiode_measurements,        // 1xNUM_DIODES, raw readings
          const float* posn_update,                    // 1x6, may be NULLPTR
          int unix_time,                               // unix time
          int jd_scalar,
          float jd_frac,
          float dt,                                    // time since last update
          float* output_currents);                     // 1x3

// Clears all filter / mode state (power-on state)
void body_reset(void);
bool body_is_pointing(void);
// Current body->ECI estimate; returns false if the attitude hasn't been initialized
bool body_get_attitude(float* q_body_to_eci);
void body_get_gyro_bias(float* bias);
// Number of times the filter has been declared failed and reset to detumbling
int body_get_filter_resets(void);
/* Apply B* and UTC epoch to the next non-NULL posn_update; consumed once. */
bool body_set_tle_metadata(float bstar, double epoch_jd);
#endif
