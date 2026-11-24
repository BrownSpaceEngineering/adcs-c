#ifndef SENSOR_FRAMES_H
#define SENSOR_FRAMES_H
/*
 * Sensor axes -> satellite body frame.
 *
 * Each device is mounted so its own axes are permuted / flipped relative to the body frame.
 * Hardware alignment table (body axis <- device axis):
 *
 *   body | magnetometer | gyro | magnetorquer
 *   -----+--------------+------+-------------
 *    X   |   U8         |  -Y  |  C2
 *    Y   |   U9         |  -X  |  C1
 *    Z   |  -U10        |  -Z  |  Air
 *   (S-band antenna on +Z)
 *
 * Magnetometer readings arrive as [U8, U9, U10]; gyro readings in the gyro's own [x, y, z].
 * The magnetorquers are already labelled and wired as body X/Y/Z, so currents need no mapping,
 * and the PHOTODIODES normals are already in body axes.
 *
 * Each mapping is body[i] = SIGN[i] * sensor[AXIS[i]] (a signed permutation), so the inverse is
 * sensor[AXIS[i]] = SIGN[i] * body[i]. The arithmetic is exact and NaNs stay on their own axis.
 *
 * The table above is the hardware as built (the *_SENSOR_TO_BODY defaults). The flight path uses an
 * active copy of each map that the ground can replace (e.g. a sensor found mounted or wired the wrong
 * way after launch). The ground / test path always uses the as-built defaults: it models the physical
 * hardware, which an uplink can't change.
 *
 * Uplink encoding: one int8 per body axis, +-1 / +-2 / +-3 = +- sensor x / y / z. E.g. the
 * magnetometer default is {1, 2, -3}, the gyro default {-2, -1, -3}. Valid codes name each sensor
 * axis exactly once (48 maps: 24 rotations, 24 reflections). A gyro map must also be a proper rotation
 * (det +1): a single 3-axis chip can't be mirrored. The magnetometer's three separate coils can.
 */
#include <stdbool.h>
#include <stdint.h>

typedef struct {
    int axis[3];   // sensor axis feeding body axis i (a permutation of 0, 1, 2)
    float sign[3]; // +1 or -1
} sensor_axis_map_t;

// As-built hardware maps (the defaults)
extern const sensor_axis_map_t MAG_SENSOR_TO_BODY;
extern const sensor_axis_map_t GYRO_SENSOR_TO_BODY;

// Flight path: raw sensor reading -> body frame, with the ACTIVE maps. `sensor` and `body_vec` must not alias.
void mag_sensor_to_body(const float* sensor, float* body_vec);
void gyro_sensor_to_body(const float* sensor, float* body_vec);

// Ground / test path: body-frame truth -> what the sensor would read, with the AS-BUILT maps. Must not alias.
void mag_body_to_sensor(const float* body_vec, float* sensor);
void gyro_body_to_sensor(const float* body_vec, float* sensor);

// Uplink codes <-> map. sensor_map_from_codes returns false (and leaves *out alone) unless the codes
// name each sensor axis exactly once.
bool sensor_map_from_codes(const int8_t* codes, sensor_axis_map_t* out);
void sensor_map_to_codes(const sensor_axis_map_t* m, int8_t* codes);
// Determinant of a valid map: +1 rotation, -1 reflection
float sensor_map_det(const sensor_axis_map_t* m);

// Active maps. The setters validate (and the gyro one requires det +1); on false nothing changes.
// Use body_set_sensor_maps() in flight: it also resets the state that depends on the frames.
bool mag_set_sensor_map(const sensor_axis_map_t* m);
bool gyro_set_sensor_map(const sensor_axis_map_t* m);
void mag_get_sensor_map(sensor_axis_map_t* m);
void gyro_get_sensor_map(sensor_axis_map_t* m);
// Restore both active maps to the as-built defaults
void sensor_maps_reset(void);

#endif
