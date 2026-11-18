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
 */

typedef struct {
    int axis[3];   // sensor axis feeding body axis i (a permutation of 0, 1, 2)
    float sign[3]; // +1 or -1
} sensor_axis_map_t;

extern const sensor_axis_map_t MAG_SENSOR_TO_BODY;
extern const sensor_axis_map_t GYRO_SENSOR_TO_BODY;

// Flight path: raw sensor reading -> body frame. `sensor` and `body_vec` must not alias.
void mag_sensor_to_body(const float* sensor, float* body_vec);
void gyro_sensor_to_body(const float* sensor, float* body_vec);

// Ground / test path: body-frame truth -> what the sensor would read. Must not alias.
void mag_body_to_sensor(const float* body_vec, float* sensor);
void gyro_body_to_sensor(const float* body_vec, float* sensor);

#endif
