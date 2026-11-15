#include "include/sensor_frames.h"

// See the alignment table in sensor_frames.h.
// clang-format off
const sensor_axis_map_t MAG_SENSOR_TO_BODY = {
    .axis = {0,    1,    2},     // body X = U8, body Y = U9, body Z = U10 ...
    .sign = {1.0f, 1.0f, -1.0f}, // ... with Z flipped: body Z = -U10
};

const sensor_axis_map_t GYRO_SENSOR_TO_BODY = {
    .axis = {1,     0,     2},     // body X = gyro Y, body Y = gyro X, body Z = gyro Z ...
    .sign = {-1.0f, -1.0f, -1.0f}, // ... all flipped: (-gy, -gx, -gz)
};
// clang-format on

static void to_body(const sensor_axis_map_t* m, const float* sensor, float* body_vec) {
    for (int i = 0; i < 3; i++) {
        body_vec[i] = m->sign[i] * sensor[m->axis[i]];
    }
}

static void to_sensor(const sensor_axis_map_t* m, const float* body_vec, float* sensor) {
    for (int i = 0; i < 3; i++) {
        sensor[m->axis[i]] = m->sign[i] * body_vec[i];
    }
}

void mag_sensor_to_body(const float* sensor, float* body_vec) { to_body(&MAG_SENSOR_TO_BODY, sensor, body_vec); }
void gyro_sensor_to_body(const float* sensor, float* body_vec) { to_body(&GYRO_SENSOR_TO_BODY, sensor, body_vec); }
void mag_body_to_sensor(const float* body_vec, float* sensor) { to_sensor(&MAG_SENSOR_TO_BODY, body_vec, sensor); }
void gyro_body_to_sensor(const float* body_vec, float* sensor) { to_sensor(&GYRO_SENSOR_TO_BODY, body_vec, sensor); }
