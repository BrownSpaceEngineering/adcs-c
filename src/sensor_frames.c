#include "include/sensor_frames.h"

// See the alignment table in sensor_frames.h.
// clang-format off
#define MAG_AS_BUILT {                                                                        \
    .axis = {0,    1,    2},     /* body X = U8, body Y = U9, body Z = U10 ... */           \
    .sign = {1.0f, 1.0f, -1.0f}, /* ... with Z flipped: body Z = -U10 */                    \
}
#define GYRO_AS_BUILT {                                                                       \
    .axis = {1,     0,     2},     /* body X = gyro Y, body Y = gyro X, body Z = gyro Z ... */ \
    .sign = {-1.0f, -1.0f, -1.0f}, /* ... all flipped: (-gy, -gx, -gz) */                    \
}
// clang-format on

const sensor_axis_map_t MAG_SENSOR_TO_BODY = MAG_AS_BUILT;
const sensor_axis_map_t GYRO_SENSOR_TO_BODY = GYRO_AS_BUILT;

// Active maps used by the flight path; start as built, replaced from the ground
static sensor_axis_map_t mag_active = MAG_AS_BUILT;
static sensor_axis_map_t gyro_active = GYRO_AS_BUILT;

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

void mag_sensor_to_body(const float* sensor, float* body_vec) { to_body(&mag_active, sensor, body_vec); }
void gyro_sensor_to_body(const float* sensor, float* body_vec) { to_body(&gyro_active, sensor, body_vec); }
void mag_body_to_sensor(const float* body_vec, float* sensor) { to_sensor(&MAG_SENSOR_TO_BODY, body_vec, sensor); }
void gyro_body_to_sensor(const float* body_vec, float* sensor) { to_sensor(&GYRO_SENSOR_TO_BODY, body_vec, sensor); }

/* ---------------- uplink codes ---------------- */

bool sensor_map_from_codes(const int8_t* codes, sensor_axis_map_t* out) {
    if (codes == 0 || out == 0) return false;
    sensor_axis_map_t m;
    bool used[3] = {false, false, false};
    for (int i = 0; i < 3; i++) {
        int a = codes[i] < 0 ? -codes[i] : codes[i];
        if (a < 1 || a > 3 || used[a - 1]) return false; // not a signed permutation
        used[a - 1] = true;
        m.axis[i] = a - 1;
        m.sign[i] = codes[i] < 0 ? -1.0f : 1.0f;
    }
    *out = m;
    return true;
}

void sensor_map_to_codes(const sensor_axis_map_t* m, int8_t* codes) {
    for (int i = 0; i < 3; i++) {
        int8_t c = (int8_t)(m->axis[i] + 1);
        codes[i] = m->sign[i] < 0.0f ? (int8_t)-c : c;
    }
}

// Valid: axis is a permutation of 0, 1, 2 and every sign is exactly +-1
static bool map_valid(const sensor_axis_map_t* m) {
    bool used[3] = {false, false, false};
    for (int i = 0; i < 3; i++) {
        int a = m->axis[i];
        if (a < 0 || a > 2 || used[a]) return false;
        used[a] = true;
        float s = m->sign[i] < 0.0f ? -m->sign[i] : m->sign[i];
        if (!(s > 0.999f && s < 1.001f)) return false; // exactly +-1 (written as a range: -Wfloat-equal)
    }
    return true;
}

float sensor_map_det(const sensor_axis_map_t* m) {
    // Sign of the permutation (count inversions) times the product of the signs
    int inversions = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            if (m->axis[i] > m->axis[j]) inversions++;
        }
    }
    float det = (inversions % 2) ? -1.0f : 1.0f;
    for (int i = 0; i < 3; i++) det *= m->sign[i];
    return det;
}

/* ---------------- active maps ---------------- */

bool mag_set_sensor_map(const sensor_axis_map_t* m) {
    if (m == 0 || !map_valid(m)) return false;
    mag_active = *m;
    for (int i = 0; i < 3; i++) mag_active.sign[i] = m->sign[i] < 0.0f ? -1.0f : 1.0f; // exactly +-1
    return true;
}

bool gyro_set_sensor_map(const sensor_axis_map_t* m) {
    if (m == 0 || !map_valid(m) || !(sensor_map_det(m) > 0.0f)) return false; // det +1 only
    gyro_active = *m;
    for (int i = 0; i < 3; i++) gyro_active.sign[i] = m->sign[i] < 0.0f ? -1.0f : 1.0f; // exactly +-1
    return true;
}

void mag_get_sensor_map(sensor_axis_map_t* m) { *m = mag_active; }
void gyro_get_sensor_map(sensor_axis_map_t* m) { *m = gyro_active; }

void sensor_maps_reset(void) {
    mag_active = MAG_SENSOR_TO_BODY;
    gyro_active = GYRO_SENSOR_TO_BODY;
}
