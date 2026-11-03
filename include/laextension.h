#ifndef LAEXTENSION
#define LAEXTENSION
#include "arm_math.h"
#include <stdbool.h>

void cross(const float32_t* A, const float32_t* B, float32_t* C);
float32_t dot3(const float32_t* A, const float32_t* B);
float32_t min_arr(const float32_t* A, int elements);
float32_t max_arr(const float32_t* A, int elements);
float32_t trace(const arm_matrix_instance_f32* A);
float32_t l2_norm(const float32_t* A, int size);
// Normalizes a vector in place; returns false (and leaves it untouched) if its norm is ~0
bool normalize_vec(float32_t* A, int size);
bool all_finite(const float32_t* A, int size);
// Not "eye": PVDXos's src/linalg already defines a global eye()
void eye_f32(float32_t* A, int n);
#endif
