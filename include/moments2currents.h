#ifndef MOMENTS2CURRENTS
#define MOMENTS2CURRENTS
#include "arm_math.h"
// I = m / k per axis, clamped to +-Imax[i]. Imax == NULL uses hardware defaults; negative Imax[i] means "no limit".
void moment2current3axis(const float32_t* m, const float32_t* Imax, float32_t* I_out);
#endif
