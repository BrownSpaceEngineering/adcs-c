#include "include/laextension.h"
#include "Include/arm_math_types.h"
#include "Include/dsp/matrix_functions.h"
#include "Include/dsp/statistics_functions.h"
#include "math.h"
#include "arm_math.h"

/**
 * \fn cross
 *
 * \brief Compute the cross product of two 1x3 vectors. C may alias A or B.
 *
 * \param[in] A: First vector
 * \param[in] B: Second vector
 * \param[out] C: Resultant vector
 */
void cross(const float32_t* A, const float32_t* B, float32_t* C) {
    float32_t c0 = A[1] * B[2] - A[2] * B[1];
    float32_t c1 = A[2] * B[0] - A[0] * B[2];
    float32_t c2 = A[0] * B[1] - A[1] * B[0];
    C[0] = c0;
    C[1] = c1;
    C[2] = c2;
}

/**
 * \fn dot3
 * \brief Dot product of two 1x3 vectors
 */
float32_t dot3(const float32_t* A, const float32_t* B) {
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

/**
 * \fn max_arr
 *
 * \brief Find the maximum value in an array (0 for an empty array)
 */
float32_t max_arr(const float32_t* A, int elements) {
    float32_t max;
    if (elements <= 0) {
        return 0;
    }
    arm_max_no_idx_f32(A, (uint32_t)elements, &max);
    return max;
}

/**
 * \fn min_arr
 *
 * \brief returns the minimum value in a flat array (0 for an empty array)
 */
float32_t min_arr(const float32_t* A, int elements) {
    float32_t min;
    if (elements <= 0) {
        return 0;
    }
    arm_min_no_idx_f32(A, (uint32_t)elements, &min);
    return min;
}

/**
 * \fn trace
 *
 * \brief Compute the trace of a square matrix (0 if the matrix is not square)
 */
float32_t trace(const arm_matrix_instance_f32* A) {
    if (A->numRows != A->numCols) {
        return 0;
    }

    int sz = A->numRows;
    const float32_t* data = A->pData;
    float32_t tr = 0;

    for (int i = 0; i < sz; i++) {
        tr += data[i * sz + i];
    }

    return tr;
}

/**
 * \fn l2_norm
 *
 * \brief Compute the L2 norm of a vector
 */
float32_t l2_norm(const float32_t* A, int n) {
    float32_t sum_squares = 0.0f;
    for (int i = 0; i < n; i++) {
        sum_squares += A[i] * A[i];
    }
    return sqrtf(sum_squares);
}

/**
 * \fn normalize_vec
 *
 * \brief Normalizes A in place. Returns false and leaves A untouched if |A| is ~0 or not finite.
 */
bool normalize_vec(float32_t* A, int n) {
    float32_t norm = l2_norm(A, n);
    if (!(norm > 1e-30f) || !isfinite(norm)) {
        return false;
    }
    for (int i = 0; i < n; i++) {
        A[i] /= norm;
    }
    return true;
}

/**
 * \fn all_finite
 * \brief True if no element is NaN or +-inf
 */
bool all_finite(const float32_t* A, int n) {
    for (int i = 0; i < n; i++) {
        if (!isfinite(A[i])) {
            return false;
        }
    }
    return true;
}

/**
 * \fn eye_f32
 *
 * \brief Fills a square matrix with the identity matrix
 *
 * \param[out] A: Square matrix to fill, n x n
 * \param[in] n: Size of the matrix
 */
void eye_f32(float32_t* A, int n) {
    for (int i = 0; i < n * n; i++) {
        A[i] = 0.0f;
    }
    for (int i = 0; i < n; i++) {
        A[i * n + i] = 1.0f;
    }
}
