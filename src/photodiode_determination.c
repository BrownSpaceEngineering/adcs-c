
#include "include/photodiode_determination.h"
#include "Include/dsp/matrix_functions.h"
#include "arm_math.h"
#include "math.h"

#define NUM_PAIRS 9
#define NUM_SELECTED_DIODES 5
#define BROKEN_DIODE_INDEX 10

// Nominal brightness for generating test measurements; the fit uses raw readings.
const float MAX_READING = 1.7f;
const float MIN_READING = 0.17f;

// Photodiode normals directly from CAD
// Index 0 = physical pin 1, index 17 = physical pin 18.
const float PHOTODIODES[18][3] = {
    {-0.41445981f, 0.27968387f, 0.86602540f},   // Pin 1  (B)
    {0.86602540f, -0.48736734f, -0.11168291f},  // Pin 2  (R)
    {0.86602540f, 0.07387881f, -0.49451180f},   // Pin 3  (Q)
    {0.86602540f, -0.32616791f, 0.37896503f},   // Pin 4  (P)
    {0.41445981f, -0.27968387f, -0.86602540f},  // Pin 5  (D)
    {-0.36819159f, 0.86602540f, -0.33828235f},  // Pin 6  (J)
    {0.48601695f, 0.86602540f, 0.11742030f},    // Pin 7  (K)
    {-0.11168291f, -0.86602540f, -0.48736734f}, // Pin 8  (H)
    {-0.00482655f, 0.49997670f, -0.86602540f},  // Pin 9  (E)
    {-0.47111455f, -0.16748457f, -0.86602540f}, // Pin 10 (F)
    {-0.37896503f, -0.86602540f, 0.32616791f},  // Pin 11 (G, broken)
    {-0.86602540f, -0.33828235f, -0.36819159f}, // Pin 12 (M)
    {-0.86602540f, 0.11742030f, 0.48601695f},   // Pin 13 (N)
    {-0.86602540f, 0.44090093f, -0.23581003f},  // Pin 14 (O)
    {0.49451221f, -0.86602540f, -0.07387608f},  // Pin 15 (I)
    {-0.23581003f, 0.86602540f, 0.44090093f},   // Pin 16 (L)
    {0.00482655f, -0.49997670f, 0.86602540f},   // Pin 17 (A)
    {0.47111455f, 0.16748457f, 0.86602540f}     // Pin 18 (C)
};

// Closest-to-opposite photodiode pairs
// Pairs are chosen to maximize total angular separation,
// with each diode appearing in exactly one pair.
const int DIODE_PAIRS[NUM_PAIRS][2] = {
    {0, 4},  // Pins 1, 5
    {8, 16}, // Pins 9, 17
    {9, 17}, // Pins 10, 18
    {2, 12}, // Pins 3, 13
    {1, 13}, // Pins 2, 14
    {7, 15}, // Pins 8, 16
    {5, 14}, // Pins 6, 15
    {6, 10}, // Pins 7, 11 (pin 11 broken)
    {3, 11}  // Pins 4, 12
};

// Assumes readings are in physical pin order (1 through 18).
bool get_vec_from_photodiode_readings(const float* photodiode_readings,
                                      float* estimated_sun_vector) {
    int valid_readings = 0;
    int selected_indices[NUM_SELECTED_DIODES];
    float selected_readings[NUM_SELECTED_DIODES];

    for (int i = 0; i < NUM_SELECTED_DIODES; i++) {
        selected_indices[i] = -1;
        selected_readings[i] = -1.0f;
    }

    // For each pair, choose the brighter functioning diode.
    // Keep only the NUM_SELECTED_DIODES brightest pair-winners.
    for (int pair = 0; pair < NUM_PAIRS; pair++) {
        int i0 = DIODE_PAIRS[pair][0];
        int i1 = DIODE_PAIRS[pair][1];

        int brighter_index;
        float brighter_reading;

        bool ok0 = (i0 != BROKEN_DIODE_INDEX) && isfinite(photodiode_readings[i0]);
        bool ok1 = (i1 != BROKEN_DIODE_INDEX) && isfinite(photodiode_readings[i1]);

        if (ok0 && ok1) {
            if (photodiode_readings[i0] >= photodiode_readings[i1]) {
                brighter_index = i0;
                brighter_reading = photodiode_readings[i0];
            } else {
                brighter_index = i1;
                brighter_reading = photodiode_readings[i1];
            }
        } else if (ok0 || ok1) {
            // Only one diode is functional.
            // Allow it to compete for the top readings.
            brighter_index = ok0 ? i0 : i1;
            brighter_reading = photodiode_readings[brighter_index];
        } else {
            continue; // Both diodes unavailable
        }

        // Reject weak readings so dark noise cannot produce a false sun vector.
        if (!(brighter_reading > MIN_READING)) {
            continue;
        }

        valid_readings += 1;

        int min_pos = 0;
        for (int i = 1; i < NUM_SELECTED_DIODES; i++) {
            if (selected_readings[i] < selected_readings[min_pos]) {
                min_pos = i;
            }
        }

        if (brighter_reading > selected_readings[min_pos]) {
            selected_readings[min_pos] = brighter_reading;
            selected_indices[min_pos] = brighter_index;
        }
    }

    if (valid_readings < 3) {
        return false;
    }

    // Skipped pairs can leave empty slots; compact the selection.
    int n_sel = 0;
    for (int i = 0; i < NUM_SELECTED_DIODES; i++) {
        if (selected_indices[i] >= 0) {
            selected_indices[n_sel] = selected_indices[i];
            selected_readings[n_sel] = selected_readings[i];
            n_sel++;
        }
    }

    if (n_sel < 3) {
        return false;
    }

    // Construct system of equations
    float32_t A[NUM_SELECTED_DIODES * 3];
    float32_t b[NUM_SELECTED_DIODES];

    for (int i = 0; i < n_sel; i++) {
        int diode_index = selected_indices[i];

        A[i * 3 + 0] = PHOTODIODES[diode_index][0];
        A[i * 3 + 1] = PHOTODIODES[diode_index][1];
        A[i * 3 + 2] = PHOTODIODES[diode_index][2];

        b[i] = selected_readings[i];
    }

    // Least-squares solve via normal equations:
    // x = (A^T A)^-1 A^T b
    arm_matrix_instance_f32 A_mat = {(uint16_t)n_sel, 3, A};

    float32_t A_T[3 * NUM_SELECTED_DIODES];
    arm_matrix_instance_f32 A_T_mat = {3, (uint16_t)n_sel, A_T};
    arm_mat_trans_f32(&A_mat, &A_T_mat);

    float32_t AtA[3 * 3];
    arm_matrix_instance_f32 AtA_mat = {3, 3, AtA};
    arm_mat_mult_f32(&A_T_mat, &A_mat, &AtA_mat);

    float32_t AtA_inv[3 * 3];
    arm_matrix_instance_f32 AtA_inv_mat = {3, 3, AtA_inv};

    if (arm_mat_inverse_f32(&AtA_mat, &AtA_inv_mat) != ARM_MATH_SUCCESS) {
        return false;
    }

    float32_t Atb[3];
    arm_matrix_instance_f32 Atb_mat = {3, 1, Atb};
    arm_matrix_instance_f32 b_mat = {(uint16_t)n_sel, 1, b};
    arm_mat_mult_f32(&A_T_mat, &b_mat, &Atb_mat);

    arm_matrix_instance_f32 result_mat = {3, 1, estimated_sun_vector};
    arm_mat_mult_f32(&AtA_inv_mat, &Atb_mat, &result_mat);

    // Normalize output vector
    float mag = sqrtf(estimated_sun_vector[0] * estimated_sun_vector[0]
                      + estimated_sun_vector[1] * estimated_sun_vector[1]
                      + estimated_sun_vector[2] * estimated_sun_vector[2]);

    if (!(mag > 1e-6f) || !isfinite(mag)) {
        return false;
    }

    estimated_sun_vector[0] /= mag;
    estimated_sun_vector[1] /= mag;
    estimated_sun_vector[2] /= mag;

    return true;
}
