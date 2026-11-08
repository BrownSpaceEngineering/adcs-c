#include "include/photodiode_determination.h"
#include "Include/dsp/matrix_functions.h"
#include "arm_math.h"
#include "math.h"

// #defines rather than `const int` so the arrays below aren't variable-length arrays
#define NUM_PAIRS 9
#define NUM_SELECTED_DIODES 5


const float MAX_READING = 1.7f;
const float MIN_READING = 0.17f;
const float PHOTODIODES[18][3] = {
    {  0.8660254f,  -0.27968387f, -0.41445981f },
    { -0.8660254f,   0.27968387f,  0.41445981f },

    {  0.8660254f,  -0.16748457f,  0.47111455f },
    { -0.8660254f,   0.16748457f, -0.47111455f },

    {  0.8660254f,   0.4999767f,   0.00482655f },
    { -0.8660254f,  -0.4999767f,  -0.00482655f },

    { -0.07387608f,  0.8660254f,   0.49451221f },
    {  0.07387608f, -0.8660254f,  -0.49451221f },

    { -0.48736734f,  0.8660254f,  -0.11168291f },
    {  0.48736734f, -0.8660254f,   0.11168291f },

    {  0.32616791f,  0.8660254f,  -0.37896503f },
    { -0.32616791f, -0.8660254f,   0.37896503f },

    { -0.36819159f, -0.33828235f, -0.8660254f  },
    {  0.36819159f,  0.33828235f,  0.8660254f  },

    {  0.48601695f,  0.1174203f,  -0.8660254f  },
    { -0.48601695f, -0.1174203f,   0.8660254f  },

    { -0.23581003f,  0.44090093f, -0.8660254f  },
    {  0.23581003f, -0.44090093f,  0.8660254f  }
};

// Assumes photodiode_readings are passed in the same order as PHOTODIODES.
// Pairs are 0-1, 2-3, ..., 16-17.
bool get_vec_from_photodiode_readings(const float* photodiode_readings,
                                      float* estimated_sun_vector) {
    int valid_readings = 0;
    int selected_indices[NUM_SELECTED_DIODES];
    float selected_readings[NUM_SELECTED_DIODES];

    for (int i = 0; i < NUM_SELECTED_DIODES; i++) {
        selected_indices[i] = -1;
        selected_readings[i] = -1.0f;
    }

    // For each pair, choose the brighter diode.
    // Keep only the NUM_SELECTED_DIODES brightest pair-winners
    for (int pair = 0; pair < NUM_PAIRS; pair++) {
        int i0 = 2 * pair;
        int i1 = 2 * pair + 1;

        int brighter_index;
        float brighter_reading;
        bool ok0 = isfinite(photodiode_readings[i0]);
        bool ok1 = isfinite(photodiode_readings[i1]);

        if (ok0 && ok1) {
            if (photodiode_readings[i0] >= photodiode_readings[i1]) {
                brighter_index = i0;
                brighter_reading = photodiode_readings[i0];
            } else {
                brighter_index = i1;
                brighter_reading = photodiode_readings[i1];
            }
        } else if (ok0 || ok1) {
            // One diode of the pair is bad. The pair faces opposite ways, so if the good one is lit
            // it must be the sunward one. If it's dark we can't tell (the bad one may be the lit
            // one), so drop the pair rather than feed a wrong zero into the fit.
            brighter_index = ok0 ? i0 : i1;
            brighter_reading = photodiode_readings[brighter_index];
            if (!(brighter_reading > MIN_READING)) continue;
        } else {
            continue; // both bad: skip the pair
        }

        if(brighter_reading > MIN_READING){
            valid_readings += 1;
        }

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

    if(valid_readings < 3){
        return false;
    }

    // Skipped pairs can leave empty slots; compact the selection and fit with what's there
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

        b[i] = selected_readings[i] / MAX_READING;
    }

    // Least-squares solve via normal equations: x = (A^T A)^-1 A^T b
    arm_matrix_instance_f32 A_mat = {(uint16_t)n_sel, 3, A};

    float32_t A_T[3 * NUM_SELECTED_DIODES];
    arm_matrix_instance_f32 A_T_mat = {3, (uint16_t)n_sel, A_T};
    arm_mat_trans_f32(&A_mat, &A_T_mat);

    float32_t AtA[3 * 3];
    arm_matrix_instance_f32 AtA_mat = {3, 3, AtA};
    arm_mat_mult_f32(&A_T_mat, &A_mat, &AtA_mat);

    float32_t AtA_inv[3 * 3];
    arm_matrix_instance_f32 AtA_inv_mat = {3, 3, AtA_inv};
    // (arm_mat_inverse_f32 overwrites AtA, which isn't needed afterwards)
    if (arm_mat_inverse_f32(&AtA_mat, &AtA_inv_mat) != ARM_MATH_SUCCESS) {
        return false; // selected normals are coplanar
    }

    float32_t Atb[3];
    arm_matrix_instance_f32 Atb_mat = {3, 1, Atb};
    arm_matrix_instance_f32 b_mat = {(uint16_t)n_sel, 1, b};
    arm_mat_mult_f32(&A_T_mat, &b_mat, &Atb_mat);

    arm_matrix_instance_f32 result_mat = {3, 1, estimated_sun_vector};
    arm_mat_mult_f32(&AtA_inv_mat, &Atb_mat, &result_mat);

    // Normalize output vector
    float mag = sqrtf(
        estimated_sun_vector[0] * estimated_sun_vector[0] +
        estimated_sun_vector[1] * estimated_sun_vector[1] +
        estimated_sun_vector[2] * estimated_sun_vector[2]
    );

    if (!(mag > 1e-6f) || !isfinite(mag)) {
        return false;
    }

    estimated_sun_vector[0] /= mag;
    estimated_sun_vector[1] /= mag;
    estimated_sun_vector[2] /= mag;

    return true;
}
