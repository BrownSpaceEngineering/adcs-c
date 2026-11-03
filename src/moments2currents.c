#include "include/moments2currents.h"
#include "arm_math.h"
#include <stddef.h>

/**
 * \fn moment2current3axis
 * \brief Converts desired moment into currents for 3 orthogonal magnetorquers. X and Y are solid rods, and Z is air-core
 *
 * Inputs:
 *  \param[in] m       : 3-element array, magnetic moment [mx, my, mz] [A·m2]
 *  \param[in] Imax    : OPTIONAL 3-element array of current limits [A]: Imax == NULL, use hardware defaults; Imax[i] < 0, no limit on that axis
 *  \param[out] I_out   : 3-element output array for coil currents [Ix, Iy, Iz] [A]
 */
void moment2current3axis(const float32_t *m, const float32_t *Imax, float32_t *I_out) {

    /* Linear torquer gains k = m / I [A·m2 per A] (from hardware team, subject to change)
     *   solid core: 2.6, air core: 1.02
     */

    const float32_t k[3] = {2.6f, 2.6f, 1.02f};   // X,Y Solid, Z Air
    static const float32_t Imax_default[3] = {0.156f, 0.156f, 0.104f};   // A, per hardware team (subject to change)
    
    /* Convert moments to currents: I = m / k */
    for (int i = 0; i < 3; i++) {
        I_out[i] = m[i] / k[i];
    }

    /* Imax == NULL, use hardware defaults; Imax[i] < 0, no limit on that axis */
    const float32_t *lim = (Imax != NULL) ? Imax : Imax_default;

    /* Saturation clamp */
    for (int i = 0; i < 3; i++) {
        if (lim[i] < 0.0f) {
            continue;
        }
        else if (I_out[i] >  lim[i]) {
            I_out[i] =  lim[i];
        }
        else if (I_out[i] < -lim[i]) {
            I_out[i] = -lim[i];
        }
    }
}
