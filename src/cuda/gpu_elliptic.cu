/**
 * @file gpu_elliptic.cu
 * @brief CUDA kernel entry points for elliptic integral computation.
 * @author Winston Meursault
 */

#include "coilgun/physics/elliptic.cuh"

namespace coilgun::physics {

/**
 * @brief Kernel: compute K(m) and E(m) on device.
 * @param m Elliptic parameter.
 * @param out_k Output pointer for K(m).
 * @param out_e Output pointer for E(m).
 */
__global__ void elliptic_kernel(double m, double* out_k, double* out_e) {
    *out_k = elliptic_k(m);
    *out_e = elliptic_e(m);
}

/**
 * @brief Kernel: compute elliptic modulus k on device.
 * @param a Radius of first loop, m.
 * @param b Radius of second loop, m.
 * @param sep Axial separation, m.
 * @param out Output pointer for modulus k.
 */
__global__ void elliptic_modulus_kernel(double a, double b, double sep, double* out) {
    *out = elliptic_modulus(a, b, sep);
}

} // namespace coilgun::physics
