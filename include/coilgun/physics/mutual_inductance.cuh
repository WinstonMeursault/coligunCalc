/**
 * @file mutual_inductance.cuh
 * @brief Device-side filament-level mutual inductance (CUDA path).
 * @author Winston Meursault
 */

#pragma once

#ifdef __CUDACC__

#include "elliptic.cuh"
#include <algorithm>
#include <cmath>

namespace coilgun::physics {

constexpr double kMU0_device = 1.2566370614359173e-06; // 4*pi*1e-7

/**
 * @brief Mutual inductance between two coaxial loops — device-compatible.
 * @param radius_a Radius of the first loop, m.
 * @param radius_b Radius of the second loop, m.
 * @param separation Axial separation between loop planes, m.
 * @return Mutual inductance, H.
 */
__host__ __device__ inline double mutual_inductance_filament_device(
        double radius_a, double radius_b, double separation) {
    double k = elliptic_modulus(radius_a, radius_b, separation);
    k = fmax(1e-12, fmin(k, 1.0 - 1e-12));
    double m = k * k;
    double K = elliptic_k(m);
    double E = elliptic_e(m);
    double sqrt_ab = std::sqrt(radius_a * radius_b);
    return kMU0_device * sqrt_ab * ((2.0 / k - k) * K - (2.0 / k) * E);
}

/**
 * @brief Spatial gradient dM/dz of mutual inductance — device-compatible.
 * @param radius_a Radius of the first loop, m.
 * @param radius_b Radius of the second loop, m.
 * @param separation Axial separation between loop planes, m.
 * @return dM/dz gradient, H/m.
 */
__host__ __device__ inline double mutual_inductance_gradient_filament_device(
        double radius_a, double radius_b, double separation) {
    if (std::abs(separation) < 1e-16) return 0.0;
    double k = elliptic_modulus(radius_a, radius_b, std::abs(separation));
    k = fmax(1e-12, fmin(k, 1.0 - 1e-12));
    double m = k * k;
    double one_minus_m = 1.0 - m;
    double K = elliptic_k(m);
    double E = elliptic_e(m);
    double sqrt_m = k;
    double sqrt_ab = std::sqrt(radius_a * radius_b);
    double bracket = 2.0 * one_minus_m * K - (2.0 - m) * E;
    double sign = (separation > 0.0) ? 1.0 : -1.0;
    return sign * kMU0_device * sqrt_m * std::abs(separation)
           / (4.0 * one_minus_m * sqrt_ab) * bracket;
}

} // namespace coilgun::physics

#endif // __CUDACC__
