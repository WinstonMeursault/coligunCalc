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

template <typename T>
struct FilamentMutualResult {
    T mutual;
    T gradient;
};

__host__ __device__ inline FilamentMutualResult<double>
mutual_inductance_filament_pair_device(
        double radius_a, double radius_b, double separation) {
    const double abs_separation = std::abs(separation);
    double k = elliptic_modulus(radius_a, radius_b, abs_separation);
    k = fmax(1e-12, fmin(k, 1.0 - 1e-12));
    const double m = k * k;
    const double one_minus_m = 1.0 - m;
    const double K = elliptic_k(m);
    const double E = elliptic_e(m);
    const double sqrt_ab = std::sqrt(radius_a * radius_b);
    const double mutual = kMU0_device * sqrt_ab
        * ((2.0 / k - k) * K - (2.0 / k) * E);

    double gradient = 0.0;
    if (abs_separation >= 1e-16) {
        const double bracket = 2.0 * one_minus_m * K - (2.0 - m) * E;
        const double sign = separation > 0.0 ? 1.0 : -1.0;
        gradient = sign * kMU0_device * k * abs_separation
            / (4.0 * one_minus_m * sqrt_ab) * bracket;
    }
    return {mutual, gradient};
}

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

// --- FP32 filament mutual inductance (for GpuOptLevel::Aggressive) ---

constexpr float kMU0_f32 = 1.2566370614359173e-06f;

__host__ __device__ inline FilamentMutualResult<float>
mutual_inductance_filament_pair_f32(
        float radius_a, float radius_b, float separation) {
    const float abs_separation = fabsf(separation);
    const float numerator = 4.0f * radius_a * radius_b;
    const float radius_sum = radius_a + radius_b;
    const float denominator = radius_sum * radius_sum
        + abs_separation * abs_separation;
    float k = sqrtf(numerator / denominator);
    k = fmaxf(1e-12f, fminf(k, 1.0f - 1e-10f));
    const float m = k * k;
    const float one_minus_m = 1.0f - m;
    const float K = elliptic_k_f32(k);
    const float E = elliptic_e_f32(k);
    const float sqrt_ab = sqrtf(radius_a * radius_b);
    const float mutual = kMU0_f32 * sqrt_ab
        * ((2.0f / k - k) * K - (2.0f / k) * E);

    float gradient = 0.0f;
    if (abs_separation >= 1e-16f) {
        const float bracket = 2.0f * one_minus_m * K - (2.0f - m) * E;
        const float sign = separation > 0.0f ? 1.0f : -1.0f;
        gradient = sign * kMU0_f32 * k * abs_separation
            / (4.0f * one_minus_m * sqrt_ab) * bracket;
    }
    return {mutual, gradient};
}

/// FP32 mutual inductance between two coaxial loops — device-compatible.
__host__ __device__ inline float mutual_inductance_filament_f32(
        float radius_a, float radius_b, float separation) {
    float num = 4.0f * radius_a * radius_b;
    float den = (radius_a + radius_b) * (radius_a + radius_b) + separation * separation;
    float k = sqrtf(num / den);
    k = fmaxf(1e-12f, fminf(k, 1.0f - 1e-10f));
    float K = elliptic_k_f32(k);
    float E = elliptic_e_f32(k);
    float sqrt_ab = sqrtf(radius_a * radius_b);
    return kMU0_f32 * sqrt_ab * ((2.0f / k - k) * K - (2.0f / k) * E);
}

/// FP32 mutual inductance gradient — device-compatible.
__host__ __device__ inline float mutual_inductance_gradient_filament_f32(
        float radius_a, float radius_b, float separation) {
    if (fabsf(separation) < 1e-16f) return 0.0f;
    float abs_sep = fabsf(separation);
    float num = 4.0f * radius_a * radius_b;
    float den = (radius_a + radius_b) * (radius_a + radius_b) + abs_sep * abs_sep;
    float k = sqrtf(num / den);
    k = fmaxf(1e-12f, fminf(k, 1.0f - 1e-10f));
    float m = k * k;
    float one_minus_m = 1.0f - m;
    float K = elliptic_k_f32(k);
    float E = elliptic_e_f32(k);
    float sqrt_ab = sqrtf(radius_a * radius_b);
    float bracket = 2.0f * one_minus_m * K - (2.0f - m) * E;
    float sign = (separation > 0.0f) ? 1.0f : -1.0f;
    return sign * kMU0_f32 * k * abs_sep / (4.0f * one_minus_m * sqrt_ab) * bracket;
}

} // namespace coilgun::physics

#endif // __CUDACC__
