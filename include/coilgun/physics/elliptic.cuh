/**
 * @file elliptic.cuh
 * @brief Device-side complete elliptic integrals (CUDA path).
 * @author Winston Meursault
 *
 * Only compiled when __CUDACC__ is defined.
 * Uses Boost.Math with BOOST_MATH_ENABLE_CUDA for __device__ support.
 */

#pragma once

#ifdef __CUDACC__

#include <algorithm>
#include <cmath>

#define BOOST_MATH_ENABLE_CUDA
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

namespace coilgun::physics {

/**
 * @brief Complete elliptic integral of the first kind K(m) — device-compatible.
 * @param m Elliptic parameter m = k^2, 0 <= m < 1.
 * @return K(m) value.
 */
__host__ __device__ inline double elliptic_k(double m) {
    return boost::math::ellint_1(std::sqrt(m));
}

/**
 * @brief Complete elliptic integral of the second kind E(m) — device-compatible.
 * @param m Elliptic parameter m = k^2, 0 <= m < 1.
 * @return E(m) value.
 */
__host__ __device__ inline double elliptic_e(double m) {
    return boost::math::ellint_2(std::sqrt(m));
}

/**
 * @brief Elliptic modulus k for two coaxial circular loops — device-compatible.
 * @param radius_a Radius of the first loop, m.
 * @param radius_b Radius of the second loop, m.
 * @param separation Axial separation between loop planes, m.
 * @return Modulus k (not m = k^2), clamped to avoid singularities.
 *
 * @f$ k = \sqrt{\frac{4ab}{(a+b)^2 + h^2}} @f$
 */
__host__ __device__ inline double elliptic_modulus(
        double radius_a, double radius_b, double separation) {
    double numerator   = 4.0 * radius_a * radius_b;
    double denominator = (radius_a + radius_b) * (radius_a + radius_b)
                       + separation * separation;
    double k = std::sqrt(numerator / denominator);
    k = fmax(1e-12, fmin(k, 1.0 - 1e-10));
    return k;
}

} // namespace coilgun::physics

#endif // __CUDACC__
