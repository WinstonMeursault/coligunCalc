/**
 * @file elliptic.hpp
 * @brief Complete elliptic integral wrappers and related utilities.
 * @author Winston Meursault
 */

#pragma once

namespace coilgun::physics {

/**
 * @brief Complete elliptic integral of the first kind K(m).
 * @param m Elliptic parameter m = k^2, 0 <= m < 1.
 * @return K(m) value.
 *
 * Internally calls boost::math::ellint_1(sqrt(m)).
 */
double elliptic_k(double m);

/**
 * @brief Complete elliptic integral of the second kind E(m).
 * @param m Elliptic parameter m = k^2, 0 <= m < 1.
 * @return E(m) value.
 *
 * Internally calls boost::math::ellint_2(sqrt(m)).
 */
double elliptic_e(double m);

/**
 * @brief Elliptic modulus k for two coaxial circular loops.
 * @param radius_a Radius of the first loop, m.
 * @param radius_b Radius of the second loop, m.
 * @param separation Axial separation between loop planes, m.
 * @return Modulus k (not m = k^2), clamped to avoid singularities.
 *
 * Eq. (4.7) in NumericalModel.md:
 * @f$ k = \sqrt{\frac{4ab}{(a+b)^2 + h^2}} @f$
 */
double elliptic_modulus(double radius_a, double radius_b, double separation);

} // namespace coilgun::physics
