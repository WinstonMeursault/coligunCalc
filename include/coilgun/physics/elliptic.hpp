/**
 * @file elliptic.hpp
 * @brief Complete elliptic integral wrappers and related utilities.
 * @author Winston Meursault
 */

#pragma once

namespace coilgun::physics {

/**
 * @brief Complete elliptic integrals of the first and second kind together.
 */
struct EllipticKe {
    double k; ///< K(m), complete elliptic integral of the first kind.
    double e; ///< E(m), complete elliptic integral of the second kind.
};

/**
 * @brief Complete elliptic integrals K(m) and E(m) from one AGM evaluation.
 * @param m Elliptic parameter m = k^2, 0 <= m < 1 (same convention as
 *          elliptic_k/elliptic_e).
 * @return Struct with `k` = K(m) and `e` = E(m).
 *
 * Uses the arithmetic-geometric mean (DLMF 19.8), which produces both
 * integrals from one shared iteration at roughly the cost of a single
 * ellint_1 evaluation. Results agree with the Boost wrappers to within a few
 * ulps; prefer elliptic_k/elliptic_e when a bit-exact Boost reference is
 * required.
 */
EllipticKe elliptic_ke(double m);

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
