/**
 * @file struve.hpp
 * @brief Struve functions H0(x) and H1(x) for self-inductance computation.
 * @author Winston Meursault
 */

#pragma once

namespace coilgun::physics {

/**
 * @brief Struve function of order 0, H0(x).
 * @param x Argument (any real number).
 * @return H0(x) value.
 *
 * H0 is odd: H0(-x) = -H0(x) (DLMF 11.2.1, (z/2)^{0+1} = z factor).
 *
 * Three-regime strategy for |x|:
 * - |x| < 8: power series
 * - |x| >= 20: asymptotic via K0 integral + Y0 Bessel
 * - 8 <= |x| < 20: both methods computed, result cross-checked to 1e-6 rel
 */
double struve_h0(double x);

/**
 * @brief Struve function of order 1, H1(x).
 * @param x Argument (any real number).
 * @return H1(x) value.
 *
 * H1 is even: H1(-x) = H1(x) (DLMF 11.2.1, (z/2)^{1+1} = (z/2)^2 factor).
 *
 * Same three-regime strategy as H0, applied to |x|.
 */
double struve_h1(double x);

} // namespace coilgun::physics
