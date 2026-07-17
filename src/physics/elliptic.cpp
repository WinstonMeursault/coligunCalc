/**
 * @file elliptic.cpp
 * @brief Complete elliptic integral wrappers and geometric modulus.
 * @author Winston Meursault
 *
 * Wraps Boost.Math ellint_1 / ellint_2. Uses the @em parameter convention
 * (m = k^2), not the modulus convention.
 *
 * @see NumericalModel Sec.4.2, Eq.(4.7).
 */

#include "coilgun/physics/elliptic.hpp"

#include <algorithm>
#include <cmath>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

namespace coilgun::physics {

double elliptic_k(double m) {
    return boost::math::ellint_1(std::sqrt(m));
}

double elliptic_e(double m) {
    return boost::math::ellint_2(std::sqrt(m));
}

double elliptic_modulus(double radius_a, double radius_b, double separation) {
    double numerator   = 4.0 * radius_a * radius_b;
    double denominator = (radius_a + radius_b) * (radius_a + radius_b)
                       + separation * separation;
    double k = std::sqrt(numerator / denominator);

    // Clamp to avoid singularities in K/E evaluation
    k = std::clamp(k, 1e-12, 1.0 - 1e-10);

    return k;
}

} // namespace coilgun::physics
