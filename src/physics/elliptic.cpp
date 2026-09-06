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
#include <limits>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

namespace coilgun::physics {

EllipticKe elliptic_ke(double m) {
    // Arithmetic-geometric mean (DLMF 19.8) on the modulus k = sqrt(m):
    // a_0 = 1, b_0 = sqrt(1 - m), c_0 = k; K(m) = pi / (2 a_inf) and
    // E(m) = K(m) * (1 - sum), sum = 0.5 c_0^2 + sum 2^(n-1) c_n^2.
    double a = 1.0;
    double b = std::sqrt(1.0 - m);
    double c = std::sqrt(m);
    double sum = 0.5 * c * c;
    double power = 0.5;
    for (int iteration = 0;
         iteration < 64 && c > std::numeric_limits<double>::epsilon() * a;
         ++iteration) {
        const double a_next = 0.5 * (a + b);
        const double b_next = std::sqrt(a * b);
        c = 0.5 * (a - b);
        a = a_next;
        b = b_next;
        power += power;
        sum += power * c * c;
    }
    const double first = 1.57079632679489661923 / a; // pi/2 / a_inf
    return {first, first * (1.0 - sum)};
}

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
