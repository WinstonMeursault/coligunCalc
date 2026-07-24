/**
 * @file self_inductance.cpp
 * @brief Self-inductance of hollow cylindrical coils (exact and table-based).
 * @author Winston Meursault
 *
 * Implements the Bessel/Struve kernel integration for T(q,p) and the
 * auto-select dispatch between table-lookup and exact computation.
 *
 * @see NumericalModel Sec.4.1, Eq.(4.1)-(4.6).
 */

#include "coilgun/physics/self_inductance.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/lookup_table_data.hpp"
#include "coilgun/physics/lookup_tables.hpp"
#include "coilgun/physics/quadrature.hpp"
#include "coilgun/physics/struve.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <cmath>

namespace coilgun::physics {

namespace {

constexpr double k_pi = 3.14159265358979323846;

// Bessel/Struve kernel  U(x)  for the T(q,p) integral.
// Python prototype closed form:
//   U(x) = pi * (-J1(x)*H0(x) + p*J1(px)*H0(px) + J0(x)*H1(x) - p*J0(px)*H1(px)) / (2*x^2)
//
// Returns 0 for x near zero where the formula is indeterminate.
double u_kernel(double x, double p) {
    if (x < 1e-14) return 0.0;

    using boost::math::cyl_bessel_j;
    const double px = p * x;
    const double numerator = k_pi * (
        -cyl_bessel_j(1, x) * struve_h0(x)
        + p * cyl_bessel_j(1, px) * struve_h0(px)
        + cyl_bessel_j(0, x) * struve_h1(x)
        - p * cyl_bessel_j(0, px) * struve_h1(px)
    );
    return numerator / (2.0 * x * x);
}

// Integrand  f(x) = U(x)^2 * (q*x + e^{-q*x} - 1)    (Eq.4.3)
double t_integrand(double x, double q, double p) {
    const double u = u_kernel(x, p);
    if (!std::isfinite(u)) return 0.0;
    const double bracket = q * x + std::exp(-q * x) - 1.0;
    return u * u * bracket;
}

// Transform [0, 1] -> [0, infty):  x = (1 - u) / u,  dx = -du / u^2
// integral_0^infty f(x) dx = integral_0^1 f((1-u)/u) / u^2  du
double transformed_integrand(double u, double q, double p) {
    if (u <= 0.0 || u >= 1.0) return 0.0;
    const double x = (1.0 - u) / u;
    return t_integrand(x, q, p) / (u * u);
}

// Integrate T(q,p) = int_0^infty  U(x)^2 * (q*x + e^{-q*x} - 1)  dx
//
// We transform [0, infty) -> [0, 1) via x = (1-u)/u  and apply composite
// Gauss-Legendre on [0, 1].  The tail beyond the last sub-interval is
// estimated analytically as ~ q/(4 * x_last^4) from the asymptotic
// behaviour U(x) ~ 1/x^3.
double integrate_T_kernel(double q, double p) {
    const auto& gl = gauss_legendre_cached(16);
    const int n_sub = 2048;
    const double du = 1.0 / n_sub;

    double total = 0.0;

    for (int s = 0; s < n_sub; ++s) {
        const double a = s * du;
        const double b = a + du;
        const double mid = 0.5 * (a + b);
        const double half = 0.5 * du;

        double sub = 0.0;
        for (std::size_t i = 0; i < gl.nodes.size(); ++i) {
            const double u = mid + half * gl.nodes[i];
            sub += gl.weights[i] * transformed_integrand(u, q, p);
        }
        total += half * sub;

        if (std::abs(half * sub) < 1e-18 * std::abs(total) && s > 32)
            break;
    }

    return total;
}

} // anonymous namespace

// ---- public API ----

static double self_inductance_impl(double ri, double re, double length,
                                   double nc, const double T) {
    return 2.0 * k_pi * MU0 * nc * nc * std::pow(ri, 5) * T;
}

static bool in_table_range(double q, double p) {
    using detail::k_table_q0;
    using detail::k_table_p0;
    using detail::k_table_q_max;
    using detail::k_table_p_max;
    return q >= k_table_q0 && q <= k_table_q_max &&
           p >= k_table_p0 && p <= k_table_p_max;
}

double self_inductance(double ri, double re, double length, double nc,
                       bool force_exact) {
    const double p = re / ri;
    const double q = length / ri;

    if (!force_exact && in_table_range(q, p)) {
        return self_inductance_impl(ri, re, length, nc, inductance_shape_factor(q, p));
    }
    return self_inductance_impl(ri, re, length, nc, integrate_T_kernel(q, p));
}

double self_inductance_exact(double ri, double re, double length, double nc) {
    const double p = re / ri;
    const double q = length / ri;

    return self_inductance_impl(ri, re, length, nc, integrate_T_kernel(q, p));
}

double inductance_shape_factor_reference(double q, double p) {
    return integrate_T_kernel(q, p);
}

} // namespace coilgun::physics
