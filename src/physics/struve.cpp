/**
 * @file struve.cpp
 * @brief Struve function implementations @f$ H_0(x) @f$, @f$ H_1(x) @f$.
 * @author Winston Meursault
 *
 * Three-regime strategy (power series, asymptotic, cross-checked) matching
 * SciPy's approach.
 *
 * @see NumericalModel Sec.4.1, Eq.(4.3)-(4.5).
 */

#include "coilgun/physics/struve.hpp"
#include "coilgun/physics/quadrature.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <algorithm>
#include <cmath>

namespace coilgun::physics {

namespace {

constexpr double k_pi = 3.14159265358979323846;

// ---- power-series branch ----

// H0(x) = (2/pi) * sum_{n=0}^{inf} (-1)^n * x^(2n+1) / (1^2 * 3^2 * ... * (2n+1)^2)
// Term recurrence: term_{n+1} = -term_n * x^2 / (2n+3)^2
//
// Numerically stable for x < ~25.  Beyond that individual terms overflow
// double before the alternating series begins to cancel.
constexpr int    k_max_series_terms        = 2000;
constexpr double k_series_tol             = 1e-16;

double struve_h0_power_series(double x) {
    double term = 2.0 * x / k_pi;   // n = 0
    double sum  = term;

    for (int n = 1; n < k_max_series_terms; ++n) {
        const double denom = 2.0 * n + 1.0;
        term *= -x * x / (denom * denom);
        sum  += term;

        if (std::abs(term) <= k_series_tol * std::abs(sum))
            break;
    }
    return sum;
}

// H1(x) = (2x^2/pi) * sum_{n=0}^{inf} (-1)^n * x^(2n) / ((2n+1)(2n+3))^2
double struve_h1_power_series(double x) {
    const double x2 = x * x;
    double term = 2.0 * x2 / (3.0 * k_pi);   // n = 0: 2x^2/(pi*3)
    double sum  = term;

    for (int n = 1; n < k_max_series_terms; ++n) {
        const double d1 = 2.0 * n + 1.0;
        const double d2 = 2.0 * n + 3.0;
        term *= -x2 / (d1 * d2);
        sum  += term;

        if (std::abs(term) <= k_series_tol * std::abs(sum))
            break;
    }
    return sum;
}

// ---- asymptotic branch via K_0 / K_1 integrals ----
//
// K_ν(z) ≡ H_ν(z) - Y_ν(z) (DLMF 11.2.5).  The K integrals have
// non-oscillatory, exponentially decaying integrands suitable for
// Gauss-Laguerre quadrature.
//
// Reliable for x >= 0.5; below that the integrand becomes singular at
// u=0 and Gauss-Laguerre loses precision.

// K0(x) = (2/pi) * integral_0^inf  e^{-u} / sqrt(u^2 + x^2)  du
double k0_integral(double x) {
    const auto& gl = gauss_laguerre_cached(30);
    double result = 0.0;
    for (std::size_t i = 0; i < gl.nodes.size(); ++i) {
        const double u = gl.nodes[i];
        result += gl.weights[i] / std::sqrt(u * u + x * x);
    }
    return (2.0 / k_pi) * result;
}

// K1(x) = (2/pi) * integral_0^inf  e^{-u} * sqrt(1 + (u/x)^2)  du
double k1_integral(double x) {
    const auto& gl = gauss_laguerre_cached(30);
    double result = 0.0;
    for (std::size_t i = 0; i < gl.nodes.size(); ++i) {
        const double u = gl.nodes[i];
        result += gl.weights[i] * std::sqrt(1.0 + (u * u) / (x * x));
    }
    return (2.0 / k_pi) * result;
}

double struve_h0_asymptotic(double x) {
    using boost::math::cyl_neumann;
    return cyl_neumann(0, x) + k0_integral(x);
}

double struve_h1_asymptotic(double x) {
    using boost::math::cyl_neumann;
    return cyl_neumann(1, x) + k1_integral(x);
}

} // anonymous namespace

// ---- public API ----
//
// H0 is odd (DLMF 11.2.1: z^{0+1}=z), H1 is even (z^{1+1}=z^2).
//
// Three-regime strategy for positive x:
//   x <  8  — power series only (asymptotic K integrals lose precision)
//   x >= 20 — asymptotic only (power series suffers cancellation)
//   8 <= x < 20 — compute both; if they agree to 1e-6 rel, use PS;
//                 otherwise asymptotic has the last word.

constexpr double k_ps_reliable_xmax  = 20.0;
constexpr double k_asym_reliable_xmin = 8.0;
constexpr double k_agreement_tol      = 1e-6;

double struve_h0(double x) {
    if (x < 0.0)
        return -struve_h0(-x);       // odd function
    if (x == 0.0) return 0.0;

    if (x >= k_ps_reliable_xmax)
        return struve_h0_asymptotic(x);

    const double ps = struve_h0_power_series(x);

    if (x < k_asym_reliable_xmin)
        return ps;

    const double asym = struve_h0_asymptotic(x);
    const double rel_diff = std::abs(ps - asym)
        / std::max(std::abs(ps), 1e-300);

    return (rel_diff < k_agreement_tol) ? ps : asym;
}

double struve_h1(double x) {
    x = std::abs(x);                 // even function
    if (x == 0.0) return 0.0;

    if (x >= k_ps_reliable_xmax)
        return struve_h1_asymptotic(x);

    const double ps = struve_h1_power_series(x);

    if (x < k_asym_reliable_xmin)
        return ps;

    const double asym = struve_h1_asymptotic(x);
    const double rel_diff = std::abs(ps - asym)
        / std::max(std::abs(ps), 1e-300);

    return (rel_diff < k_agreement_tol) ? ps : asym;
}

} // namespace coilgun::physics
