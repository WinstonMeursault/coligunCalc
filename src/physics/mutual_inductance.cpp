/**
 * @file mutual_inductance.cpp
 * @brief Mutual inductance and spatial gradient at filament and coil level.
 * @author Winston Meursault
 *
 * Filament-level: closed-form elliptic integrals with LRU caching.
 * Coil-level: 4D Gauss-Legendre tensor-product quadrature.
 *
 * @see NumericalModel Sec.4.2, Sec.4.3, Sec.4.4.
 */

#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/cache.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/elliptic.hpp"
#include "coilgun/physics/quadrature.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <tuple>

namespace coilgun::physics {

namespace {

struct TupleHash {
    std::size_t operator()(const std::tuple<double, double, double>& t) const {
        std::size_t h = 0;
        auto combine = [&h](double v) {
            h ^= std::hash<double>{}(v) + 0x9e3779b97f4a7c15ULL
                + (h << 6) + (h >> 2);
        };
        combine(std::get<0>(t));
        combine(std::get<1>(t));
        combine(std::get<2>(t));
        return h;
    }
};

using M_cache = LRUCache<std::tuple<double, double, double>, double, 4096,
                          TupleHash>;

static M_cache m_cache;
static M_cache grad_cache;

constexpr double k_min = 1e-12;
constexpr double k_max = 1.0 - 1e-12;

} // anonymous namespace

double mutual_inductance_filament(double radius_a, double radius_b,
                                  double separation) {
    auto key = std::make_tuple(radius_a, radius_b, separation);

    double cached;
    if (m_cache.get(key, cached)) {
        return cached;
    }

    double k = elliptic_modulus(radius_a, radius_b, separation);
    k = std::clamp(k, k_min, k_max);

    double m = k * k;
    double K = elliptic_k(m);
    double E = elliptic_e(m);

    double sqrt_ab = std::sqrt(radius_a * radius_b);
    double result = MU0 * sqrt_ab * ((2.0 / k - k) * K - (2.0 / k) * E);

    m_cache.put(key, result);
    return result;
}

double mutual_inductance_gradient_filament(double radius_a, double radius_b,
                                           double separation) {
    auto key = std::make_tuple(radius_a, radius_b, separation);

    double cached;
    if (grad_cache.get(key, cached)) {
        return cached;
    }

    if (std::abs(separation) < 1e-16) {
        grad_cache.put(key, 0.0);
        return 0.0;
    }

    double k = elliptic_modulus(radius_a, radius_b, std::abs(separation));
    k = std::clamp(k, k_min, k_max);

    double m = k * k;
    double one_minus_m = 1.0 - m;
    double K = elliptic_k(m);
    double E = elliptic_e(m);

    double sqrt_m = k;
    double sqrt_ab = std::sqrt(radius_a * radius_b);

    double bracket = 2.0 * one_minus_m * K - (2.0 - m) * E;
    double sign = (separation > 0.0) ? 1.0 : -1.0;
    double result = sign * MU0 * sqrt_m * std::abs(separation)
                    / (4.0 * one_minus_m * sqrt_ab) * bracket;

    grad_cache.put(key, result);
    return result;
}

// ---- coil-level 4D Gauss-Legendre tensor-product integration ----

namespace {

// Map [-1,1] to physical radial / axial coordinate.
inline double map_coord(double centre, double half_width, double xi) {
    return centre + half_width * xi;
}

// 4D tensor-product Gauss-Legendre quadrature on the unit hypercube.
template <typename F>
double integrate_4d(const F& integrand, int n_nodes = 9) {
    const auto& gl = gauss_legendre(n_nodes);
    double result = 0.0;

    for (std::size_t i1 = 0; i1 < gl.nodes.size(); ++i1) {
        const double w1 = gl.weights[i1];
        for (std::size_t j1 = 0; j1 < gl.nodes.size(); ++j1) {
            const double w2 = w1 * gl.weights[j1];
            for (std::size_t i2 = 0; i2 < gl.nodes.size(); ++i2) {
                const double w3 = w2 * gl.weights[i2];
                for (std::size_t j2 = 0; j2 < gl.nodes.size(); ++j2) {
                    result += w3 * gl.weights[j2] * integrand(
                        gl.nodes[i1], gl.nodes[j1],
                        gl.nodes[i2], gl.nodes[j2]);
                }
            }
        }
    }
    return result;
}

} // anonymous namespace

double mutual_inductance_coil(double rai, double rae, double la, int na,
                               double rbi, double rbe, double lb, int nb,
                               double separation, int n_nodes) {
    const double ra_mid = 0.5 * (rae + rai);
    const double ra_half = 0.5 * (rae - rai);
    const double rb_mid = 0.5 * (rbe + rbi);
    const double rb_half = 0.5 * (rbe - rbi);
    const double la_half = 0.5 * la;
    const double lb_half = 0.5 * lb;

    // Prefactor (Na * Nb / 16) per NumericalModel Eq.4.13.
    // The coordinate Jacobian is absorbed into this factor.

    auto kernel = [&](double r1, double z1, double r2, double z2) -> double {
        const double ra = map_coord(ra_mid, ra_half, r1);
        const double rb = map_coord(rb_mid, rb_half, r2);
        const double za = map_coord(0.0, la_half, z1);
        const double zb = map_coord(separation, lb_half, z2);
        return mutual_inductance_filament(ra, rb, std::abs(zb - za));
    };

    return (na * nb / 16.0) * integrate_4d(kernel, n_nodes);
}

double mutual_inductance_gradient_coil(double rai, double rae, double la,
                                        int na, double rbi, double rbe,
                                        double lb, int nb, double separation,
                                        int n_nodes) {
    const double ra_mid = 0.5 * (rae + rai);
    const double ra_half = 0.5 * (rae - rai);
    const double rb_mid = 0.5 * (rbe + rbi);
    const double rb_half = 0.5 * (rbe - rbi);
    const double la_half = 0.5 * la;
    const double lb_half = 0.5 * lb;

    auto kernel = [&](double r1, double z1, double r2, double z2) -> double {
        const double ra = map_coord(ra_mid, ra_half, r1);
        const double rb = map_coord(rb_mid, rb_half, r2);
        const double za = map_coord(0.0, la_half, z1);
        const double zb = map_coord(separation, lb_half, z2);
        return mutual_inductance_gradient_filament(ra, rb, zb - za);
    };

    return (na * nb / 16.0) * integrate_4d(kernel, n_nodes);
}

} // namespace coilgun::physics
