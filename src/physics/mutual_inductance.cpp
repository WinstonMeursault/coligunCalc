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

mutual_detail::MutualPairResult compute_filament_pair(double radius_a,
                                               double radius_b,
                                               double separation) {
    const double abs_separation = std::abs(separation);
    double k = elliptic_modulus(radius_a, radius_b, abs_separation);
    k = std::clamp(k, k_min, k_max);
    const double m = k * k;
    const double one_minus_m = 1.0 - m;
    const double K = elliptic_k(m);
    const double E = elliptic_e(m);
    const double sqrt_ab = std::sqrt(radius_a * radius_b);
    const double mutual = MU0 * sqrt_ab
        * ((2.0 / k - k) * K - (2.0 / k) * E);

    double gradient = 0.0;
    if (abs_separation >= 1e-16) {
        const double bracket = 2.0 * one_minus_m * K - (2.0 - m) * E;
        const double sign = separation > 0.0 ? 1.0 : -1.0;
        gradient = sign * MU0 * k * abs_separation
                   / (4.0 * one_minus_m * sqrt_ab) * bracket;
    }
    return {mutual, gradient};
}

mutual_detail::MutualPairResult filament_pair(double radius_a, double radius_b,
                                       double separation, bool use_cache) {
    if (!use_cache) return compute_filament_pair(radius_a, radius_b, separation);

    const auto key = std::make_tuple(radius_a, radius_b, separation);
    double cached_mutual = 0.0;
    double cached_gradient = 0.0;
    const bool mutual_hit = m_cache.get(key, cached_mutual);
    const bool gradient_hit = grad_cache.get(key, cached_gradient);
    if (mutual_hit && gradient_hit) return {cached_mutual, cached_gradient};

    const auto pair = compute_filament_pair(radius_a, radius_b, separation);
    if (!mutual_hit) m_cache.put(key, pair.mutual);
    if (!gradient_hit) grad_cache.put(key, pair.gradient);
    return pair;
}

} // anonymous namespace

double mutual_inductance_filament(double radius_a, double radius_b,
                                  double separation) {
    return mutual_inductance_filament(radius_a, radius_b, separation, false);
}

double mutual_inductance_filament(double radius_a, double radius_b,
                                  double separation, bool use_cache) {
    return filament_pair(radius_a, radius_b, separation, use_cache).mutual;
}

double mutual_inductance_gradient_filament(double radius_a, double radius_b,
                                           double separation) {
    return mutual_inductance_gradient_filament(radius_a, radius_b, separation, false);
}

double mutual_inductance_gradient_filament(double radius_a, double radius_b,
                                           double separation, bool use_cache) {
    return filament_pair(radius_a, radius_b, separation, use_cache).gradient;
}

// ---- coil-level 4D Gauss-Legendre tensor-product integration ----

namespace {

// Map [-1,1] to physical radial / axial coordinate.
inline double map_coord(double centre, double half_width, double xi) {
    return centre + half_width * xi;
}

template <typename F>
mutual_detail::MutualPairResult integrate_4d_pair(const F& integrand, int n_nodes) {
    const auto& gl = gauss_legendre_cached(n_nodes);
    mutual_detail::MutualPairResult result{0.0, 0.0};

    for (std::size_t i1 = 0; i1 < gl.nodes.size(); ++i1) {
        const double w1 = gl.weights[i1];
        for (std::size_t j1 = 0; j1 < gl.nodes.size(); ++j1) {
            const double w2 = w1 * gl.weights[j1];
            for (std::size_t i2 = 0; i2 < gl.nodes.size(); ++i2) {
                const double w3 = w2 * gl.weights[i2];
                for (std::size_t j2 = 0; j2 < gl.nodes.size(); ++j2) {
                    const auto pair = integrand(gl.nodes[i1], gl.nodes[j1],
                                                gl.nodes[i2], gl.nodes[j2]);
                    const double weight = w3 * gl.weights[j2];
                    result.mutual += weight * pair.mutual;
                    result.gradient += weight * pair.gradient;
                }
            }
        }
    }
    return result;
}

} // anonymous namespace

namespace mutual_detail {

MutualPairResult mutual_inductance_coil_pair(
        double rai, double rae, double la, int na,
        double rbi, double rbe, double lb, int nb,
        double separation, int n_nodes, bool use_cache) {
    const double ra_mid = 0.5 * (rae + rai);
    const double ra_half = 0.5 * (rae - rai);
    const double rb_mid = 0.5 * (rbe + rbi);
    const double rb_half = 0.5 * (rbe - rbi);
    const double la_half = 0.5 * la;
    const double lb_half = 0.5 * lb;

    auto kernel = [&](double r1, double z1, double r2, double z2) {
        const double ra = map_coord(ra_mid, ra_half, r1);
        const double rb = map_coord(rb_mid, rb_half, r2);
        const double za = map_coord(0.0, la_half, z1);
        const double zb = map_coord(separation, lb_half, z2);
        return filament_pair(ra, rb, zb - za, use_cache);
    };

    const auto integral = integrate_4d_pair(kernel, n_nodes);
    const double prefactor = (na * nb) / 16.0;
    return {prefactor * integral.mutual, prefactor * integral.gradient};
}

} // namespace mutual_detail

double mutual_inductance_coil(double rai, double rae, double la, int na,
                               double rbi, double rbe, double lb, int nb,
                               double separation, int n_nodes) {
    return mutual_inductance_coil(rai, rae, la, na, rbi, rbe, lb, nb,
                                  separation, n_nodes, false);
}

double mutual_inductance_coil(double rai, double rae, double la, int na,
                               double rbi, double rbe, double lb, int nb,
                               double separation, int n_nodes, bool use_cache) {
    return mutual_detail::mutual_inductance_coil_pair(
        rai, rae, la, na, rbi, rbe, lb, nb, separation, n_nodes, use_cache)
        .mutual;
}

double mutual_inductance_gradient_coil(double rai, double rae, double la,
                                        int na, double rbi, double rbe,
                                        double lb, int nb, double separation,
                                        int n_nodes) {
    return mutual_inductance_gradient_coil(rai, rae, la, na, rbi, rbe, lb, nb,
                                           separation, n_nodes, false);
}

double mutual_inductance_gradient_coil(double rai, double rae, double la,
                                        int na, double rbi, double rbe,
                                        double lb, int nb, double separation,
                                        int n_nodes, bool use_cache) {
    return mutual_detail::mutual_inductance_coil_pair(
        rai, rae, la, na, rbi, rbe, lb, nb, separation, n_nodes, use_cache)
        .gradient;
}

} // namespace coilgun::physics
