#include <doctest/doctest.h>
#include "coilgun/physics/mutual_inductance.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>

using coilgun::physics::mutual_inductance_filament;
using coilgun::physics::mutual_inductance_gradient_filament;
using coilgun::physics::mutual_inductance_coil;
using coilgun::physics::mutual_inductance_gradient_coil;
using coilgun::physics::mutual_detail::mutual_inductance_coil_pair;

namespace {

struct GoldenPair {
    double separation;
    int n_nodes;
    double mutual;
    double gradient;
};

// Values were generated from the pre-change implementation and are kept as
// independent regression oracles for the fused production path.
constexpr std::array<GoldenPair, 15> golden_pairs = {{
    { 0.00015,  4, 5.477318812775153e-07, -6.7419257289816249e-07},
    { 0.00015,  9, 5.0811225229173607e-07, -3.2631236370671433e-07},
    { 0.00015, 16, 5.0442376957113391e-07, -1.0014878721561817e-07},
    { 0.03,     4, 3.2755008138911007e-07, -9.3691675159691946e-06},
    { 0.03,     9, 3.4887016983446409e-07, -6.80520777077855e-06},
    { 0.03,    16, 3.4725384266923269e-07, -7.7234593242791517e-06},
    { 0.12,     4, 8.4682091412026196e-09, -2.4641585311625103e-07},
    { 0.12,     9, 8.4684218474156013e-09, -2.4643537073954618e-07},
    { 0.12,    16, 8.4684218474146153e-09, -2.4643537073925347e-07},
    {-0.03,     4, 3.2755008138911007e-07,  9.3691675159691982e-06},
    {-0.03,     9, 3.4887016983446431e-07,  6.8052077707785542e-06},
    {-0.03,    16, 3.4725384266923052e-07,  7.723459324279089e-06},
    {-0.12,     4, 8.4682091412026147e-09,  2.4641585311625093e-07},
    {-0.12,     9, 8.4684218474156063e-09,  2.4643537073954676e-07},
    {-0.12,    16, 8.4684218474146335e-09,  2.4643537073925268e-07},
}};

struct ReferenceErrorStats {
    double max_abs_mutual = 0.0;
    double max_abs_gradient = 0.0;
    double max_relative_mutual = 0.0;
    double max_relative_gradient = 0.0;
};

ReferenceErrorStats reference_errors;

GoldenPair golden_pair(double separation, int n_nodes) {
    for (const auto& golden : golden_pairs) {
        if (golden.separation == separation && golden.n_nodes == n_nodes)
            return golden;
    }
    throw std::logic_error("missing mutual-inductance golden value");
}

void check_pair_against_reference(double separation, int n_nodes,
                                  bool use_cache) {
    constexpr double rai = 0.015125;
    constexpr double rae = 0.029125;
    constexpr double la = 0.063;
    constexpr int na = 63;
    constexpr double rbi = 0.012;
    constexpr double rbe = 0.015;
    constexpr double lb = 0.07;
    constexpr int nb = 1;

    const auto expected = golden_pair(separation, n_nodes);
    const auto actual = mutual_inductance_coil_pair(
        rai, rae, la, na, rbi, rbe, lb, nb, separation, n_nodes, use_cache);
    const double mutual_error = std::abs(actual.mutual - expected.mutual);
    const double gradient_error = std::abs(actual.gradient - expected.gradient);
    reference_errors.max_abs_mutual = std::max(
        reference_errors.max_abs_mutual, mutual_error);
    reference_errors.max_abs_gradient = std::max(
        reference_errors.max_abs_gradient, gradient_error);
    reference_errors.max_relative_mutual = std::max(
        reference_errors.max_relative_mutual,
        mutual_error / std::max(std::abs(expected.mutual), 1e-300));
    reference_errors.max_relative_gradient = std::max(
        reference_errors.max_relative_gradient,
        gradient_error / std::max(std::abs(expected.gradient), 1e-300));
    CHECK(actual.mutual == doctest::Approx(expected.mutual).epsilon(1e-12));
    CHECK(actual.gradient == doctest::Approx(expected.gradient).epsilon(1e-12));

    CHECK(mutual_inductance_coil(
        rai, rae, la, na, rbi, rbe, lb, nb, separation, n_nodes, use_cache)
          == doctest::Approx(expected.mutual).epsilon(1e-12));
    CHECK(mutual_inductance_gradient_coil(
        rai, rae, la, na, rbi, rbe, lb, nb, separation, n_nodes, use_cache)
          == doctest::Approx(expected.gradient).epsilon(1e-12));
}

} // namespace

using coilgun::physics::mutual_inductance_filament;
using coilgun::physics::mutual_inductance_gradient_filament;

TEST_CASE("mutual_inductance_filament — symmetry") {
    CHECK(mutual_inductance_filament(0.02, 0.03, 0.005)
          == doctest::Approx(mutual_inductance_filament(0.03, 0.02, 0.005))
                 .epsilon(1e-12));
}

TEST_CASE("mutual_inductance_filament — reference values") {
    CHECK(mutual_inductance_filament(0.02, 0.03, 0.005)
          == doctest::Approx(2.959493865841560e-08).epsilon(1e-10));

    CHECK(mutual_inductance_filament(0.0135, 0.0135, 0.001)
          == doctest::Approx(4.557715786266084e-08).epsilon(1e-10));

    CHECK(mutual_inductance_filament(0.01, 0.01, 0.001)
          == doctest::Approx(3.002876303701474e-08).epsilon(1e-10));

    CHECK(mutual_inductance_filament(0.01, 0.01, 0.01)
          == doctest::Approx(4.940784630798257e-09).epsilon(1e-10));

    CHECK(mutual_inductance_filament(0.02, 0.01, 0.01)
          == doctest::Approx(6.987324633639440e-09).epsilon(1e-10));
}

TEST_CASE("mutual_inductance_gradient_filament — antisymmetry") {
    CHECK(mutual_inductance_gradient_filament(0.01, 0.01, 0.001)
          == doctest::Approx(
                 -mutual_inductance_gradient_filament(0.01, 0.01, -0.001))
                 .epsilon(1e-12));

    CHECK(mutual_inductance_gradient_filament(0.02, 0.01, 0.01)
          == doctest::Approx(
                 -mutual_inductance_gradient_filament(0.02, 0.01, -0.01))
                 .epsilon(1e-12));
}

TEST_CASE("mutual_inductance_gradient_filament — zero at d=0") {
    CHECK(mutual_inductance_gradient_filament(0.01, 0.01, 0.0)
          == doctest::Approx(0.0).epsilon(1e-12));

    CHECK(mutual_inductance_gradient_filament(0.02, 0.03, 0.0)
          == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("mutual_inductance_gradient_filament — reference values") {
    CHECK(mutual_inductance_gradient_filament(0.01, 0.01, 0.001)
          == doctest::Approx(-1.239937005403076e-05).epsilon(1e-10));

    CHECK(mutual_inductance_gradient_filament(0.02, 0.01, 0.01)
          == doctest::Approx(-5.079612386972773e-07).epsilon(1e-10));

    CHECK(mutual_inductance_gradient_filament(0.03, 0.015, 0.02)
          == doctest::Approx(-4.205063421913522e-07).epsilon(1e-10));
}

TEST_CASE("mutual_inductance_filament — cache hit") {
    double m1 = mutual_inductance_filament(0.025, 0.025, 0.003);
    double m2 = mutual_inductance_filament(0.025, 0.025, 0.003);
    CHECK(m1 == doctest::Approx(m2).epsilon(1e-15));
}

TEST_CASE("mutual_inductance_gradient_filament — cache hit") {
    double g1 = mutual_inductance_gradient_filament(0.025, 0.025, 0.003);
    double g2 = mutual_inductance_gradient_filament(0.025, 0.025, 0.003);
    CHECK(g1 == doctest::Approx(g2).epsilon(1e-15));
}

TEST_CASE("mutual_inductance_filament — M even in d") {
    CHECK(mutual_inductance_filament(0.01, 0.02, 0.005)
          == doctest::Approx(
                 mutual_inductance_filament(0.01, 0.02, -0.005))
                 .epsilon(1e-12));
}

TEST_CASE("mutual_inductance_coil_pair — returns fused values") {
    const auto pair = mutual_inductance_coil_pair(
        0.015125, 0.029125, 0.063, 63,
        0.012,    0.015,    0.07,   1,
        0.03, 4, false);

    CHECK(pair.mutual == doctest::Approx(mutual_inductance_coil(
        0.015125, 0.029125, 0.063, 63,
        0.012,    0.015,    0.07,   1,
        0.03, 4, false)).epsilon(1e-12));
    CHECK(pair.gradient == doctest::Approx(mutual_inductance_gradient_coil(
        0.015125, 0.029125, 0.063, 63,
        0.012,    0.015,    0.07,   1,
        0.03, 4, false)).epsilon(1e-12));
}

TEST_CASE("mutual_inductance_coil_pair — reference separation and node matrix") {
    for (const double separation : {0.00015, 0.03, 0.12, -0.03, -0.12}) {
        for (const int n_nodes : {4, 9, 16}) {
            check_pair_against_reference(separation, n_nodes, false);
            check_pair_against_reference(separation, n_nodes, true);
        }
    }
    std::cout << "B1-T1 max_abs_mutual=" << reference_errors.max_abs_mutual
              << " max_abs_gradient=" << reference_errors.max_abs_gradient
              << " max_rel_mutual=" << reference_errors.max_relative_mutual
              << " max_rel_gradient=" << reference_errors.max_relative_gradient
              << '\n';
}

TEST_CASE("mutual_inductance_coil_pair — cutoff boundary") {
    for (const double separation : {0.0, 0.5e-16, 2.0e-16, -2.0e-16}) {
        constexpr double expected_mutual = 1.0972358946947983e-08;
        const double expected_gradient = separation == 0.0
            || std::abs(separation) < 1e-16
            ? 0.0
            : (separation > 0.0 ? -2.5468514453368917e-20
                                : 2.5468514453368917e-20);
        const auto actual = mutual_inductance_coil_pair(
            0.02, 0.02, 0.0, 1,
            0.01, 0.01, 0.0, 1,
            separation, 4, false);
        CHECK(actual.mutual == doctest::Approx(expected_mutual).epsilon(1e-12));
        CHECK(actual.gradient
              == doctest::Approx(expected_gradient).epsilon(1e-12));
        CHECK(mutual_inductance_gradient_coil(
            0.02, 0.02, 0.0, 1,
            0.01, 0.01, 0.0, 1,
            separation, 4, true)
              == doctest::Approx(expected_gradient).epsilon(1e-12));
    }
}

TEST_CASE("mutual_inductance_coil_pair — repeated timing evidence") {
    constexpr int repetitions = 12;
    constexpr int n_nodes = 9;
    constexpr double rai = 0.015125;
    constexpr double rae = 0.029125;
    constexpr double la = 0.063;
    constexpr int na = 63;
    constexpr double rbi = 0.012;
    constexpr double rbe = 0.015;
    constexpr double lb = 0.07;
    constexpr int nb = 1;
    constexpr double separation = 0.03;

    volatile double separate_sink = 0.0;
    const auto separate_start = std::chrono::steady_clock::now();
    for (int i = 0; i < repetitions; ++i) {
        separate_sink += mutual_inductance_coil(
            rai, rae, la, na, rbi, rbe, lb, nb, separation, n_nodes, false);
        separate_sink += mutual_inductance_gradient_coil(
            rai, rae, la, na, rbi, rbe, lb, nb, separation, n_nodes, false);
    }
    const auto separate_end = std::chrono::steady_clock::now();

    volatile double fused_sink = 0.0;
    const auto fused_start = std::chrono::steady_clock::now();
    for (int i = 0; i < repetitions; ++i) {
        const auto pair = mutual_inductance_coil_pair(
            rai, rae, la, na, rbi, rbe, lb, nb, separation, n_nodes, false);
        fused_sink += pair.mutual + pair.gradient;
    }
    const auto fused_end = std::chrono::steady_clock::now();

    const double separate_us = std::chrono::duration<double, std::micro>(
        separate_end - separate_start).count();
    const double fused_us = std::chrono::duration<double, std::micro>(
        fused_end - fused_start).count();
    std::cout << "B1-T1 benchmark separate_us=" << separate_us
              << " fused_us=" << fused_us
              << " ratio=" << separate_us / fused_us << '\n';
    CHECK(std::isfinite(static_cast<double>(separate_sink)));
    CHECK(std::isfinite(static_cast<double>(fused_sink)));
}

// ---- coil-level (4D integration) ----

TEST_CASE("mutual_inductance_coil — single-turn degenerates to filament") {
    // A single-turn coil with negligible radial thickness should match
    // filament-level M closely.
    const double ra = 0.0135;  // mean radius
    const double dr = 0.0001;  // very thin radial width
    const double dl = 0.0001;  // very thin axial slice
    const double d  = 0.005;

    double m_fil = mutual_inductance_filament(ra, ra, d);
    double m_coil = mutual_inductance_coil(
        ra - dr, ra + dr, dl, 1,
        ra - dr, ra + dr, dl, 1, d);

    // Should agree to within a few percent for thin coils
    double rel = std::abs(m_fil - m_coil) / std::max(std::abs(m_fil), 1e-20);
    CHECK(rel < 0.05);
}

TEST_CASE("mutual_inductance_coil — driving coil to filament benchmark") {
    // TestA params: DC (ri=0.015125, re=0.029125, l=0.063, N=63)
    //               to single armature filament (ri=0.012, re=0.015, l=0.07)
    // Reference from scipy nquad via NumericalModel Eq.4.13.

    double M = mutual_inductance_coil(
        0.015125, 0.029125, 0.063, 63,
        0.012,    0.015,    0.07,   1,
        0.03);

    // Reference: 3.46866e-07 H (±5% for GL9 vs adaptive)
    CHECK(M == doctest::Approx(3.469e-07).epsilon(0.05));
}

TEST_CASE("mutual_inductance_gradient_coil — benchmark") {
    double dM = mutual_inductance_gradient_coil(
        0.015125, 0.029125, 0.063, 63,
        0.012,    0.015,    0.07,   1,
        0.03);

    // Reference: -2.16676e-05 H/m
    CHECK(dM == doctest::Approx(-2.167e-05).epsilon(0.05));
}
