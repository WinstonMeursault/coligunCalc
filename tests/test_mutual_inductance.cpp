#include <doctest/doctest.h>
#include "coilgun/physics/mutual_inductance.hpp"

#include <algorithm>
#include <cmath>

using coilgun::physics::mutual_inductance_filament;
using coilgun::physics::mutual_inductance_gradient_filament;
using coilgun::physics::mutual_inductance_coil;
using coilgun::physics::mutual_inductance_gradient_coil;

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
