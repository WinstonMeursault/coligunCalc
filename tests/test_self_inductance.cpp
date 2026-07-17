#include <doctest/doctest.h>
#include "coilgun/physics/self_inductance.hpp"

#include <cmath>

using coilgun::physics::self_inductance;
using coilgun::physics::self_inductance_exact;

// Reference values from scipy-based Python prototype (limit=500).
// Tolerance: 1e-10 relative for exact.

TEST_CASE("self_inductance — table path (in range)") {
    // TestB: q=1.25, p=1.5625 — well within table bounds
    double nc = 40.0 / ((0.05 - 0.032) * 0.04);
    double L = self_inductance(0.032, 0.05, 0.04, nc);
    CHECK(L == doctest::Approx(1.079449316808e-04).epsilon(5e-4));
}

TEST_CASE("self_inductance — exact fallback (q below range)") {
    // q = 0.07 / 0.01 = 7 / 20 = 0.35, p = 1.5 — within table
    // Let's force a fallback by using q=0.02 (below 0.05)
    double nc = 1.0 / ((0.03 - 0.01) * 0.0002);
    double L = self_inductance(0.01, 0.03, 0.0002, nc);
    double L_exact = self_inductance_exact(0.01, 0.03, 0.0002, nc);
    CHECK(L == doctest::Approx(L_exact).epsilon(1e-10));
}

TEST_CASE("self_inductance_exact — driving coil TestA") {
    // ri=0.015125, re=0.029125, l=0.063, N=63
    double nc = 63.0 / ((0.029125 - 0.015125) * 0.063);
    double L = self_inductance_exact(0.015125, 0.029125, 0.063, nc);
    CHECK(L == doctest::Approx(7.191793975090e-05).epsilon(1e-10));
}

TEST_CASE("self_inductance_exact — driving coil TestB") {
    double nc = 40.0 / ((0.05 - 0.032) * 0.04);
    double L = self_inductance_exact(0.032, 0.05, 0.04, nc);
    CHECK(L == doctest::Approx(1.079449316808e-04).epsilon(1e-10));
}

TEST_CASE("self_inductance_exact — driving coil TestC") {
    double nc = 66.0 / ((0.05 - 0.037) * 0.052);
    double L = self_inductance_exact(0.037, 0.05, 0.052, nc);
    CHECK(L == doctest::Approx(3.032879963879e-04).epsilon(1e-10));
}

TEST_CASE("self_inductance_exact — armature filament") {
    // Single current filament: ri=0.012, re=0.015, l=0.07
    double nc = 1.0 / (0.07 * (0.015 - 0.012));
    double L = self_inductance_exact(0.012, 0.015, 0.07, nc);
    CHECK(L == doctest::Approx(8.074681826873e-09).epsilon(1e-10));
}
