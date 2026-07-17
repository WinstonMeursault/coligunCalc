#include <doctest/doctest.h>
#include "coilgun/components/armature.hpp"

#include <cmath>

using coilgun::components::Armature;

namespace {

// TestA: rai=0.012, rae=0.015, la=0.07, rho=1.75e-8, density=2700,
//        v0=0, ma=6.384, m=1, n=1, x0=0.0505
constexpr double k_ri      = 0.012;
constexpr double k_re      = 0.015;
constexpr double k_l       = 0.07;
constexpr double k_rho     = 1.75e-8;
constexpr double k_density = 2700.0;
constexpr double k_v0      = 0.0;
constexpr double k_ma      = 6.384;
constexpr int    k_m       = 1;
constexpr int    k_n       = 1;
constexpr double k_x0      = 0.0505;

} // namespace

TEST_CASE("Armature — filament geometry (m=1, n=1)") {
    Armature arm(k_ri, k_re, k_l, k_rho, k_density, k_v0, k_ma, k_m, k_n, k_x0);

    CHECK(arm.inner_radius() == k_ri);
    CHECK(arm.outer_radius() == k_re);
    CHECK(arm.length() == k_l);
    CHECK(arm.axial_filaments() == k_m);
    CHECK(arm.radial_filaments() == k_n);
    CHECK(arm.total_filaments() == k_m * k_n);

    CHECK(arm.filament_inner_radius(1) == doctest::Approx(k_ri));
    CHECK(arm.filament_outer_radius(1) == doctest::Approx(k_re));
    // mean radius: ri + (re-ri)*(1-0.5)/n = ri + 0.5*dr = 0.012 + 0.0015 = 0.0135
    CHECK(arm.filament_mean_radius(1) == doctest::Approx(0.0135));
    // axial position(1): x - l/2 + (1-0.5)*dl = 0.0505 - 0.035 + 0.5*0.07 = 0.0505
    CHECK(arm.filament_axial_position(1) == doctest::Approx(0.0505));
}

TEST_CASE("Armature — resistance (TestA, m=1, n=1)") {
    Armature arm(k_ri, k_re, k_l, k_rho, k_density, k_v0, k_ma, k_m, k_n, k_x0);

    // Eq.6.7: R = 2*pi*rho * [m/(2*l) + m*n*ri/(l*(re-ri)) + (j-1)*m/l]
    double expected = 2.0 * M_PI * k_rho
        * (k_m / (2.0 * k_l) + static_cast<double>(k_m * k_n) * k_ri / (k_l * (k_re - k_ri))
           + (1 - 1) * k_m / k_l);

    CHECK(arm.resistances().size() == k_m * k_n);
    CHECK(arm.resistances()[0] == doctest::Approx(expected));
}

TEST_CASE("Armature — inductance (TestA, m=1, n=1)") {
    Armature arm(k_ri, k_re, k_l, k_rho, k_density, k_v0, k_ma, k_m, k_n, k_x0);

    // From test_self_inductance: nc = 1/(0.07*(0.015-0.012)), L ≈ 8.07468e-09 H
    CHECK(arm.inductances().size() == k_m * k_n);
    CHECK(arm.inductances()[0] == doctest::Approx(8.074681826873e-09).epsilon(1e-8));
}

TEST_CASE("Armature — filament mass (TestA, m=1, n=1)") {
    Armature arm(k_ri, k_re, k_l, k_rho, k_density, k_v0, k_ma, k_m, k_n, k_x0);

    // mass = density * pi * (re² - ri²) * dl
    double expected = k_density * M_PI * (k_re * k_re - k_ri * k_ri) * (k_l / k_m);
    CHECK(arm.masses().size() == k_m * k_n);
    CHECK(arm.masses()[0] == doctest::Approx(expected));
}

TEST_CASE("Armature — motion state") {
    Armature arm(k_ri, k_re, k_l, k_rho, k_density, k_v0, k_ma, k_m, k_n, k_x0);

    CHECK(arm.position() == doctest::Approx(k_x0));
    CHECK(arm.velocity() == doctest::Approx(k_v0));
    CHECK(arm.mass() == doctest::Approx(k_ma));

    arm.update_position(0.01);
    CHECK(arm.position() == doctest::Approx(k_x0 + 0.01));

    arm.set_velocity(42.0);
    CHECK(arm.velocity() == doctest::Approx(42.0));
}

TEST_CASE("Armature — filament indexing (m=2, n=3)") {
    // Multi-filament test: verify row-major ordering
    Armature arm(k_ri, k_re, k_l, k_rho, k_density, k_v0, k_ma, 2, 3, k_x0);

    CHECK(arm.axial_filaments() == 2);
    CHECK(arm.radial_filaments() == 3);
    CHECK(arm.total_filaments() == 6);
    CHECK(arm.resistances().size() == 6);
    CHECK(arm.inductances().size() == 6);
    CHECK(arm.masses().size() == 6);

    // Resistance should be same for same radial layer (j)
    double dr = (k_re - k_ri) / 3.0;
    for (int j = 1; j <= 3; ++j) {
        for (int i = 1; i <= 2; ++i) {
            double r_inner = k_ri + (j - 1) * dr;
            double r_outer = k_ri + j * dr;
            CHECK(arm.filament_inner_radius(j) == doctest::Approx(r_inner));
            CHECK(arm.filament_outer_radius(j) == doctest::Approx(r_outer));
        }
    }

    // Insist that filaments in the same radial layer have the same resistance and L
    for (int j = 0; j < 3; ++j) {
        double val_r = arm.resistances()[j];
        double val_l = arm.inductances()[j];
        double val_m = arm.masses()[j];
        for (int i = 1; i < 2; ++i) {
            CHECK(arm.resistances()[i * 3 + j] == doctest::Approx(val_r));
            CHECK(arm.inductances()[i * 3 + j] == doctest::Approx(val_l));
            CHECK(arm.masses()[i * 3 + j] == doctest::Approx(val_m));
        }
    }
}

TEST_CASE("Armature — axial position updates with update_position") {
    Armature arm(k_ri, k_re, k_l, k_rho, k_density, k_v0, k_ma, k_m, k_n, k_x0);

    double x1 = arm.filament_axial_position(1);
    arm.update_position(0.005);
    double x2 = arm.filament_axial_position(1);
    CHECK(x2 == doctest::Approx(x1 + 0.005));
}
