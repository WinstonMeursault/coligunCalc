#include <doctest/doctest.h>
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"

#include <cmath>

using coilgun::components::DrivingCoil;
using coilgun::components::Armature;

// Integration tests matching the Python prototype TestA/B/C benchmarks.
// Each test constructs a DrivingCoil + Armature pair and verifies that
// all computed properties match the scipy-based Python reference values.

TEST_CASE("Integration TestA — driving coil properties") {
    // Python prototype TestA: drivingCoil(0.015125, 0.029125, 0.063, 63,
    //                                     1.75e-8, 2.8e-5, 0.7)
    DrivingCoil dc(0.015125, 0.029125, 0.063, 63,
                   1.75e-8, 2.8e-5, 0.7);

    CHECK(dc.inner_radius()  == doctest::Approx(0.015125).epsilon(1e-12));
    CHECK(dc.outer_radius()  == doctest::Approx(0.029125).epsilon(1e-12));
    CHECK(dc.length()        == doctest::Approx(0.063).epsilon(1e-12));
    CHECK(dc.turns()         == 63);
    CHECK(dc.mean_radius()   == doctest::Approx(0.022125).epsilon(1e-12));

    // nc = 63 / ((0.029125 - 0.015125) * 0.063) = 71428.57 turns/m^2
    CHECK(dc.turns_density() == doctest::Approx(71428.57).epsilon(1e-3));

    // From Python scipy: L = 7.19179e-05 H
    CHECK(dc.self_inductance() == doctest::Approx(7.19179e-05).epsilon(1e-10));

    // Default position: l/2
    CHECK(dc.position() == doctest::Approx(0.0315).epsilon(1e-12));
}

TEST_CASE("Integration TestA — armature properties (m=1, n=1)") {
    // Python prototype TestA: armature(0.012, 0.015, 0.07, 1.75e-8,
    //                                  0, 6.384, 1, 1, 0.0505)
    // density of aluminum = 2700 kg/m^3
    Armature arm(0.012, 0.015, 0.07, 1.75e-8, 2700.0,
                 0.0, 6.384, 1, 1, 0.0505);

    CHECK(arm.inner_radius()   == doctest::Approx(0.012).epsilon(1e-12));
    CHECK(arm.outer_radius()   == doctest::Approx(0.015).epsilon(1e-12));
    CHECK(arm.length()         == doctest::Approx(0.07).epsilon(1e-12));
    CHECK(arm.axial_filaments()  == 1);
    CHECK(arm.radial_filaments() == 1);
    CHECK(arm.total_filaments()  == 1);
    CHECK(arm.position() == doctest::Approx(0.0505).epsilon(1e-12));
    CHECK(arm.velocity() == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(arm.mass()     == doctest::Approx(6.384).epsilon(1e-12));

    // Filament geometry (1-indexed)
    CHECK(arm.filament_inner_radius(1)  == doctest::Approx(0.012).epsilon(1e-12));
    CHECK(arm.filament_outer_radius(1)  == doctest::Approx(0.015).epsilon(1e-12));
    CHECK(arm.filament_mean_radius(1)   == doctest::Approx(0.0135).epsilon(1e-12));

    // Axial position: x - l/2 + (i-0.5)*dl = 0.0505 - 0.035 + 0.5*0.07 = 0.0505
    CHECK(arm.filament_axial_position(1) == doctest::Approx(0.0505).epsilon(1e-12));

    // Self-inductance: from Python scipy = 8.07468e-09 H
    CHECK(arm.inductances()[0] == doctest::Approx(8.07468e-09).epsilon(1e-9));

    // Update position
    arm.update_position(0.001);
    CHECK(arm.position() == doctest::Approx(0.0515).epsilon(1e-12));
    CHECK(arm.filament_axial_position(1) == doctest::Approx(0.0515).epsilon(1e-12));
}

TEST_CASE("Integration TestB") {
    // TestB: DC (0.032, 0.05, 0.04, 40, 1.75e-8, 2.8e-5, 0.7)
    DrivingCoil dc(0.032, 0.05, 0.04, 40, 1.75e-8, 2.8e-5, 0.7);

    CHECK(dc.turns_density() == doctest::Approx(55555.56).epsilon(1e-3));
    CHECK(dc.self_inductance() == doctest::Approx(1.07945e-04).epsilon(5e-4));
}

TEST_CASE("Integration TestC") {
    // TestC: DC (0.037, 0.05, 0.052, 66, 1.75e-8, 2.8e-5, 0.7)
    DrivingCoil dc(0.037, 0.05, 0.052, 66, 1.75e-8, 2.8e-5, 0.7);

    CHECK(dc.turns_density() == doctest::Approx(97633.14).epsilon(1e-3));
    CHECK(dc.self_inductance() == doctest::Approx(3.03288e-04).epsilon(5e-4));
}
