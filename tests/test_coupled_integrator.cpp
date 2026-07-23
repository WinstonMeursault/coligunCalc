#include <doctest/doctest.h>

#include "coilgun/coilgun.hpp"
#include "test_fixtures.hpp"

#include <memory>

using namespace coilgun;

TEST_CASE("single-stage Euler uses pre-step current and records post-step time") {
    auto fixture = test::make_test_coil_and_armature();
    auto excitation = std::make_unique<simulation::CapacitorExcitation>(100.0, 1.0);
    simulation::SingleStageSim<simulation::EulerStepper> sim(
        fixture.coil, fixture.armature, std::move(excitation), 1e-4, true);
    const auto& first = sim.step();
    CHECK(first.time == doctest::Approx(1e-4));
    CHECK(first.cap_voltage == doctest::Approx(100.0));
    for (const double temperature : first.filament_temperatures)
        CHECK(temperature == doctest::Approx(physics::T_REFERENCE));
}

TEST_CASE("single-stage RK4 advances complete coupled state") {
    auto fixture = test::make_test_coil_and_armature();
    auto excitation = std::make_unique<simulation::WaveformExcitation>(
        [](double) { return 10.0; });
    simulation::SingleStageSim<simulation::RK4Stepper> sim(
        fixture.coil, fixture.armature, std::move(excitation), 1e-6, true);
    const auto& step = sim.step();
    CHECK(step.time == doctest::Approx(1e-6));
    CHECK(std::isfinite(step.coil_current));
    CHECK(std::isfinite(step.arm_velocity));
}

TEST_CASE("excitation completion does not immediately remove residual current") {
    auto sim = test::make_waveform_stage_with_nonzero_current();
    const auto& step = sim.step();
    CHECK(sim.stage_state(0).excitation_finished);
    CHECK(sim.stage_state(0).circuit_active ==
          (std::abs(step.coil_currents[0]) >= 1e-6));
}
