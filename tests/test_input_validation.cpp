#include <doctest/doctest.h>

#include "coilgun/coilgun.hpp"

#include <limits>
#include <stdexcept>

using namespace coilgun;

TEST_CASE("DrivingCoil validates geometry and preserves explicit negative position") {
    CHECK_THROWS_AS(components::DrivingCoil(
        0.03, 0.01, 0.05, 100, physics::COPPER.resistivity_ref, 1e-6, 0.7),
        std::invalid_argument);
    CHECK_THROWS_AS(components::DrivingCoil(
        0.01, 0.03, 0.0, 100, physics::COPPER.resistivity_ref, 1e-6, 0.7),
        std::invalid_argument);
    components::DrivingCoil coil(
        0.01, 0.03, 0.05, 100, physics::COPPER.resistivity_ref, 1e-6, 0.7,
        -0.25);
    CHECK(coil.position() == doctest::Approx(-0.25));
}

TEST_CASE("Armature validates dimensions and mass before division") {
    CHECK_THROWS_AS(components::Armature(
        0.03, 0.01, 0.08, physics::ALUMINUM.resistivity_ref,
        physics::ALUMINUM.density, 0.0, 0.1, 5, 2, 0.0),
        std::invalid_argument);
    CHECK_THROWS_AS(components::Armature(
        0.005, 0.025, 0.08, physics::ALUMINUM.resistivity_ref,
        physics::ALUMINUM.density, 0.0, 0.0, 5, 2, 0.0),
        std::invalid_argument);
    CHECK_THROWS_AS(components::Armature(
        0.005, 0.025, 0.08, physics::ALUMINUM.resistivity_ref,
        physics::ALUMINUM.density, 0.0, 0.1, 0, 2, 0.0),
        std::invalid_argument);
}

TEST_CASE("excitation and trigger validation rejects malformed input") {
    CHECK_THROWS_AS(simulation::CapacitorExcitation(100.0, 0.0),
                    std::invalid_argument);
    CHECK_THROWS_AS(simulation::CapacitorExcitation(0.0, 1.0),
                    std::invalid_argument);
    CHECK_THROWS_AS(simulation::WaveformExcitation({}), std::invalid_argument);
    simulation::WaveformExcitation waveform([](double) { return 1.0; });
    CHECK_THROWS_AS(waveform.set_end_time(-1.0), std::invalid_argument);
    CHECK_NOTHROW(waveform.set_end_time(std::numeric_limits<double>::infinity()));
}
