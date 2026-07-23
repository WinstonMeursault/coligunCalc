#pragma once

#include "coilgun/coilgun.hpp"

#include <memory>
#include <vector>

namespace coilgun::test {

struct SingleStageFixture {
    components::DrivingCoil coil;
    components::Armature armature;
};

inline components::DrivingCoil make_test_coil() {
    return {0.005, 0.010, 0.010, 10, physics::COPPER.resistivity_ref,
            1e-6, 0.7};
}

inline components::Armature make_thermal_armature() {
    return {0.002, 0.008, 0.010, physics::ALUMINUM.resistivity_ref,
            physics::ALUMINUM.density, 0.0, 0.005, 1, 1, 0.008};
}

inline SingleStageFixture make_test_coil_and_armature() {
    return {make_test_coil(), make_thermal_armature()};
}

inline simulation::MultiStageSim<simulation::EulerStepper>
make_waveform_stage_with_nonzero_current() {
    auto waveform = std::make_unique<simulation::WaveformExcitation>(
        [](double) { return 100.0; });
    waveform->set_end_time(1e-6);
    std::vector<std::unique_ptr<simulation::Excitation>> excitations;
    excitations.push_back(std::move(waveform));
    return {{make_test_coil()}, make_thermal_armature(), std::move(excitations),
            {}, 1e-6};
}

inline simulation::MultiStageSim<simulation::EulerStepper>
make_stage_that_decays_below_threshold() {
    auto waveform = std::make_unique<simulation::WaveformExcitation>(
        [](double) { return 1.0; });
    waveform->set_end_time(1e-6);
    std::vector<std::unique_ptr<simulation::Excitation>> excitations;
    excitations.push_back(std::move(waveform));
    return {{make_test_coil()}, make_thermal_armature(), std::move(excitations),
            {}, 1e-6};
}

inline simulation::MultiStageSim<simulation::EulerStepper>
make_two_stage_position_trigger(double position) {
    auto first = make_test_coil();
    auto second = make_test_coil();
    first.set_position(0.0);
    second.set_position(0.03);
    std::vector<std::unique_ptr<simulation::Excitation>> excitations;
    excitations.push_back(
        std::make_unique<simulation::CapacitorExcitation>(100.0, 1e-3));
    excitations.push_back(
        std::make_unique<simulation::CapacitorExcitation>(100.0, 1e-3));
    return {{first, second}, make_thermal_armature(), std::move(excitations),
            {{simulation::TriggerMode::Position, position}}, 1e-6};
}

} // namespace coilgun::test
