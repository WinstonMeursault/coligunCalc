/**
 * @file test_gpu_batch.cpp
 * @brief Integration test: 10-point parameter sweep — GPU vs CPU.
 * @author Winston Meursault
 */

#include <doctest/doctest.h>
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/physics/constants.hpp"
#include <memory>
#include <vector>

using namespace coilgun::physics;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::simulation::MultiStageSim;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::TriggerConfig;
using coilgun::simulation::TriggerMode;
using coilgun::simulation::cuda::GpuMultiStageSim;
using coilgun::simulation::cuda::GpuBackend;
using coilgun::simulation::cuda::GpuOptLevel;

struct SweepPoint { double voltage; };

TEST_CASE("GPU batch — 10 parameter sweep vs CPU") {
    std::vector<SweepPoint> sweeps = {
        {300.0}, {320.0}, {340.0}, {360.0}, {380.0},
        {400.0}, {420.0}, {440.0}, {460.0}, {480.0}
    };
    REQUIRE(sweeps.size() == 10);

    auto make_coils = [](double pos_c2) -> std::vector<DrivingCoil> {
        DrivingCoil c1(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
        DrivingCoil c2(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, pos_c2);
        return {c1, c2};
    };

    auto make_arm = []() {
        return Armature(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref,
                        ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    };

    std::vector<double> v_cpu(10), v_gpu(10);

    for (int i = 0; i < 10; ++i) {
        auto coils = make_coils(0.10);
        auto arm   = make_arm();
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
        MultiStageSim<EulerStepper> sim(coils, arm, std::move(excs), triggers, 1e-6);
        sim.run();
        v_cpu[i] = sim.result().summary.muzzle_velocity;
    }

    for (int i = 0; i < 10; ++i) {
        auto coils = make_coils(0.10);
        auto arm   = make_arm();
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
        GpuBackend be;
        be.use_persistent = false;
        GpuMultiStageSim<EulerStepper> sim(coils, arm, std::move(excs), triggers, 1e-6,
                                            false, GpuOptLevel::Full, be);
        sim.run();
        v_gpu[i] = sim.result().summary.muzzle_velocity;
    }

    for (int i = 0; i < 10; ++i) {
        INFO("Sweep " << i << " voltage=" << sweeps[i].voltage);
        CHECK(v_gpu[i] == doctest::Approx(v_cpu[i]).epsilon(1e-5));
    }
}
