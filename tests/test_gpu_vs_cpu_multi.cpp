/**
 * @file test_gpu_vs_cpu_multi.cpp
 * @brief End-to-end test: GPU vs CPU multi-stage simulation.
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

static double run_cpu_multi() {
    DrivingCoil c1(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);
    std::vector<DrivingCoil> coils = {c1, c2};
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
    MultiStageSim<EulerStepper> sim(std::move(coils), arm, std::move(excs), triggers, 1e-6);
    sim.run();
    return sim.result().summary.muzzle_velocity;
}

static double run_gpu_multi() {
    DrivingCoil c1(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);
    std::vector<DrivingCoil> coils = {c1, c2};
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
    GpuBackend be;
    be.use_persistent = false;  // fallback for deterministic GPU vs CPU comparison
    GpuMultiStageSim<EulerStepper> sim(std::move(coils), arm, std::move(excs), triggers,
                                       1e-6, false, GpuOptLevel::Full, be);
    sim.run();
    return sim.result().summary.muzzle_velocity;
}

TEST_CASE("GPU multi-stage muzzle velocity matches CPU") {
    double v_cpu = run_cpu_multi();
    double v_gpu = run_gpu_multi();
    CHECK(v_gpu == doctest::Approx(v_cpu).epsilon(5e-3));
}
