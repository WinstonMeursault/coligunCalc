#include <doctest/doctest.h>
#include "coilgun/simulation/single_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/physics/constants.hpp"
#include <memory>

using namespace coilgun::physics;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::SingleStageSim;
using coilgun::simulation::cuda::GpuSingleStageSim;

static double run_cpu_single() {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false);
    sim.run();
    return sim.result().summary.muzzle_velocity;
}

static double run_gpu_single() {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    coilgun::simulation::cuda::GpuBackend be;
    be.use_persistent = false;  // fallback for deterministic GPU vs CPU comparison
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false,
                                        coilgun::simulation::cuda::GpuOptLevel::Full, be);
    sim.run();
    return sim.result().summary.muzzle_velocity;
}

TEST_CASE("GPU single-stage muzzle velocity matches CPU") {
    double v_cpu = run_cpu_single();
    double v_gpu = run_gpu_single();
    CHECK(v_gpu == doctest::Approx(v_cpu).epsilon(5e-3));
}

TEST_CASE("GPU single-stage muzzle velocity direction consistent") {
    double v_gpu = run_gpu_single();
    CHECK(v_gpu > 0.0);
}
