#include <doctest/doctest.h>
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/simulation/cuda/gpu_execution_report.hpp"
#include "coilgun/simulation/cuda/gpu_mutual_pipeline.hpp"
#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/single_stage_sim.hpp"
#include "gpu_numerical_tolerances.hpp"
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using coilgun::physics::mutual_inductance_coil;
using coilgun::physics::mutual_inductance_gradient_coil;
using coilgun::simulation::cuda::CoilGeo;
using coilgun::simulation::cuda::FilGeo;
using coilgun::simulation::cuda::GpuOptLevel;
using coilgun::simulation::cuda::GpuBackend;
using coilgun::simulation::cuda::BackendMode;
using coilgun::simulation::cuda::ExecutionReport;
using coilgun::simulation::cuda::GpuSingleStageSim;
using coilgun::simulation::cuda::MutualPipelineView;
using coilgun::simulation::cuda::launch_mutual_pipeline;
using coilgun::simulation::cuda::mutual_pipeline_index;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::SingleStageSim;

namespace {

using gpu_test::numerically_equal;
using gpu_test::tolerance_for;

double relative_error(double actual, double expected) {
    const double scale = std::max(std::abs(actual), std::abs(expected));
    return scale > 0.0 ? std::abs(actual - expected) / scale : 0.0;
}

struct SingleStageRun {
    coilgun::simulation::SimResult cpu;
    coilgun::simulation::SimResult gpu;
    coilgun::simulation::SimState gpu_state;
    std::vector<double> gpu_resistances;
    std::vector<double> cpu_resistances;
    ExecutionReport report;
};

SingleStageRun run_single_stage(GpuOptLevel mode, bool thermal = false) {
    using coilgun::components::Armature;
    using coilgun::components::DrivingCoil;
    using namespace coilgun::physics;

    DrivingCoil cpu_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature cpu_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                     0.0, 0.120, 5, 2, 0.05);
    const auto reference_resistances = cpu_arm.resistances();
    auto cpu_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> cpu(std::move(cpu_coil), std::move(cpu_arm),
                                     std::move(cpu_exc), 1e-6, thermal);

    DrivingCoil gpu_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature gpu_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                     0.0, 0.120, 5, 2, 0.05);
    auto gpu_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> gpu(std::move(gpu_coil), std::move(gpu_arm),
                                        std::move(gpu_exc), 1e-6, thermal, mode,
                                        backend);

    coilgun::simulation::TerminationPolicy policy;
    policy.max_steps = 20;
    policy.enable_velocity_check = false;
    policy.enable_bound_check = false;
    cpu.run(policy);
    gpu.run(policy);

    std::vector<double> expected_resistances;
    if (thermal) {
        expected_resistances.reserve(reference_resistances.size());
        const auto& final_temperatures = cpu.result().history.back().filament_temperatures;
        for (std::size_t i = 0; i < reference_resistances.size(); ++i) {
            expected_resistances.push_back(reference_resistances[i] *
                (1.0 + ALUMINUM.temp_coefficient *
                    (final_temperatures[i] - coilgun::physics::T_REFERENCE)));
        }
    }
    return {cpu.result(), gpu.result(), gpu.state(), gpu.filament_resistances(),
            std::move(expected_resistances), gpu.execution_report()};
}

} // namespace

TEST_CASE("GPU mutual pipeline uses fixed B-S-F indexing and active masks") {
    CHECK(mutual_pipeline_index(1, 2, 3, 4, 5) == 33);
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count == 0) {
        MESSAGE("CUDA device unavailable; skipping mutual-kernel test");
        return;
    }

    const CoilGeo coils[] = {{0.01, 0.03, 0.05, 0.0, 150},
                             {0.01, 0.03, 0.05, 0.0, 150}};
    const FilGeo filaments[] = {{0.005, 0.01, 0.016}};
    const double separations[] = {0.03, 0.04, 0.03, 0.04};
    const std::uint8_t active[] = {1, 0, 1, 1};
    CoilGeo* d_coils = nullptr;
    FilGeo* d_filaments = nullptr;
    double* d_separations = nullptr;
    std::uint8_t* d_active = nullptr;
    double* d_mutual = nullptr;
    double* d_gradient = nullptr;
    cudaMalloc(&d_coils, sizeof(coils));
    cudaMalloc(&d_filaments, sizeof(filaments));
    cudaMalloc(&d_separations, sizeof(separations));
    cudaMalloc(&d_active, sizeof(active));
    cudaMalloc(&d_mutual, sizeof(separations));
    cudaMalloc(&d_gradient, sizeof(separations));
    cudaMemcpy(d_coils, coils, sizeof(coils), cudaMemcpyHostToDevice);
    cudaMemcpy(d_filaments, filaments, sizeof(filaments), cudaMemcpyHostToDevice);
    cudaMemcpy(d_separations, separations, sizeof(separations), cudaMemcpyHostToDevice);
    cudaMemcpy(d_active, active, sizeof(active), cudaMemcpyHostToDevice);

    MutualPipelineView view{d_coils, d_filaments, d_separations, d_active,
                            d_mutual, d_gradient, 2, 2, 1, 9};
    launch_mutual_pipeline(view, GpuOptLevel::Standard);
    REQUIRE(cudaDeviceSynchronize() == cudaSuccess);

    double mutual[4] = {};
    double gradient[4] = {};
    cudaMemcpy(mutual, d_mutual, sizeof(mutual), cudaMemcpyDeviceToHost);
    cudaMemcpy(gradient, d_gradient, sizeof(gradient), cudaMemcpyDeviceToHost);
    CHECK(mutual[1] == 0.0);
    CHECK(gradient[1] == 0.0);
    CHECK(mutual[0] == doctest::Approx(mutual_inductance_coil(
        0.01, 0.03, 0.05, 150, 0.005, 0.01, 0.016, 1, 0.03, 9, false)).epsilon(5e-7));
    CHECK(gradient[0] == doctest::Approx(mutual_inductance_gradient_coil(
        0.01, 0.03, 0.05, 150, 0.005, 0.01, 0.016, 1, 0.03, 9, false)).epsilon(5e-7));

    cudaFree(d_coils);
    cudaFree(d_filaments);
    cudaFree(d_separations);
    cudaFree(d_active);
    cudaFree(d_mutual);
    cudaFree(d_gradient);
}

TEST_CASE("GPU mutual pipeline exposes all precision modes") {
    CHECK(GpuOptLevel::Standard != GpuOptLevel::Full);
    CHECK(GpuOptLevel::Full != GpuOptLevel::Aggressive);
}

TEST_CASE("GPU mutual pipeline matches CPU M and dM for every precision mode" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    const CoilGeo coils[] = {{0.01, 0.03, 0.05, 0.0, 150}};
    const FilGeo filaments[] = {{0.005, 0.01, 0.016}};
    const double separations[] = {0.03};
    const std::uint8_t active[] = {1};
    const double expected_mutual = mutual_inductance_coil(
        0.01, 0.03, 0.05, 150, 0.005, 0.01, 0.016, 1, 0.03, 9, false);
    const double expected_gradient = mutual_inductance_gradient_coil(
        0.01, 0.03, 0.05, 150, 0.005, 0.01, 0.016, 1, 0.03, 9, false);

    for (const auto mode : {GpuOptLevel::Standard, GpuOptLevel::Full,
                            GpuOptLevel::Aggressive}) {
        double* d_mutual = nullptr;
        double* d_gradient = nullptr;
        CoilGeo* d_coils = nullptr;
        FilGeo* d_filaments = nullptr;
        double* d_separations = nullptr;
        std::uint8_t* d_active = nullptr;
        REQUIRE(cudaMalloc(&d_coils, sizeof(coils)) == cudaSuccess);
        REQUIRE(cudaMalloc(&d_filaments, sizeof(filaments)) == cudaSuccess);
        REQUIRE(cudaMalloc(&d_separations, sizeof(separations)) == cudaSuccess);
        REQUIRE(cudaMalloc(&d_active, sizeof(active)) == cudaSuccess);
        REQUIRE(cudaMalloc(&d_mutual, sizeof(double)) == cudaSuccess);
        REQUIRE(cudaMalloc(&d_gradient, sizeof(double)) == cudaSuccess);
        cudaMemcpy(d_coils, coils, sizeof(coils), cudaMemcpyHostToDevice);
        cudaMemcpy(d_filaments, filaments, sizeof(filaments), cudaMemcpyHostToDevice);
        cudaMemcpy(d_separations, separations, sizeof(separations), cudaMemcpyHostToDevice);
        cudaMemcpy(d_active, active, sizeof(active), cudaMemcpyHostToDevice);

        MutualPipelineView view{d_coils, d_filaments, d_separations, d_active,
                                d_mutual, d_gradient, 1, 1, 1, 9};
        launch_mutual_pipeline(view, mode);
        REQUIRE(cudaDeviceSynchronize() == cudaSuccess);
        double actual_mutual = 0.0;
        double actual_gradient = 0.0;
        REQUIRE(cudaMemcpy(&actual_mutual, d_mutual, sizeof(double),
                           cudaMemcpyDeviceToHost) == cudaSuccess);
        REQUIRE(cudaMemcpy(&actual_gradient, d_gradient, sizeof(double),
                           cudaMemcpyDeviceToHost) == cudaSuccess);

        const auto tolerance = tolerance_for(mode);
        MESSAGE("precision_id=" << static_cast<int>(mode)
             << " M_relative_error=" << relative_error(actual_mutual, expected_mutual)
             << " dM_relative_error=" << relative_error(actual_gradient, expected_gradient));
        CHECK(numerically_equal(actual_mutual, expected_mutual, tolerance));
        CHECK(numerically_equal(actual_gradient, expected_gradient, tolerance));

        cudaFree(d_coils);
        cudaFree(d_filaments);
        cudaFree(d_separations);
        cudaFree(d_active);
        cudaFree(d_mutual);
        cudaFree(d_gradient);
    }
}

TEST_CASE("GPU precision modes match CPU stepwise and final single-stage reference" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    for (const auto mode : {GpuOptLevel::Standard, GpuOptLevel::Full,
                            GpuOptLevel::Aggressive}) {
        const auto run = run_single_stage(mode);
        REQUIRE(run.report.gpu_executed);
        CHECK(run.report.backend == BackendMode::Direct);
        CHECK(run.report.requested_precision ==
              static_cast<coilgun::simulation::cuda::PrecisionMode>(mode));
        CHECK(run.report.precision ==
              static_cast<coilgun::simulation::cuda::PrecisionMode>(mode));
        REQUIRE(run.gpu.history.size() == run.cpu.history.size());

        const auto tolerance = tolerance_for(mode);
        INFO("precision mode=" << coilgun::simulation::cuda::to_string(
            static_cast<coilgun::simulation::cuda::PrecisionMode>(mode)));
        double max_coil_error = 0.0;
        double max_position_error = 0.0;
        double max_velocity_error = 0.0;
        double max_force_error = 0.0;
        double max_filament_error = 0.0;
        for (std::size_t i = 0; i < run.cpu.history.size(); ++i) {
            const auto& expected = run.cpu.history[i];
            const auto& actual = run.gpu.history[i];
            max_coil_error = std::max(max_coil_error,
                relative_error(actual.coil_current, expected.coil_current));
            max_position_error = std::max(max_position_error,
                relative_error(actual.arm_position, expected.arm_position));
            max_velocity_error = std::max(max_velocity_error,
                relative_error(actual.arm_velocity, expected.arm_velocity));
            max_force_error = std::max(max_force_error,
                relative_error(actual.force, expected.force));
            REQUIRE(actual.filament_currents.size() == expected.filament_currents.size());
            for (std::size_t k = 0; k < expected.filament_currents.size(); ++k)
                max_filament_error = std::max(max_filament_error,
                    relative_error(actual.filament_currents[k], expected.filament_currents[k]));
        }
        MESSAGE("precision_id=" << static_cast<int>(mode)
             << " max_relative_error coil=" << max_coil_error
             << " position=" << max_position_error
             << " velocity=" << max_velocity_error
             << " force=" << max_force_error
             << " filament=" << max_filament_error);
        CHECK(max_coil_error <= 1.0);
        CHECK(max_position_error <= 1.0);
        CHECK(max_velocity_error <= 1.0);
        CHECK(max_force_error <= 1.0);
        CHECK(max_filament_error <= 1.0);
        CHECK(numerically_equal(run.gpu.summary.muzzle_velocity,
                                run.cpu.summary.muzzle_velocity, tolerance));
        CHECK(numerically_equal(run.gpu.summary.max_force,
                                run.cpu.summary.max_force, tolerance));
    }
}

TEST_CASE("GPU precision modes preserve thermal temperature and resistance thresholds" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    for (const auto mode : {GpuOptLevel::Standard, GpuOptLevel::Full,
                            GpuOptLevel::Aggressive}) {
        const auto run = run_single_stage(mode, true);
        REQUIRE(run.report.gpu_executed);
        REQUIRE(run.gpu.history.back().filament_temperatures.size() ==
                run.cpu.history.back().filament_temperatures.size());
        const auto tolerance = tolerance_for(mode);
        for (std::size_t i = 0; i < run.gpu.history.back().filament_temperatures.size(); ++i) {
            CHECK(numerically_equal(run.gpu.history.back().filament_temperatures[i],
                                    run.cpu.history.back().filament_temperatures[i], tolerance));
        }
        CHECK(run.gpu_resistances.size() == run.gpu_state.currents.size() - 1);
        REQUIRE(run.gpu_resistances.size() == run.cpu_resistances.size());
        for (std::size_t i = 0; i < run.gpu_resistances.size(); ++i) {
            CHECK(std::isfinite(run.gpu_resistances[i]));
            CHECK(numerically_equal(run.gpu_resistances[i], run.cpu_resistances[i], tolerance));
        }
    }
}
