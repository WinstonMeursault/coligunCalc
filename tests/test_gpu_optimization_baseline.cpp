/**
 * @file test_gpu_optimization_baseline.cpp
 * @brief Real-device CPU/CUDA numerical and batch-timing baseline.
 */

#include <doctest/doctest.h>

#include "coilgun/components/armature.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
#include "coilgun/simulation/cuda/sim_batch.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/termination.hpp"
#include "gpu_numerical_tolerances.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

using coilgun::components::Armature;
using coilgun::components::DrivingCoil;
using coilgun::physics::ALUMINUM;
using coilgun::physics::COPPER;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::MultiStageResult;
using coilgun::simulation::MultiStageSim;
using coilgun::simulation::OptimizationLevel;
using coilgun::simulation::TriggerConfig;
using coilgun::simulation::TriggerMode;
using coilgun::simulation::TerminationPolicy;
using coilgun::simulation::cuda::BackendMode;
using coilgun::simulation::cuda::GpuBackend;
using coilgun::simulation::cuda::GpuMultiStageSim;
using coilgun::simulation::cuda::GpuOptLevel;
using coilgun::simulation::cuda::SimBatch;

struct Candidate {
    int id;
    double stage0_voltage;
    double stage1_voltage;
};

struct Metrics {
    double muzzle_velocity = 0.0;
    double peak_coil_current = 0.0;
    double peak_cap_voltage = 0.0;
    double max_filament_temperature = 0.0;
};

struct TimedBatch {
    double host_wall_ms = 0.0;
    std::vector<Metrics> metrics;
    coilgun::simulation::cuda::ExecutionReport report;
};

constexpr double kCudaPeakCurrentBaselineMinRelativeShortfall = 0.05;
constexpr double kCudaPeakCurrentBaselineMaxRelativeShortfall = 0.08;

std::vector<DrivingCoil> make_coils() {
    return {
        DrivingCoil(0.010, 0.030, 0.050, 150, COPPER.resistivity_ref,
                    1e-6, 0.7, 0.015),
        DrivingCoil(0.010, 0.030, 0.050, 150, COPPER.resistivity_ref,
                    1e-6, 0.7, 0.085),
    };
}

Armature make_armature() {
    return Armature(0.005, 0.025, 0.080, ALUMINUM.resistivity_ref,
                    ALUMINUM.density, 0.0, 0.120, 2, 2, 0.000);
}

std::vector<Candidate> make_candidates() {
    std::vector<Candidate> candidates;
    candidates.reserve(128);
    for (int id = 0; id < 128; ++id) {
        // Distinct values make accidental row permutation observable while
        // keeping every row on the same fixed geometry and trigger schedule.
        candidates.push_back({id, 280.0 + 1.25 * id, 215.0 + 0.875 * id});
    }
    return candidates;
}

std::vector<std::unique_ptr<coilgun::simulation::Excitation>>
make_excitations(const Candidate& candidate) {
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excitations;
    excitations.push_back(std::make_unique<CrowbarExcitation>(
        candidate.stage0_voltage, 0.0008));
    excitations.push_back(std::make_unique<CrowbarExcitation>(
        candidate.stage1_voltage, 0.0008));
    return excitations;
}

std::vector<TriggerConfig> make_triggers() {
    return {{TriggerMode::TimeDelay, 8e-6}};
}

TerminationPolicy probe_policy() {
    auto policy = TerminationPolicy::defaults();
    policy.max_steps = 8;
    policy.enable_velocity_check = false;
    policy.enable_bound_check = false;
    return policy;
}

double maximum_cap_voltage(const MultiStageResult& result,
                           const Candidate& candidate) {
    double peak = std::max(candidate.stage0_voltage, candidate.stage1_voltage);
    for (const auto& step : result.history) {
        for (double voltage : step.cap_voltages) peak = std::max(peak, voltage);
    }
    return peak;
}

double maximum_filament_temperature(const MultiStageResult& result) {
    double peak = 0.0;
    for (const auto& step : result.history) {
        for (double temperature : step.state.filament_temperatures)
            peak = std::max(peak, temperature);
    }
    return peak;
}

Metrics metrics(const MultiStageResult& result, const Candidate& candidate) {
    double peak_current = 0.0;
    for (const auto& step : result.history) {
        for (double current : step.coil_currents)
            peak_current = std::max(peak_current, std::abs(current));
    }
    return {
        result.summary.muzzle_velocity,
        peak_current,
        maximum_cap_voltage(result, candidate),
        maximum_filament_temperature(result),
    };
}

Metrics run_cpu(const Candidate& candidate, OptimizationLevel level,
                bool thermal) {
    MultiStageSim<EulerStepper> simulation(
        make_coils(), make_armature(), make_excitations(candidate), make_triggers(),
        1e-6, thermal, level);
    simulation.run(probe_policy());
    return metrics(simulation.result(), candidate);
}

Metrics run_gpu_single(const Candidate& candidate, bool thermal) {
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> simulation(
        make_coils(), make_armature(), make_excitations(candidate), make_triggers(),
        1e-6, thermal, GpuOptLevel::Full, backend);
    const auto start = std::chrono::steady_clock::now();
    simulation.run(probe_policy());
    const auto stop = std::chrono::steady_clock::now();
    REQUIRE(simulation.execution_report().gpu_executed);
    REQUIRE(simulation.execution_report().backend != BackendMode::Fallback);
    const auto& report = simulation.execution_report();
    std::cout << "BASELINE_THERMAL_REPORT row=" << candidate.id
              << " host_wall_ms="
              << std::chrono::duration<double, std::milli>(stop - start).count()
              << " gpu_time_ms=" << report.gpu_time_ms
              << " transfer_time_ms=" << report.transfer_time_ms
              << " solver_time_ms=" << report.solver_time_ms
              << " thermal_time_ms=" << report.thermal_time_ms
              << " backend=" << report.backend
              << " solver=" << report.solver
              << " precision=" << report.precision
              << " thermal=" << report.thermal
              << " device_id=" << report.device_id
              << " gpu_executed=" << (report.gpu_executed ? "true" : "false") << '\n';
    return metrics(simulation.result(), candidate);
}

TimedBatch run_gpu_batch(const std::vector<Candidate>& candidates) {
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    backend.enable_profiling = true;
    SimBatch<EulerStepper> batch(
        make_coils(), make_armature(), static_cast<int>(candidates.size()), 1e-6,
        backend);
    for (std::size_t row = 0; row < candidates.size(); ++row) {
        batch.set_excitations(static_cast<int>(row),
                              make_excitations(candidates[row]), make_triggers());
    }
    const auto start = std::chrono::steady_clock::now();
    batch.run(probe_policy());
    const auto stop = std::chrono::steady_clock::now();

    const auto& report = batch.execution_report();
    REQUIRE(report.gpu_executed);
    REQUIRE(report.backend == BackendMode::Direct);
    REQUIRE(report.backend != BackendMode::Fallback);
    REQUIRE(report.gpu_time_ms > 0.0);
    REQUIRE(report.transfer_time_ms >= 0.0);

    TimedBatch output;
    output.host_wall_ms = std::chrono::duration<double, std::milli>(stop - start).count();
    output.report = report;
    output.metrics.reserve(candidates.size());
    for (std::size_t row = 0; row < candidates.size(); ++row)
        output.metrics.push_back(metrics(batch.result(static_cast<int>(row)), candidates[row]));
    return output;
}

void require_close(const char* name, double actual, double expected,
                   gpu_test::NumericalTolerance tolerance) {
    INFO(name << " actual=" << std::setprecision(17) << actual
              << " expected=" << expected);
    CHECK(gpu_test::numerically_equal(actual, expected, tolerance));
}

void require_cuda_peak_current_baseline(const char* path, int row,
                                        double actual, double expected,
                                        gpu_test::NumericalTolerance tolerance) {
    const double relative_shortfall = (expected - actual) / std::abs(expected);
    INFO(path << " row=" << row << " CUDA peak current="
         << std::setprecision(17) << actual << " CPU Full peak current=" << expected
         << " relative shortfall=" << relative_shortfall * 100.0 << "%");
    std::cout << "BASELINE_PEAK_CURRENT_CLASSIFICATION path=" << path
              << " row=" << row << " actual=" << std::setprecision(17) << actual
              << " expected=" << expected
              << " relative_shortfall=" << relative_shortfall
              << " normal_tolerance_match="
              << (gpu_test::numerically_equal(actual, expected, tolerance) ? "true" : "false")
              << " baseline_band=0.05..0.08\n";
    CHECK_FALSE(gpu_test::numerically_equal(actual, expected, tolerance));
    CHECK(relative_shortfall >= kCudaPeakCurrentBaselineMinRelativeShortfall);
    CHECK(relative_shortfall <= kCudaPeakCurrentBaselineMaxRelativeShortfall);
}

void print_metrics(const char* path, int batch_size, int row, const Candidate& candidate,
                   const Metrics& value) {
    std::cout << "BASELINE_METRIC path=" << path
              << " batch_size=" << batch_size << " row=" << row
              << " candidate_id=" << candidate.id
              << " stage0_voltage=" << std::setprecision(17) << candidate.stage0_voltage
              << " stage1_voltage=" << candidate.stage1_voltage
              << " muzzle_velocity=" << value.muzzle_velocity
              << " peak_coil_current=" << value.peak_coil_current
              << " peak_cap_voltage=" << value.peak_cap_voltage
              << " max_filament_temperature=" << value.max_filament_temperature << '\n';
}

} // namespace

TEST_CASE("CUDA optimization numerical baseline uses real SimBatch execution" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    const auto candidates = make_candidates();
    const std::vector<int> representative_rows = {0, 7, 31, 127};
    const auto tolerance = gpu_test::tolerance_for(GpuOptLevel::Full);

    std::cout << std::setprecision(17);
    std::cout << "BASELINE_SCHEMA version=1 workload=two-stage-euler-fixed-geometry"
              << " dt=1e-6 max_steps=8 geometry_rows=2 geometry_axial=2"
              << " trigger_mode=time-delay trigger=8e-6"
              << " cuda_opt_level=full backend_request=direct thermal_batch=disabled\n";

    std::vector<Metrics> cpu_reference(candidates.size());
    std::vector<Metrics> cpu_full(candidates.size());
    for (int row : representative_rows) {
        cpu_reference[static_cast<std::size_t>(row)] =
            run_cpu(candidates[static_cast<std::size_t>(row)], OptimizationLevel::Reference, false);
        cpu_full[static_cast<std::size_t>(row)] =
            run_cpu(candidates[static_cast<std::size_t>(row)], OptimizationLevel::Full, false);
        print_metrics("cpu_reference", 0, row, candidates[static_cast<std::size_t>(row)],
                      cpu_reference[static_cast<std::size_t>(row)]);
        print_metrics("cpu_full", 0, row, candidates[static_cast<std::size_t>(row)],
                      cpu_full[static_cast<std::size_t>(row)]);
        require_close("CPU Reference vs Full muzzle velocity",
                      cpu_full[static_cast<std::size_t>(row)].muzzle_velocity,
                      cpu_reference[static_cast<std::size_t>(row)].muzzle_velocity, tolerance);
        require_close("CPU Reference vs Full peak current",
                      cpu_full[static_cast<std::size_t>(row)].peak_coil_current,
                      cpu_reference[static_cast<std::size_t>(row)].peak_coil_current, tolerance);
        require_close("CPU Reference vs Full peak capacitor voltage",
                      cpu_full[static_cast<std::size_t>(row)].peak_cap_voltage,
                      cpu_reference[static_cast<std::size_t>(row)].peak_cap_voltage, tolerance);
        require_close("CPU Reference vs Full maximum filament temperature",
                      cpu_full[static_cast<std::size_t>(row)].max_filament_temperature,
                      cpu_reference[static_cast<std::size_t>(row)].max_filament_temperature,
                      tolerance);
    }

    for (int batch_size : {1, 8, 32, 128}) {
        std::vector<Candidate> batch_candidates(candidates.begin(),
                                                candidates.begin() + batch_size);
        auto batch = run_gpu_batch(batch_candidates);
        const auto& report = batch.report;
        std::cout << "BASELINE_BATCH batch_size=" << batch_size
                  << " host_wall_ms=" << batch.host_wall_ms
                  << " gpu_time_ms=" << report.gpu_time_ms
                  << " transfer_time_ms=" << report.transfer_time_ms
                  << " solver_time_ms=" << report.solver_time_ms
                  << " thermal_time_ms=" << report.thermal_time_ms
                  << " backend=" << report.backend
                  << " solver=" << report.solver
                  << " precision=" << report.precision
                  << " thermal=" << report.thermal
                  << " device_id=" << report.device_id
                  << " threads_per_block=" << report.threads_per_block
                  << " gpu_executed=" << (report.gpu_executed ? "true" : "false") << '\n';

        for (int row : representative_rows) {
            if (row >= batch_size) continue;
            const auto& gpu = batch.metrics[static_cast<std::size_t>(row)];
            print_metrics("cuda_batch_full", batch_size, row,
                          candidates[static_cast<std::size_t>(row)], gpu);
            // The CPU candidate is indexed by the same input row. This is a
            // direct ordering check in addition to the numerical comparison.
            CHECK(batch_candidates[static_cast<std::size_t>(row)].id == row);
            require_close("CUDA Full vs CPU Full muzzle velocity",
                          gpu.muzzle_velocity, cpu_full[static_cast<std::size_t>(row)].muzzle_velocity,
                          tolerance);
            require_cuda_peak_current_baseline(
                "cuda_batch_full_vs_cpu_full", row, gpu.peak_coil_current,
                cpu_full[static_cast<std::size_t>(row)].peak_coil_current, tolerance);
            require_close("CUDA Full vs CPU Full peak capacitor voltage",
                          gpu.peak_cap_voltage, cpu_full[static_cast<std::size_t>(row)].peak_cap_voltage,
                          tolerance);
            require_close("CUDA Full vs CPU Full maximum filament temperature",
                          gpu.max_filament_temperature,
                          cpu_full[static_cast<std::size_t>(row)].max_filament_temperature,
                          tolerance);
        }
    }

    std::cout << "BASELINE_THERMAL comparison=cpu_reference,cpu_full,cuda_full"
              << " execution=GpuMultiStageSim representative_rows=0,7,31,127\n";
    for (int row : representative_rows) {
        const auto candidate = candidates[static_cast<std::size_t>(row)];
        const auto cpu_ref = run_cpu(candidate, OptimizationLevel::Reference, true);
        const auto cpu_opt = run_cpu(candidate, OptimizationLevel::Full, true);
        const auto cuda_opt = run_gpu_single(candidate, true);
        print_metrics("cpu_reference_thermal", 0, row, candidate, cpu_ref);
        print_metrics("cpu_full_thermal", 0, row, candidate, cpu_opt);
        print_metrics("cuda_full_thermal", 1, row, candidate, cuda_opt);
        require_close("thermal CPU Reference vs Full muzzle velocity",
                      cpu_opt.muzzle_velocity, cpu_ref.muzzle_velocity, tolerance);
        require_close("thermal CPU Reference vs Full peak coil current",
                      cpu_opt.peak_coil_current, cpu_ref.peak_coil_current, tolerance);
        require_close("thermal CPU Reference vs Full peak capacitor voltage",
                      cpu_opt.peak_cap_voltage, cpu_ref.peak_cap_voltage, tolerance);
        require_close("thermal CPU Reference vs Full maximum filament temperature",
                      cpu_opt.max_filament_temperature, cpu_ref.max_filament_temperature,
                      tolerance);
        require_close("thermal CUDA Full vs CPU Full muzzle velocity",
                      cuda_opt.muzzle_velocity, cpu_opt.muzzle_velocity, tolerance);
        require_cuda_peak_current_baseline(
            "cuda_full_thermal_vs_cpu_full", row, cuda_opt.peak_coil_current,
            cpu_opt.peak_coil_current, tolerance);
        require_close("thermal CUDA Full vs CPU Full peak capacitor voltage",
                      cuda_opt.peak_cap_voltage, cpu_opt.peak_cap_voltage, tolerance);
        require_close("thermal CUDA Full vs CPU Full max filament temperature",
                      cuda_opt.max_filament_temperature, cpu_opt.max_filament_temperature, tolerance);
    }
}
