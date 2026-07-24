/**
 * @file bench_gpu_engine.cu
 * @brief Reproducible CUDA execution-boundary benchmark.
 *
 * This is a measurement executable, not a pass/fail performance test. It
 * reports setup, first-step/capture-inclusive, replay-only, and steady-state
 * samples separately so CPU fallback work cannot be reported as GPU speedup.
 */

#include "coilgun/components/armature.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/termination.hpp"

#include <cuda_runtime_api.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

using coilgun::components::Armature;
using coilgun::components::DrivingCoil;
using coilgun::physics::ALUMINUM;
using coilgun::physics::COPPER;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::MultiStageSim;
using coilgun::simulation::TerminationPolicy;
using coilgun::simulation::cuda::BackendMode;
using coilgun::simulation::cuda::ExecutionReport;
using coilgun::simulation::cuda::GpuEngine;
using coilgun::simulation::cuda::GpuEngineState;
using coilgun::simulation::cuda::GpuExecutionConfig;
using coilgun::simulation::cuda::GpuGeometryInput;
using coilgun::simulation::cuda::PrecisionMode;
using coilgun::simulation::cuda::SolverMode;
using coilgun::simulation::cuda::ThermalMode;

constexpr int kBenchmarkRepeats = 3;
constexpr int kWarmupSteps = 5;
constexpr int kMeasuredSteps = 10;

struct Case {
    std::string name;
    int stages;
    int axial_filaments;
    int radial_filaments;
    int steady_steps;
    bool thermal;
};

struct Request {
    BackendMode backend;
    SolverMode solver;
    ThermalMode thermal;
    PrecisionMode precision = PrecisionMode::Full;
    int batch_size;
    bool runtime_mask_change = false;
};

struct PhaseSample {
    std::string phase;
    int repeat = 0;
    int iterations = 0;
    double wall_ms = 0.0;
    double gpu_ms = 0.0;
    double solver_ms = 0.0;
    double thermal_ms = 0.0;
    double transfer_ms = 0.0;
    int graph_rebuilds = 0;
    int graph_rebuild_count = 0;
    int fallback_count = 0;
    int fallback_total = 0;
    bool gpu_executed = false;
    bool finite = false;
    ExecutionReport report;
};

struct CpuPhaseSample {
    std::string phase;
    int iterations = 0;
    double wall_ms = 0.0;
};

std::vector<DrivingCoil> make_coils(int stages) {
    std::vector<DrivingCoil> coils;
    coils.reserve(static_cast<std::size_t>(stages));
    for (int stage = 0; stage < stages; ++stage) {
        coils.emplace_back(0.01, 0.03, 0.05, 150,
                           COPPER.resistivity_ref, 1e-6, 0.7,
                           0.015 + 0.07 * stage);
    }
    return coils;
}

Armature make_armature(const Case& test_case) {
    return Armature(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref,
                    ALUMINUM.density, 0.0, 0.120,
                    test_case.axial_filaments, test_case.radial_filaments, 0.0);
}

std::vector<std::unique_ptr<coilgun::simulation::Excitation>> make_excitations(int stages) {
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excitations;
    excitations.reserve(static_cast<std::size_t>(stages));
    for (int stage = 0; stage < stages; ++stage)
        excitations.push_back(std::make_unique<CrowbarExcitation>(
            480.0 - 40.0 * static_cast<double>(stage), 0.001));
    return excitations;
}

TerminationPolicy benchmark_policy(int max_steps) {
    TerminationPolicy policy;
    policy.max_steps = max_steps;
    policy.enable_velocity_check = false;
    policy.enable_bound_check = false;
    return policy;
}

double elapsed_ms(std::chrono::steady_clock::time_point start,
                  std::chrono::steady_clock::time_point stop) {
    return std::chrono::duration<double, std::milli>(stop - start).count();
}

GpuGeometryInput make_geometry(const Case& test_case, const Armature& armature,
                               bool thermal_enabled) {
    const auto coils = make_coils(test_case.stages);
    const auto stages = static_cast<std::size_t>(test_case.stages);
    const auto filaments = static_cast<std::size_t>(armature.total_filaments());
    GpuGeometryInput geometry;
    geometry.n_stages = stages;
    geometry.n_filaments = filaments;
    geometry.thermal_enabled = thermal_enabled;
    geometry.stage_geometry.resize(stages);
    geometry.stage_inner_radii.resize(stages);
    geometry.stage_outer_radii.resize(stages);
    geometry.stage_lengths.resize(stages);
    geometry.stage_turns.resize(stages);
    geometry.stage_positions.resize(stages);
    geometry.stage_resistances.resize(stages);
    geometry.stage_inductances.resize(stages);
    for (std::size_t stage = 0; stage < stages; ++stage) {
        const auto& coil = coils[stage];
        geometry.stage_geometry[stage] = coil.mean_radius();
        geometry.stage_inner_radii[stage] = coil.inner_radius();
        geometry.stage_outer_radii[stage] = coil.outer_radius();
        geometry.stage_lengths[stage] = coil.length();
        geometry.stage_turns[stage] = coil.turns();
        geometry.stage_positions[stage] = coil.position();
        geometry.stage_resistances[stage] = coil.resistance();
        geometry.stage_inductances[stage] = coil.self_inductance();
    }

    geometry.filament_geometry.resize(filaments);
    geometry.filament_inner_radii.resize(filaments);
    geometry.filament_outer_radii.resize(filaments);
    geometry.filament_lengths.resize(filaments);
    geometry.filament_positions.resize(filaments);
    const auto filament_length = armature.length() / armature.axial_filaments();
    for (std::size_t filament = 0; filament < filaments; ++filament) {
        const int radial = static_cast<int>(filament % armature.radial_filaments()) + 1;
        const int axial = static_cast<int>(filament / armature.radial_filaments()) + 1;
        geometry.filament_geometry[filament] = armature.filament_mean_radius(radial);
        geometry.filament_inner_radii[filament] = armature.filament_inner_radius(radial);
        geometry.filament_outer_radii[filament] = armature.filament_outer_radius(radial);
        geometry.filament_lengths[filament] = filament_length;
        geometry.filament_positions[filament] = armature.filament_axial_position(axial);
    }
    geometry.filament_resistances = armature.resistances();
    geometry.filament_inductances = armature.inductances();
    return geometry;
}

GpuEngineState make_state(const Case& test_case, const Armature& armature,
                          int batch_size, bool thermal_enabled) {
    const auto stages = static_cast<std::size_t>(test_case.stages);
    const auto filaments = static_cast<std::size_t>(armature.total_filaments());
    const auto batch = static_cast<std::size_t>(batch_size);
    const auto dimension = stages + filaments;
    GpuEngineState state;
    state.currents.assign(batch * dimension, 0.0);
    state.m1.assign(batch * stages * filaments, 0.0);
    state.dm1.assign(batch * stages * filaments, 0.0);
    state.active_mask.assign(batch, 1);
    state.trigger_mask.assign(batch * stages, 1);
    state.stage_mask.assign(batch * stages, 1);
    state.mutual_stage_mask.assign(batch * stages, 1);
    state.stage_voltages.resize(batch * stages);
    for (std::size_t b = 0; b < batch; ++b) {
        for (std::size_t stage = 0; stage < stages; ++stage)
            state.stage_voltages[b * stages + stage] =
                480.0 - 40.0 * static_cast<double>(stage);
    }
    state.velocity.assign(batch, -armature.velocity());
    state.position.assign(batch, 0.0);
    state.dt = 1e-6;
    state.mass = armature.mass();

    if (thermal_enabled) {
        const auto& masses = armature.masses();
        const auto& resistances = armature.resistances();
        state.temperatures.assign(batch * filaments, coilgun::physics::T_REFERENCE);
        state.filament_masses.reserve(batch * filaments);
        state.reference_resistances.reserve(batch * filaments);
        state.resistances.reserve(batch * filaments);
        state.filament_materials.assign(batch * filaments, 0);
        for (std::size_t b = 0; b < batch; ++b) {
            state.filament_masses.insert(state.filament_masses.end(), masses.begin(), masses.end());
            state.reference_resistances.insert(state.reference_resistances.end(),
                                               resistances.begin(), resistances.end());
            state.resistances.insert(state.resistances.end(), resistances.begin(), resistances.end());
        }
        state.material_density = ALUMINUM.density;
    }
    return state;
}

const char* phase_name(int phase) {
    switch (phase) {
    case 0: return "setup";
    case 1: return "first-step/capture-inclusive";
    case 2: return "replay-only";
    case 3: return "steady-state";
    default: return "unknown";
    }
}

PhaseSample make_sample(const std::string& phase, int repeat, int iterations, double wall_ms,
                        const ExecutionReport& before, const ExecutionReport& after) {
    PhaseSample sample;
    sample.phase = phase;
    sample.repeat = repeat;
    sample.iterations = iterations;
    sample.wall_ms = wall_ms;
    sample.gpu_ms = after.gpu_time_ms - before.gpu_time_ms;
    sample.solver_ms = after.solver_time_ms - before.solver_time_ms;
    sample.thermal_ms = after.thermal_time_ms - before.thermal_time_ms;
    sample.transfer_ms = after.transfer_time_ms - before.transfer_time_ms;
    sample.graph_rebuilds = after.graph_rebuild_count - before.graph_rebuild_count;
    sample.graph_rebuild_count = after.graph_rebuild_count;
    sample.fallback_count = after.fallback_count - before.fallback_count;
    sample.fallback_total = after.fallback_count;
    sample.gpu_executed = after.gpu_executed &&
                          after.backend != BackendMode::Fallback &&
                          sample.fallback_count == 0;
    sample.finite = std::isfinite(sample.wall_ms) && std::isfinite(sample.gpu_ms) &&
                    std::isfinite(sample.solver_ms) && std::isfinite(sample.thermal_ms) &&
                    std::isfinite(sample.transfer_ms);
    sample.report = after;
    return sample;
}

std::vector<PhaseSample> benchmark_gpu(const Case& test_case, const Request& request) {
    std::vector<PhaseSample> samples;
    samples.reserve(static_cast<std::size_t>(kBenchmarkRepeats) * 5);
    const bool thermal_enabled = request.thermal != ThermalMode::Disabled;

    for (int repeat = 0; repeat < kBenchmarkRepeats; ++repeat) {
        const auto setup_start = std::chrono::steady_clock::now();
        const auto armature = make_armature(test_case);
        GpuExecutionConfig config;
        config.backend = request.backend;
        config.solver = request.solver;
        config.precision = request.precision;
        config.thermal = request.thermal;
        config.deterministic = request.backend != BackendMode::Persistent;
        auto engine = std::make_unique<GpuEngine>(
            make_geometry(test_case, armature, thermal_enabled),
            make_state(test_case, armature, request.batch_size, thermal_enabled), config);
        const auto setup_stop = std::chrono::steady_clock::now();

        const ExecutionReport empty_report;
        samples.push_back(make_sample("setup", repeat, 0,
                                      elapsed_ms(setup_start, setup_stop),
                                      empty_report, engine->report()));

        std::vector<std::uint8_t> enabled_mask(
            static_cast<std::size_t>(request.batch_size * test_case.stages), 1);
        std::vector<std::uint8_t> disabled_mask = enabled_mask;
        for (int batch = 0; batch < request.batch_size; ++batch)
            disabled_mask[static_cast<std::size_t>(batch * test_case.stages +
                                                   test_case.stages - 1)] = 0;
        int mask_iteration = 0;
        auto update_runtime_mask = [&]() {
            if (!request.runtime_mask_change) return;
            engine->set_mutual_stage_mask(
                (mask_iteration++ % 2 == 0) ? disabled_mask : enabled_mask);
        };
        auto run_phase = [&](const char* name, int iterations, bool update_mask) {
            const auto before = engine->report();
            const auto start = std::chrono::steady_clock::now();
            for (int iteration = 0; iteration < iterations; ++iteration) {
                if (update_mask) update_runtime_mask();
                engine->step();
            }
            const auto stop = std::chrono::steady_clock::now();
            samples.push_back(make_sample(name, repeat, iterations,
                                          elapsed_ms(start, stop), before,
                                          engine->report()));
        };

        run_phase(phase_name(1), 1, false);
        run_phase(phase_name(2), 1, request.runtime_mask_change);
        run_phase("warm-up", kWarmupSteps, request.runtime_mask_change);
        run_phase(phase_name(3), kMeasuredSteps, request.runtime_mask_change);
    }
    return samples;
}

std::vector<CpuPhaseSample> benchmark_cpu(const Case& test_case) {
    const auto setup_start = std::chrono::steady_clock::now();
    auto simulation = std::make_unique<MultiStageSim<EulerStepper>>(
        make_coils(test_case.stages), make_armature(test_case),
        make_excitations(test_case.stages),
        std::vector<coilgun::simulation::TriggerConfig>(
            static_cast<std::size_t>(test_case.stages - 1),
            {coilgun::simulation::TriggerMode::TimeDelay, 0.0}),
        1e-6, test_case.thermal,
        coilgun::simulation::OptimizationLevel::Reference);
    const auto setup_stop = std::chrono::steady_clock::now();

    std::vector<CpuPhaseSample> samples;
    samples.push_back({"setup", 0, elapsed_ms(setup_start, setup_stop)});
    auto run_phase = [&](const char* name, int max_steps, int iterations) {
        const auto start = std::chrono::steady_clock::now();
        simulation->run(benchmark_policy(max_steps));
        const auto stop = std::chrono::steady_clock::now();
        samples.push_back({name, iterations, elapsed_ms(start, stop)});
    };
    run_phase(phase_name(1), 1, 1);
    run_phase(phase_name(2), 2, 1);
    run_phase(phase_name(3), 2 + test_case.steady_steps,
              test_case.steady_steps);
    return samples;
}

void print_cpu_baseline(const Case& test_case, const std::vector<CpuPhaseSample>& samples) {
    for (const auto& sample : samples) {
        const double per_step = sample.iterations > 0
            ? sample.wall_ms / sample.iterations : sample.wall_ms;
        std::cout << "| " << test_case.name << " | cpu-reference | "
                  << (test_case.thermal ? "thermal" : "cold") << " | " << sample.phase
                  << " | " << sample.iterations << " | " << std::fixed << std::setprecision(3)
                  << sample.wall_ms << " | " << per_step << " |\n";
    }
}

const CpuPhaseSample* find_cpu_phase(const std::vector<CpuPhaseSample>& samples,
                                     const std::string& phase) {
    for (const auto& sample : samples)
        if (sample.phase == phase) return &sample;
    return nullptr;
}

void print_gpu_row(const Case& test_case, const Request& request,
                   const PhaseSample& sample, const std::vector<CpuPhaseSample>& cpu) {
    const auto* cpu_sample = find_cpu_phase(cpu, sample.phase);
    const double per_step = sample.iterations > 0
        ? sample.wall_ms / sample.iterations : sample.wall_ms;
    const bool comparable = sample.gpu_executed && cpu_sample != nullptr && per_step > 0.0;
    const double cpu_per_step = cpu_sample != nullptr && cpu_sample->iterations > 0
        ? cpu_sample->wall_ms / cpu_sample->iterations
        : (cpu_sample != nullptr ? cpu_sample->wall_ms : 0.0);
    const double speedup = comparable
        ? (cpu_per_step * request.batch_size) / per_step : 0.0;
    const double steps_per_second = sample.iterations > 0 && per_step > 0.0
        ? 1000.0 / per_step : 0.0;
    const auto& report = sample.report;
    std::cout << "| " << test_case.name << " | "
              << coilgun::simulation::cuda::to_string(request.backend) << " | "
              << coilgun::simulation::cuda::to_string(request.solver) << " | "
              << coilgun::simulation::cuda::to_string(request.thermal) << " | "
              << request.batch_size << " | "
              << (request.runtime_mask_change ? "yes" : "no") << " | " << sample.repeat
              << " | " << sample.phase << " | " << sample.iterations
              << " | " << coilgun::simulation::cuda::to_string(report.backend) << " | "
              << coilgun::simulation::cuda::to_string(report.solver) << " | "
              << coilgun::simulation::cuda::to_string(report.precision) << " | "
              << coilgun::simulation::cuda::to_string(report.thermal) << " | "
              << std::fixed << std::setprecision(3) << sample.wall_ms << " | " << per_step
              << " | " << steps_per_second
              << " | " << sample.gpu_ms << " | " << sample.solver_ms << " | "
              << sample.thermal_ms << " | " << sample.transfer_ms << " | "
              << sample.graph_rebuilds << " | " << sample.graph_rebuild_count
              << " | " << sample.fallback_count << " | " << sample.fallback_total
              << " | " << (comparable ? std::to_string(speedup) : "n/a")
              << " | " << (sample.gpu_executed ? "yes" : "no") << " | "
              << (sample.finite ? "yes" : "no") << " | "
              << (report.fallback_reason.empty() ? "none" : report.fallback_reason) << " |\n";
}

void print_runtime_metadata(bool cuda_available, int device_count) {
    std::cout << "CUDA available: " << (cuda_available ? "yes" : "no")
              << " (device_count=" << device_count << ")\n";
    if (cuda_available) {
        cudaDeviceProp properties{};
        const auto property_status = cudaGetDeviceProperties(&properties, 0);
        if (property_status == cudaSuccess)
            std::cout << "GPU: " << properties.name << " (compute "
                      << properties.major << "." << properties.minor << ")\n";
    }
    int driver_version = 0;
    int runtime_version = 0;
    const auto driver_status = cudaDriverGetVersion(&driver_version);
    const auto runtime_status = cudaRuntimeGetVersion(&runtime_version);
    if (driver_status == cudaSuccess)
        std::cout << "CUDA driver API: " << driver_version / 1000 << "."
                  << (driver_version % 1000) / 10 << "\n";
    else
        std::cout << "CUDA driver API: unavailable ("
                  << cudaGetErrorString(driver_status) << ")\n";
    if (runtime_status == cudaSuccess)
        std::cout << "CUDA runtime: " << runtime_version / 1000 << "."
                  << (runtime_version % 1000) / 10 << "\n";
    else
        std::cout << "CUDA runtime: unavailable ("
                  << cudaGetErrorString(runtime_status) << ")\n";
    std::cout << "Build contract: CUDA Release; setup, first-step/capture-inclusive, replay-only, "
                 "warm-up, steady-state\n";
    std::cout << "Benchmark repeats: " << kBenchmarkRepeats
              << "; warm-up steps: " << kWarmupSteps
              << "; measured steady-state steps: " << kMeasuredSteps << "\n\n";
}

} // namespace

int main() {
    int device_count = 0;
    const auto device_status = cudaGetDeviceCount(&device_count);
    const bool cuda_available = device_status == cudaSuccess && device_count > 0;
    print_runtime_metadata(cuda_available, device_count);

    const std::vector<Case> cases = {
        {"baseline", 2, 8, 4, kMeasuredSteps, false},
        {"thermal", 2, 8, 4, kMeasuredSteps, true},
    };
    std::cout << "Fixed workloads:\n";
    for (const auto& test_case : cases) {
        std::cout << "- " << test_case.name << ": stages=" << test_case.stages
                  << ", axial_filaments=" << test_case.axial_filaments
                  << ", radial_filaments=" << test_case.radial_filaments
                  << ", total_filaments="
                  << test_case.axial_filaments * test_case.radial_filaments
                  << ", thermal=" << (test_case.thermal ? "gpu" : "disabled") << "\n";
    }
    std::cout << "\n";

    std::vector<std::vector<CpuPhaseSample>> cpu_baselines;
    cpu_baselines.reserve(cases.size());
    std::cout << "CPU reference phases:\n";
    std::cout << "| Workload | Requested | Thermal | Phase | Iterations | Wall ms | Per-step ms |\n";
    std::cout << "|---|---|---|---|---:|---:|---:|\n";
    for (const auto& test_case : cases) {
        cpu_baselines.push_back(benchmark_cpu(test_case));
        print_cpu_baseline(test_case, cpu_baselines.back());
    }

    std::cout << "\nGPU execution phases:\n";
    std::cout << "| Workload | Requested backend | Requested solver | Requested thermal | Batch | Runtime mask change | Repeat | Phase | Iterations | Resolved backend | Resolved solver | Precision | Resolved thermal | Wall ms | Per-step ms | Steps/s | GPU ms | Solver ms | Thermal ms | Transfer ms | Graph rebuild delta | Graph rebuild total | Fallback delta | Fallback total | CPU/GPU speedup | GPU executed | Finite | Fallback reason |\n";
    std::cout << "|---|---|---|---|---:|---|---:|---|---:|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|\n";

    const std::vector<Request> baseline_requests = {
        {BackendMode::Direct, SolverMode::Eigen, ThermalMode::Disabled,
         PrecisionMode::Full, 32},
        {BackendMode::Graph, SolverMode::Batched, ThermalMode::Disabled,
         PrecisionMode::Full, 32},
        {BackendMode::Persistent, SolverMode::Eigen, ThermalMode::Disabled,
         PrecisionMode::Full, 1},
        {BackendMode::Fallback, SolverMode::Eigen, ThermalMode::Disabled,
         PrecisionMode::Full, 1},
        {BackendMode::Graph, SolverMode::Batched, ThermalMode::Disabled,
         PrecisionMode::Full, 32, true},
    };
    for (const auto& request : baseline_requests) {
        for (const auto& sample : benchmark_gpu(cases.front(), request))
            print_gpu_row(cases.front(), request, sample, cpu_baselines.front());
    }

    const Request thermal_request{BackendMode::Graph, SolverMode::Batched,
                                  ThermalMode::Gpu, PrecisionMode::Full, 32};
    for (const auto& sample : benchmark_gpu(cases.back(), thermal_request))
        print_gpu_row(cases.back(), thermal_request, sample, cpu_baselines.back());
    return 0;
}
