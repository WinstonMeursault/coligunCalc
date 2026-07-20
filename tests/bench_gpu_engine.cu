/**
 * @file bench_gpu_engine.cu
 * @brief CPU/GPU wall-time and execution-report benchmark for the unified engine.
 *
 * This is a measurement executable, not a pass/fail performance test. It emits
 * Markdown-table rows so runs can be copied into the benchmark record without
 * treating one machine's timings as universal guarantees.
 */

#include "coilgun/components/armature.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/simulation/cuda/sim_batch.hpp"
#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/termination.hpp"

#include <cuda_runtime_api.h>

#include <chrono>
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
using coilgun::simulation::cuda::GpuBackend;
using coilgun::simulation::cuda::GpuEngine;
using coilgun::simulation::cuda::GpuEngineState;
using coilgun::simulation::cuda::GpuExecutionConfig;
using coilgun::simulation::cuda::GpuGeometryInput;
using coilgun::simulation::cuda::GpuOptLevel;
using coilgun::simulation::cuda::PrecisionMode;
using coilgun::simulation::cuda::SolverMode;
using coilgun::simulation::cuda::ThermalMode;
using coilgun::simulation::cuda::SimBatch;

struct Case {
    std::string name;
    int stages;
    int axial_filaments;
    int radial_filaments;
    int steps;
    bool thermal;
};

struct Sample {
    double wall_ms = 0.0;
    double gpu_ms = 0.0;
    double solver_ms = 0.0;
    double thermal_ms = 0.0;
    double transfer_ms = 0.0;
    int graph_rebuilds = 0;
    std::string backend;
    std::string solver;
    std::string precision;
    std::string thermal;
    bool gpu_executed = false;
    std::string fallback;
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
        excitations.push_back(std::make_unique<CrowbarExcitation>(480.0 - 40.0 * stage, 0.001));
    return excitations;
}

TerminationPolicy benchmark_policy(int steps) {
    TerminationPolicy policy;
    policy.max_steps = steps;
    policy.enable_velocity_check = false;
    policy.enable_bound_check = false;
    return policy;
}

double elapsed_ms(std::chrono::steady_clock::time_point start,
                  std::chrono::steady_clock::time_point stop) {
    return std::chrono::duration<double, std::milli>(stop - start).count();
}

Sample benchmark_cpu(const Case& test_case) {
    const auto start = std::chrono::steady_clock::now();
    MultiStageSim<EulerStepper> simulation(
        make_coils(test_case.stages), make_armature(test_case),
        make_excitations(test_case.stages),
        std::vector<coilgun::simulation::TriggerConfig>(
            static_cast<std::size_t>(test_case.stages - 1),
            {coilgun::simulation::TriggerMode::TimeDelay, 0.0}),
        1e-6, test_case.thermal,
        coilgun::simulation::OptimizationLevel::Reference);
    simulation.run(benchmark_policy(test_case.steps));
    const auto stop = std::chrono::steady_clock::now();

    Sample sample;
    sample.wall_ms = elapsed_ms(start, stop);
    sample.backend = "cpu-reference";
    sample.solver = "eigen";
    sample.precision = "reference";
    sample.thermal = test_case.thermal ? "cpu" : "disabled";
    return sample;
}

Sample benchmark_gpu(const Case& test_case, BackendMode requested_backend,
                     int batch_size = 1) {
    GpuBackend backend;
    backend.backend = requested_backend;
    backend.use_persistent = requested_backend == BackendMode::Persistent;
    backend.max_batch_sims = 256;

    SimBatch<EulerStepper> batch(make_coils(test_case.stages), make_armature(test_case),
                                 batch_size, 1e-6, backend);
    const auto triggers = std::vector<coilgun::simulation::TriggerConfig>(
        static_cast<std::size_t>(test_case.stages - 1),
        {coilgun::simulation::TriggerMode::TimeDelay, 0.0});
    for (int id = 0; id < batch_size; ++id)
        batch.set_excitations(id, make_excitations(test_case.stages), triggers);

    const auto start = std::chrono::steady_clock::now();
    batch.run(benchmark_policy(test_case.steps));
    const auto stop = std::chrono::steady_clock::now();
    const auto& report = batch.execution_report();

    Sample sample;
    sample.wall_ms = elapsed_ms(start, stop);
    sample.gpu_ms = report.gpu_time_ms;
    sample.solver_ms = report.solver_time_ms;
    sample.thermal_ms = report.thermal_time_ms;
    sample.transfer_ms = report.transfer_time_ms;
    sample.graph_rebuilds = report.graph_rebuild_count;
    sample.backend = coilgun::simulation::cuda::to_string(report.backend);
    sample.solver = coilgun::simulation::cuda::to_string(report.solver);
    sample.precision = coilgun::simulation::cuda::to_string(report.precision);
    sample.thermal = coilgun::simulation::cuda::to_string(report.thermal);
    sample.gpu_executed = report.gpu_executed;
    sample.fallback = report.fallback_reason;
    return sample;
}

GpuGeometryInput make_engine_geometry(const Case& test_case, const Armature& armature) {
    GpuGeometryInput geometry;
    geometry.n_stages = static_cast<std::size_t>(test_case.stages);
    geometry.n_filaments = static_cast<std::size_t>(armature.total_filaments());
    geometry.thermal_enabled = true;
    geometry.stage_geometry.resize(1);
    geometry.stage_inner_radii = {0.01};
    geometry.stage_outer_radii = {0.03};
    geometry.stage_lengths = {0.05};
    geometry.stage_turns = {150};
    geometry.stage_positions = {0.015};
    geometry.stage_resistances = {DrivingCoil(0.01, 0.03, 0.05, 150,
                                              COPPER.resistivity_ref, 1e-6, 0.7).resistance()};
    geometry.stage_inductances = {DrivingCoil(0.01, 0.03, 0.05, 150,
                                              COPPER.resistivity_ref, 1e-6, 0.7).self_inductance()};
    const auto filaments = geometry.n_filaments;
    geometry.filament_geometry.resize(filaments);
    geometry.filament_inner_radii.resize(filaments);
    geometry.filament_outer_radii.resize(filaments);
    geometry.filament_lengths.resize(filaments);
    geometry.filament_positions.resize(filaments);
    for (std::size_t k = 0; k < filaments; ++k) {
        const int radial = static_cast<int>(k % armature.radial_filaments()) + 1;
        const int axial = static_cast<int>(k / armature.radial_filaments()) + 1;
        geometry.filament_geometry[k] = armature.filament_mean_radius(radial);
        geometry.filament_inner_radii[k] = armature.filament_inner_radius(radial);
        geometry.filament_outer_radii[k] = armature.filament_outer_radius(radial);
        geometry.filament_lengths[k] = armature.length() / armature.axial_filaments();
        geometry.filament_positions[k] = armature.filament_axial_position(axial);
    }
    geometry.filament_resistances = armature.resistances();
    geometry.filament_inductances = armature.inductances();
    return geometry;
}

GpuEngineState make_engine_thermal_state(const Case& test_case, const Armature& armature) {
    const auto stages = static_cast<std::size_t>(test_case.stages);
    const auto filaments = static_cast<std::size_t>(armature.total_filaments());
    GpuEngineState state;
    state.currents.assign(stages + filaments, 0.0);
    state.m1.assign(stages * filaments, 0.0);
    state.dm1.assign(stages * filaments, 0.0);
    state.active_mask = {1};
    state.trigger_mask.assign(stages, 1);
    state.stage_mask.assign(stages, 1);
    state.mutual_stage_mask.assign(stages, 1);
    state.stage_voltages.assign(stages, 480.0);
    state.velocity = {-armature.velocity()};
    state.position = {0.0};
    state.temperatures.assign(filaments, coilgun::physics::T_REFERENCE);
    state.filament_masses = armature.masses();
    state.reference_resistances = armature.resistances();
    state.filament_materials.assign(filaments, 0);
    state.resistances = armature.resistances();
    state.material_density = ALUMINUM.density;
    state.dt = 1e-6;
    state.mass = armature.mass();
    return state;
}

Sample benchmark_engine_thermal(const Case& test_case, ThermalMode thermal,
                                PrecisionMode precision) {
    auto armature = make_armature(test_case);
    GpuExecutionConfig config;
    config.backend = thermal == ThermalMode::Cpu ? BackendMode::Fallback : BackendMode::Direct;
    config.solver = SolverMode::Eigen;
    config.precision = precision;
    config.thermal = thermal;
    config.deterministic = true;
    GpuEngine engine(make_engine_geometry(test_case, armature),
                     make_engine_thermal_state(test_case, armature), config);
    const auto start = std::chrono::steady_clock::now();
    engine.run(static_cast<std::size_t>(test_case.steps));
    const auto stop = std::chrono::steady_clock::now();
    const auto& report = engine.report();
    Sample sample;
    sample.wall_ms = elapsed_ms(start, stop);
    sample.gpu_ms = report.gpu_time_ms;
    sample.solver_ms = report.solver_time_ms;
    sample.thermal_ms = report.thermal_time_ms;
    sample.transfer_ms = report.transfer_time_ms;
    sample.graph_rebuilds = report.graph_rebuild_count;
    sample.backend = coilgun::simulation::cuda::to_string(report.backend);
    sample.solver = coilgun::simulation::cuda::to_string(report.solver);
    sample.precision = coilgun::simulation::cuda::to_string(report.precision);
    sample.thermal = coilgun::simulation::cuda::to_string(report.thermal);
    sample.gpu_executed = report.gpu_executed;
    sample.fallback = report.fallback_reason;
    return sample;
}

void print_row(const std::string& workload, const std::string& requested,
               int batch_size, const Sample& sample, double cpu_wall_ms) {
    const double speedup = sample.wall_ms > 0.0 ? cpu_wall_ms / sample.wall_ms : 0.0;
    std::cout << "| " << workload << " | " << requested << " | " << batch_size
              << " | " << sample.backend << " | " << sample.solver
              << " | " << sample.precision << " | " << sample.thermal
              << " | " << std::fixed << std::setprecision(3)
              << sample.wall_ms << " | " << sample.gpu_ms
              << " | " << sample.solver_ms << " | " << sample.thermal_ms
              << " | " << sample.transfer_ms << " | " << sample.graph_rebuilds
              << " | " << speedup << " | " << (sample.gpu_executed ? "yes" : "no")
              << " | " << (sample.fallback.empty() ? "" : sample.fallback) << " |\n";
}

} // namespace

int main() {
    int device_count = 0;
    const bool cuda_available = cudaGetDeviceCount(&device_count) == cudaSuccess && device_count > 0;
    const std::vector<Case> cases = {
        {"small-single", 1, 5, 2, 8, false},
        {"medium-multi", 2, 8, 4, 8, false},
        {"large-single", 1, 16, 8, 4, false},
        {"thermal-single", 1, 8, 4, 8, true},
    };

    std::cout << "GPU available: " << (cuda_available ? "yes" : "no") << "\n\n";
    std::cout << "| Workload | Requested | Batch | Resolved backend | Solver | Precision | Thermal | Wall ms | GPU ms | Solver ms | Thermal ms | Transfer ms | Graph rebuilds | CPU/GPU speedup | GPU executed | Fallback |\n";
    std::cout << "|---|---|---:|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|\n";

    for (const auto& test_case : cases) {
        const auto cpu = benchmark_cpu(test_case);
        print_row(test_case.name, "CPU", 1, cpu, cpu.wall_ms);
        if (!cuda_available) continue;

        const std::vector<BackendMode> backends = {
            BackendMode::Direct, BackendMode::Graph,
            BackendMode::Persistent, BackendMode::Fallback};
        for (const auto backend : backends) {
            const auto gpu = benchmark_gpu(test_case, backend);
            print_row(test_case.name, coilgun::simulation::cuda::to_string(backend),
                      1, gpu, cpu.wall_ms);
        }
    }

    if (cuda_available) {
        const Case batch_case{"batch-medium", 2, 8, 4, 8, false};
        const auto cpu = benchmark_cpu(batch_case);
        for (const int batch_size : {1, 8, 32, 128}) {
            const auto gpu = benchmark_gpu(batch_case, BackendMode::Direct, batch_size);
            print_row(batch_case.name, "direct", batch_size, gpu, cpu.wall_ms * batch_size);
        }

        const Case thermal_case{"thermal-engine", 1, 8, 4, 8, true};
        const auto cpu_thermal = benchmark_cpu(thermal_case);
        print_row(thermal_case.name, "CPU thermal", 1, cpu_thermal, cpu_thermal.wall_ms);
        const auto full_thermal = benchmark_engine_thermal(
            thermal_case, ThermalMode::Gpu, PrecisionMode::Full);
        print_row(thermal_case.name, "GPU Full thermal", 1, full_thermal,
                  cpu_thermal.wall_ms);
        const auto aggressive_thermal = benchmark_engine_thermal(
            thermal_case, ThermalMode::Gpu, PrecisionMode::Aggressive);
        print_row(thermal_case.name, "GPU Aggressive thermal", 1,
                  aggressive_thermal, cpu_thermal.wall_ms);
    }
    return 0;
}
