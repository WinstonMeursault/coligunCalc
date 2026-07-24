/**
 * @file bench_cpu_sim.cpp
 * @brief Machine-readable CPU simulation phase benchmark.
 *
 * This is a measurement executable, not a pass/fail performance test. It
 * separates construction, first step, and warmed steady-state timing.
 */

#include "coilgun/components/armature.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/simulation/cpu_phase_timing.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/time_stepper.hpp"

#include <Eigen/Core>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using coilgun::components::Armature;
using coilgun::components::DrivingCoil;
using coilgun::physics::ALUMINUM;
using coilgun::physics::COPPER;
using coilgun::simulation::CpuPhase;
using coilgun::simulation::CpuPhaseTiming;
using coilgun::simulation::CpuPhaseTimingCollector;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::MultiStageSim;

using Clock = std::chrono::steady_clock;

struct Workload {
    const char* name;
    int stages;
    int axial_filaments;
    int radial_filaments;
    bool thermal;
    int steady_steps;
};

struct Sample {
    const char* workload;
    bool thermal;
    const char* phase;
    int iterations;
    double wall_ms;
    CpuPhaseTiming timing;
    bool finite;
    double position;
    double velocity;
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

Armature make_armature(const Workload& workload) {
    return Armature(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref,
                    ALUMINUM.density, 0.0, 0.120,
                    workload.axial_filaments, workload.radial_filaments, 0.0);
}

std::vector<std::unique_ptr<coilgun::simulation::Excitation>>
make_excitations(int stages) {
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excitations;
    excitations.reserve(static_cast<std::size_t>(stages));
    for (int stage = 0; stage < stages; ++stage) {
        excitations.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
            [](double) { return 1.0; }));
    }
    return excitations;
}

std::vector<coilgun::simulation::TriggerConfig> make_triggers(int stages) {
    return std::vector<coilgun::simulation::TriggerConfig>(
        static_cast<std::size_t>(stages - 1),
        {coilgun::simulation::TriggerMode::TimeDelay,
         std::numeric_limits<double>::infinity()});
}

double elapsed_ms(Clock::time_point start, Clock::time_point stop) {
    return std::chrono::duration<double, std::milli>(stop - start).count();
}

bool state_is_finite(const MultiStageSim<EulerStepper>& simulation) {
    const auto& state = simulation.state();
    return state.currents.allFinite() &&
           (!state.filament_temperatures.size() || state.filament_temperatures.allFinite()) &&
           std::isfinite(state.arm_position) && std::isfinite(state.arm_velocity);
}

Sample run_timed_phase(MultiStageSim<EulerStepper>& simulation,
                       const Workload& workload, const char* phase,
                       int iterations) {
    CpuPhaseTiming timing;
    const auto start = Clock::now();
    {
        CpuPhaseTimingCollector collector(timing);
        for (int iteration = 0; iteration < iterations; ++iteration)
            simulation.step();
    }
    const auto stop = Clock::now();
    return {workload.name, workload.thermal, phase, iterations,
            elapsed_ms(start, stop), timing, state_is_finite(simulation),
            simulation.state().arm_position, simulation.state().arm_velocity};
}

std::vector<Sample> benchmark(const Workload& workload) {
    const auto setup_start = Clock::now();
    auto simulation = std::make_unique<MultiStageSim<EulerStepper>>(
        make_coils(workload.stages), make_armature(workload),
        make_excitations(workload.stages), make_triggers(workload.stages),
        1e-8, workload.thermal,
        coilgun::simulation::OptimizationLevel::Reference);
    const auto setup_stop = Clock::now();

    CpuPhaseTiming empty_timing;
    std::vector<Sample> samples;
    samples.push_back({workload.name, workload.thermal, "setup", 0,
                       elapsed_ms(setup_start, setup_stop), empty_timing,
                       state_is_finite(*simulation), simulation->state().arm_position,
                       simulation->state().arm_velocity});
    samples.push_back(run_timed_phase(*simulation, workload, "first-step", 1));
    run_timed_phase(*simulation, workload, "warm-up", 1);
    samples.push_back(run_timed_phase(*simulation, workload, "steady-state",
                                      workload.steady_steps));
    return samples;
}

void print_header() {
        std::cout << "# benchmark=bench_cpu_sim\n"
              << "# clock=steady_clock\n"
              << "# integrator=Euler\n"
              << "# dt_s=1e-8\n"
              << "# excitation=constant_waveform_1V\n"
              << "# later_stage_trigger=positive_infinity\n"
              << "# timing=wall_and_optional_cpu_phase_sink\n"
              << "workload,thermal,phase,iterations,wall_ms,wall_ms_per_step,"
                 "derivative_count,orchestration_ms,mutual_ms,assembly_ms,solve_ms,"
                 "thermal_ms,phase_sum_ms,phase_gap_ms,finite,position_m,velocity_m_per_s\n";
}

void print_sample(const Sample& sample) {
    const double per_step = sample.iterations > 0
        ? sample.wall_ms / static_cast<double>(sample.iterations)
        : sample.wall_ms;
    const double mutual_ms = sample.timing.milliseconds(CpuPhase::Mutual);
    const double assembly_ms = sample.timing.milliseconds(CpuPhase::Assembly);
    const double solve_ms = sample.timing.milliseconds(CpuPhase::Solve);
    const double thermal_ms = sample.timing.milliseconds(CpuPhase::Thermal);
    const double orchestration_ms =
        sample.timing.milliseconds(CpuPhase::Orchestration);
    const double phase_sum = mutual_ms + assembly_ms + solve_ms + thermal_ms;
    const double phase_gap = orchestration_ms - phase_sum;
    std::cout << sample.workload << ',' << (sample.thermal ? "on" : "off") << ','
              << sample.phase << ',' << sample.iterations << ','
              << sample.wall_ms << ',' << per_step << ','
              << sample.timing.derivative_count << ',' << orchestration_ms << ','
              << mutual_ms << ',' << assembly_ms << ',' << solve_ms << ','
              << thermal_ms << ',' << phase_sum << ',' << phase_gap << ','
              << (sample.finite ? "true" : "false") << ','
              << sample.position << ',' << sample.velocity << '\n';
}

int parse_repeats(int argc, char** argv) {
    if (argc < 2) return 3;
    char* end = nullptr;
    const long value = std::strtol(argv[1], &end, 10);
    if (end == argv[1] || *end != '\0' || value <= 0 || value > 100)
        throw std::invalid_argument("repeats must be an integer in [1, 100]");
    return static_cast<int>(value);
}

} // namespace

int main(int argc, char** argv) {
    try {
        const int repeats = parse_repeats(argc, argv);
        const std::vector<Workload> workloads = {
            {"single-16", 1, 16, 1, false, 4},
            {"multi-32", 2, 16, 2, false, 4},
            {"multi-128", 8, 16, 8, false, 4},
            {"single-16-thermal", 1, 16, 1, true, 4},
            {"multi-32-thermal", 2, 16, 2, true, 4},
            {"multi-128-thermal", 8, 16, 8, true, 4},
        };

        std::cout << std::fixed << std::setprecision(6);
        print_header();
        for (const auto& workload : workloads) {
            for (int repeat = 0; repeat < repeats; ++repeat) {
                for (const auto& sample : benchmark(workload))
                    print_sample(sample);
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "bench_cpu_sim: " << error.what() << '\n';
        return 2;
    }
}
