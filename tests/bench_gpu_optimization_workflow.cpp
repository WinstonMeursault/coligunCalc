/**
 * @file bench_gpu_optimization_workflow.cpp
 * @brief Versioned real-device CPU/GPU optimization workflow evidence.
 */

#include "coilgun/optimization/coilgun_problem.hpp"
#include "coilgun/optimization/cuda_batch_evaluator.hpp"
#include "coilgun/optimization/genetic_optimizer.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/simulation/cuda/gpu_execution_report.hpp"

#include <cuda_runtime_api.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace coilgun::optimization;
using coilgun::components::Armature;
using coilgun::physics::ALUMINUM;
using coilgun::physics::COPPER;
using coilgun::simulation::OptimizationLevel;
using coilgun::simulation::cuda::BackendMode;
using coilgun::simulation::cuda::GpuBackend;

constexpr std::size_t kWarmups = 2;
constexpr std::size_t kRepetitions = 5;
constexpr std::uint64_t kSeed = 20260917;

struct BatchMeasurement {
    std::size_t size = 0;
    std::vector<double> cpu_ms;
    std::vector<double> gpu_ms;
    double cpu_median_ms = 0.0;
    double gpu_median_ms = 0.0;
    double max_abs_delta = 0.0;
    double max_relative_delta = 0.0;
    std::size_t representative_rows = 0;
    bool order_ok = false;
    bool numerical_ok = false;
    CudaExecutionSnapshot snapshot;
    std::uint64_t fallback_events = 0;
};

struct Hardware {
    std::string name;
    std::string driver;
    std::string driver_api;
    std::string runtime;
    std::string toolkit;
    std::string compute;
};

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char character : value) {
        switch (character) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default: out << character; break;
        }
    }
    return out.str();
}

std::string quote(const std::string& value) { return "\"" + json_escape(value) + "\""; }

std::string number(double value) {
    if (!std::isfinite(value)) throw std::runtime_error("benchmark generated a non-finite number");
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

double median(std::vector<double> values) {
    if (values.empty()) throw std::runtime_error("median requires samples");
    std::sort(values.begin(), values.end());
    const auto middle = values.size() / 2;
    if (values.size() % 2 != 0) return values[middle];
    return (values[middle - 1] + values[middle]) / 2.0;
}

CoilgunOptimizationProblem::Config workload_config(OptimizationLevel level = OptimizationLevel::Full) {
    CoilgunOptimizationProblem::Config config;
    config.coils.emplace_back(0.010, 0.030, 0.050, 150,
                              COPPER.resistivity_ref, 1e-6, 0.7, 0.015);
    config.armature = Armature(0.005, 0.025, 0.080,
                               ALUMINUM.resistivity_ref, ALUMINUM.density,
                               0.0, 0.120, 2, 2, 0.000);
    config.excitations = {{500.0, 500e-6, true}};
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    config.dt = 1e-6;
    config.enable_thermal = false;
    config.optimization_level = level;
    config.termination.max_steps = 64;
    config.termination.enable_velocity_check = false;
    config.termination.enable_bound_check = false;
    config.objective_id = "muzzle_velocity";
    config.constraints.push_back({"velocity_floor", CoilgunMetric::TerminalVelocity,
        ConstraintDefinition{"velocity_floor", ConstraintKind::Hard,
            ConstraintRelation::GreaterEqual, 0.0, 0.0, 0.001}});
    return config;
}

VariableSchema workload_schema() {
    return VariableSchema({VariableSpec::continuous("voltage", 450.0, 550.0)});
}

std::vector<CandidateVariables> candidates(std::size_t count) {
    std::vector<CandidateVariables> values;
    values.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        // A monotone rational stride guarantees 128 distinct values while
        // remaining inside the schema's 450..550 V bounds.
        const double voltage = 455.0 + static_cast<double>(index) * (95.0 / 128.0);
        values.emplace_back(std::vector<double>{voltage});
    }
    return values;
}

GpuBackend direct_backend() {
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    backend.max_batch_sims = 256;
    backend.enable_profiling = true;
    return backend;
}

void require_success(const std::vector<EvaluationResult>& results, const char* path) {
    for (std::size_t index = 0; index < results.size(); ++index) {
        if (results[index].status != EvaluationStatus::Success)
            throw std::runtime_error(std::string(path) + " row " + std::to_string(index) + " failed");
        if (results[index].objectives.size() != 1 ||
            !std::isfinite(results[index].objectives.front().value))
            throw std::runtime_error(std::string(path) + " row has malformed objective");
    }
}

void compare_results(const std::vector<EvaluationResult>& cpu,
                     const std::vector<EvaluationResult>& gpu,
                     BatchMeasurement& measurement) {
    if (cpu.size() != gpu.size()) throw std::runtime_error("CPU/GPU result count mismatch");
    measurement.order_ok = true;
    measurement.numerical_ok = true;
    const std::vector<std::size_t> rows = {0, measurement.size / 2, measurement.size - 1};
    measurement.representative_rows = std::min(rows.size(), measurement.size);
    const auto metric = [](const EvaluationResult& result, const char* key) {
        const auto found = result.metadata.find(key);
        if (found == result.metadata.end()) throw std::runtime_error(std::string("missing metric ") + key);
        return std::stod(found->second);
    };
    const auto close = [](double actual, double expected) {
        return std::abs(actual - expected) <=
            1e-9 + 1e-4 * std::max(std::abs(actual), std::abs(expected));
    };
    for (const auto row : rows) {
        const double actual = gpu[row].objectives.front().value;
        const double expected = cpu[row].objectives.front().value;
        const double absolute = std::abs(actual - expected);
        const double relative = absolute / std::max({std::abs(actual), std::abs(expected), 1e-300});
        measurement.max_abs_delta = std::max(measurement.max_abs_delta, absolute);
        measurement.max_relative_delta = std::max(measurement.max_relative_delta, relative);
        const double tolerance = 1e-9 + 1e-4 * std::max(std::abs(actual), std::abs(expected));
        if (gpu[row].objectives.front().id != cpu[row].objectives.front().id || absolute > tolerance)
            measurement.numerical_ok = false;
        for (const char* key : {"peak_voltage", "maximum_temperature"}) {
            const double gpu_metric = metric(gpu[row], key);
            const double cpu_metric = metric(cpu[row], key);
            const double metric_abs = std::abs(gpu_metric - cpu_metric);
            const double metric_tol = 1e-9 + 1e-4 * std::max(std::abs(gpu_metric), std::abs(cpu_metric));
            measurement.max_abs_delta = std::max(measurement.max_abs_delta, metric_abs);
            measurement.max_relative_delta = std::max(measurement.max_relative_delta,
                metric_abs / std::max({std::abs(gpu_metric), std::abs(cpu_metric), 1e-300}));
            if (metric_abs > metric_tol) measurement.numerical_ok = false;
        }
    }
    for (std::size_t row = 0; row < cpu.size(); ++row) {
        const double actual = gpu[row].objectives.front().value;
        const double expected = cpu[row].objectives.front().value;
        const double tolerance = 1e-9 + 1e-4 * std::max(std::abs(actual), std::abs(expected));
        if (gpu[row].objectives.front().id != cpu[row].objectives.front().id ||
            std::abs(actual - expected) > tolerance) {
            measurement.order_ok = false;
            measurement.numerical_ok = false;
        }
    }
    for (std::size_t row = 0; row < cpu.size(); ++row) {
        const auto& cpu_constraints = cpu[row].constraints;
        const auto& gpu_constraints = gpu[row].constraints;
        if (cpu_constraints.size() != gpu_constraints.size()) {
            measurement.numerical_ok = false;
            continue;
        }
        for (std::size_t constraint = 0; constraint < cpu_constraints.size(); ++constraint) {
            const auto& expected = cpu_constraints[constraint];
            const auto& actual = gpu_constraints[constraint];
            if (actual.id != expected.id || actual.kind != expected.kind ||
                actual.relation != expected.relation || actual.satisfied != expected.satisfied ||
                actual.priority != expected.priority ||
                !close(actual.value, expected.value) ||
                !close(actual.lower_bound, expected.lower_bound) ||
                !close(actual.upper_bound, expected.upper_bound) ||
                !close(actual.violation, expected.violation) ||
                !close(actual.normalized_violation, expected.normalized_violation))
                measurement.numerical_ok = false;
        }
    }
    if (!measurement.order_ok || !measurement.numerical_ok) {
        throw std::runtime_error("CPU/GPU ordering or Full numerical tolerance failed; max_abs=" +
                                 number(measurement.max_abs_delta) + "; max_relative=" +
                                 number(measurement.max_relative_delta));
    }
}

Hardware hardware_info() {
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count <= 0)
        throw std::runtime_error("no CUDA device available");
    int device = 0;
    if (cudaGetDevice(&device) != cudaSuccess) device = 0;
    cudaDeviceProp properties{};
    if (cudaGetDeviceProperties(&properties, device) != cudaSuccess)
        throw std::runtime_error("cudaGetDeviceProperties failed");
    int driver = 0;
    int runtime = 0;
    cudaDriverGetVersion(&driver);
    cudaRuntimeGetVersion(&runtime);
    Hardware result;
    result.name = properties.name;
    result.driver_api = std::to_string(driver / 1000) + "." + std::to_string((driver % 1000) / 10);
    if (FILE* pipe = popen("nvidia-smi --query-gpu=driver_version --format=csv,noheader,nounits", "r")) {
        char buffer[128]{};
        if (std::fgets(buffer, sizeof(buffer), pipe)) result.driver = buffer;
        pclose(pipe);
    }
    while (!result.driver.empty() &&
           (result.driver.back() == '\n' || result.driver.back() == '\r' || result.driver.back() == ' '))
        result.driver.pop_back();
    if (result.driver.empty()) result.driver = result.driver_api;
    result.runtime = std::to_string(runtime / 1000) + "." + std::to_string((runtime % 1000) / 10);
    result.toolkit = std::to_string(CUDART_VERSION / 1000) + "." +
        std::to_string((CUDART_VERSION % 1000) / 10);
    result.compute = std::to_string(properties.major) + "." + std::to_string(properties.minor);
    return result;
}

BatchMeasurement measure_batch(std::size_t size) {
    const auto schema = workload_schema();
    CoilgunOptimizationProblem cpu_problem(schema, workload_config());
    CoilgunOptimizationProblem gpu_problem(schema, workload_config());
    CudaBatchEvaluator gpu_evaluator(gpu_problem, direct_backend());
    const auto batch = candidates(size);
    const EvaluationContext context{kSeed, false};
    using clock = std::chrono::steady_clock;
    for (std::size_t warmup = 0; warmup < kWarmups; ++warmup) {
        (void)cpu_problem.evaluate_batch(batch, context);
        (void)gpu_evaluator.evaluate_batch(batch, context);
        const auto warmup_report = gpu_evaluator.execution_snapshot().report;
        if (!warmup_report.gpu_executed || warmup_report.backend == BackendMode::Fallback)
            throw std::runtime_error("GPU warm-up did not execute on a non-fallback backend");
    }

    BatchMeasurement measurement;
    measurement.size = size;
    measurement.cpu_ms.reserve(kRepetitions);
    measurement.gpu_ms.reserve(kRepetitions);
    std::vector<EvaluationResult> first_cpu;
    std::vector<EvaluationResult> first_gpu;
    for (std::size_t repetition = 0; repetition < kRepetitions; ++repetition) {
        const auto cpu_start = clock::now();
        auto cpu_results = cpu_problem.evaluate_batch(batch, context);
        const auto cpu_stop = clock::now();
        const auto gpu_start = clock::now();
        auto gpu_results = gpu_evaluator.evaluate_batch(batch, context);
        const auto gpu_stop = clock::now();
        require_success(cpu_results, "CPU Full");
        require_success(gpu_results, "GPU Full");
        const CudaExecutionSnapshot snapshot = gpu_evaluator.execution_snapshot();
        if (!snapshot.report.gpu_executed || snapshot.report.backend == BackendMode::Fallback ||
            snapshot.report.precision != coilgun::simulation::cuda::PrecisionMode::Full ||
            snapshot.report.requested_precision != coilgun::simulation::cuda::PrecisionMode::Full)
            throw std::runtime_error("GPU measurement did not prove Full non-fallback execution");
        if (snapshot.report.fallback_count != 0)
            throw std::runtime_error("GPU measurement reported fallback events");
        measurement.cpu_ms.push_back(std::chrono::duration<double, std::milli>(cpu_stop - cpu_start).count());
        measurement.gpu_ms.push_back(std::chrono::duration<double, std::milli>(gpu_stop - gpu_start).count());
        if (repetition == 0) {
            first_cpu = cpu_results;
            first_gpu = gpu_results;
            measurement.snapshot = snapshot;
        }
    }
    compare_results(first_cpu, first_gpu, measurement);
    const auto stats = gpu_evaluator.statistics_snapshot();
    if (!stats || stats->gpu_fallbacks != 0 || stats->gpu_failed_batches != 0)
        throw std::runtime_error("GPU evaluator statistics reported fallback/failure");
    measurement.fallback_events = stats ? stats->gpu_fallbacks : 0;
    measurement.cpu_median_ms = median(measurement.cpu_ms);
    measurement.gpu_median_ms = median(measurement.gpu_ms);
    return measurement;
}

OptimizationResult run_cpu_optimization(const VariableSchema& schema) {
    CoilgunOptimizationProblem problem(schema, workload_config());
    OptimizationConfig config;
    config.population_size = 8;
    config.max_generations = 4;
    config.elite_count = 1;
    config.crossover_rate = 0.8;
    config.mutation_rate = 0.2;
    config.random_seed = kSeed;
    return GeneticOptimizer(schema, problem, config).run();
}

OptimizationResult run_gpu_optimization(const VariableSchema& schema) {
    CoilgunOptimizationProblem problem(schema, workload_config());
    CudaBatchEvaluator evaluator(problem, direct_backend());
    OptimizationConfig config;
    config.population_size = 8;
    config.max_generations = 4;
    config.elite_count = 1;
    config.crossover_rate = 0.8;
    config.mutation_rate = 0.2;
    config.random_seed = kSeed;
    auto result = GeneticOptimizer(schema, evaluator, config).run();
    if (result.statistics.gpu_batches == 0 || result.statistics.gpu_fallbacks != 0 ||
        result.statistics.gpu_failed_batches != 0 ||
        result.statistics.gpu_executed_evaluations != result.statistics.gpu_successful_evaluations)
        throw std::runtime_error("GPU optimizer statistics failed execution/fallback gate");
    return result;
}

const Candidate& best_candidate(const OptimizationResult& result) {
    const auto found = result.best_by_objective.find("muzzle_velocity");
    if (found == result.best_by_objective.end()) throw std::runtime_error("optimizer returned no best candidate");
    return found->second;
}

void write_statistics(std::ostream& out, const OptimizationStatistics& stats) {
    out << "{\"seed\":" << stats.seed
        << ",\"evaluations\":" << stats.evaluations
        << ",\"successful_evaluations\":" << stats.successful_evaluations
        << ",\"failed_evaluations\":" << stats.failed_evaluations
        << ",\"cache_hits\":" << stats.cache_hits
        << ",\"gpu_fallbacks\":" << stats.gpu_fallbacks
        << ",\"gpu_requested_evaluations\":" << stats.gpu_requested_evaluations
        << ",\"gpu_executed_evaluations\":" << stats.gpu_executed_evaluations
        << ",\"gpu_successful_evaluations\":" << stats.gpu_successful_evaluations
        << ",\"gpu_failed_evaluations\":" << stats.gpu_failed_evaluations
        << ",\"cpu_fallback_evaluations\":" << stats.cpu_fallback_evaluations
        << ",\"gpu_batches\":" << stats.gpu_batches
        << ",\"gpu_failed_batches\":" << stats.gpu_failed_batches
        << ",\"gpu_transfer_seconds\":" << number(stats.gpu_transfer_seconds)
        << ",\"gpu_kernel_seconds\":" << number(stats.gpu_kernel_seconds)
        << ",\"gpu_elapsed_seconds\":" << number(stats.gpu_elapsed_seconds)
        << ",\"generations\":" << stats.generations
        << ",\"elapsed_seconds\":" << number(stats.elapsed_seconds) << "}";
}

} // namespace

int main(int argc, char** argv) {
    const std::string output = argc > 1 ? argv[1] : "optimization-gpu-workflow.json";
    try {
        const auto hardware = hardware_info();
        std::vector<BatchMeasurement> measurements;
        for (const std::size_t size : {std::size_t{1}, std::size_t{8}, std::size_t{32}, std::size_t{128}})
            measurements.push_back(measure_batch(size));
        const auto schema = workload_schema();
        const auto cpu = run_cpu_optimization(schema);
        const auto gpu = run_gpu_optimization(schema);
        const auto cpu_best = best_candidate(cpu);
        const auto gpu_best = best_candidate(gpu);
        auto cpu_reference_config = workload_config(OptimizationLevel::Reference);
        auto gpu_reference_config = cpu_reference_config;
        CoilgunOptimizationProblem cpu_reference_problem(schema, std::move(cpu_reference_config));
        CoilgunOptimizationProblem gpu_reference_problem(schema, std::move(gpu_reference_config));
        const auto cpu_reference_eval = cpu_reference_problem.evaluate(cpu_best.variables);
        const auto gpu_reference_eval = gpu_reference_problem.evaluate(gpu_best.variables);
        if (cpu_reference_eval.status != EvaluationStatus::Success ||
            gpu_reference_eval.status != EvaluationStatus::Success)
            throw std::runtime_error("Reference recheck failed");
        const double cpu_reference_error = std::abs(
            cpu_best.objectives.front().value - cpu_reference_eval.objectives.front().value);
        const double gpu_reference_error = std::abs(
            gpu_best.objectives.front().value - gpu_reference_eval.objectives.front().value);
        const double cpu_reference_tolerance = 5e-8 +
            1e-6 * std::abs(cpu_reference_eval.objectives.front().value);
        const double gpu_reference_tolerance = 5e-8 +
            1e-6 * std::abs(gpu_reference_eval.objectives.front().value);
        const bool cpu_reference_valid = cpu_reference_error <= cpu_reference_tolerance;
        const bool gpu_reference_valid = gpu_reference_error <= gpu_reference_tolerance;
        if (!cpu_reference_valid || !gpu_reference_valid)
            throw std::runtime_error("Reference recheck exceeded existing tolerance");
        // Write the document in one pass; the writer's implementation is kept in main to
        // ensure report output is never partially accepted after an exception.
        std::ofstream out(output);
        if (!out) throw std::runtime_error("cannot open benchmark output: " + output);
        out << std::setprecision(17)
            << "{\"schema_version\":1,\"source_revision\":" << quote(OPTIMIZATION_BENCH_SOURCE_REVISION)
            << ",\"worktree_state\":" << quote(OPTIMIZATION_BENCH_WORKTREE_STATE)
            << ",\"toolchain\":{\"compiler\":"
            << quote(OPTIMIZATION_BENCH_COMPILER)
            << ",\"cmake\":" << quote(OPTIMIZATION_BENCH_CMAKE_VERSION)
            << ",\"build_type\":" << quote(OPTIMIZATION_BENCH_BUILD_TYPE) << "}"
            << ",\"workload_provenance\":{\"name\":\"optimization-gpu-workflow\",\"configuration\":\"fixed-geometry-euler-full\",\"command\":\"bench_gpu_optimization_workflow\"}"
            << ",\"hardware\":{\"gpu_name\":" << quote(hardware.name)
            << ",\"driver\":" << quote(hardware.driver)
            << ",\"driver_api_version\":" << quote(hardware.driver_api)
            << ",\"compute_capability\":" << quote(hardware.compute)
            << ",\"cuda_runtime\":" << quote(hardware.runtime)
            << ",\"cuda_toolkit\":" << quote(hardware.toolkit) << "}"
            << ",\"workload\":{\"batch_sizes\":[1,8,32,128],\"warmups\":2,\"repetitions\":5"
            << ",\"precision\":\"full\",\"optimization_level\":\"full\",\"dt\":1e-6"
            << ",\"max_steps\":64,\"seed\":" << kSeed
            << ",\"termination\":\"max_steps=64;velocity_check=false;bound_check=false\",\"peak_current_excluded\":true}"
            << ",\"batches\":[";
        for (std::size_t index = 0; index < measurements.size(); ++index) {
            if (index != 0) out << ',';
            const auto& measurement = measurements[index];
            const auto& report = measurement.snapshot.report;
            out << "{\"batch_size\":" << measurement.size << ",\"cpu_samples_ms\":[";
            for (std::size_t i = 0; i < measurement.cpu_ms.size(); ++i) {
                if (i != 0) out << ','; out << number(measurement.cpu_ms[i]);
            }
            out << "],\"gpu_samples_ms\":[";
            for (std::size_t i = 0; i < measurement.gpu_ms.size(); ++i) {
                if (i != 0) out << ','; out << number(measurement.gpu_ms[i]);
            }
            out << "],\"cpu_median_ms\":" << number(measurement.cpu_median_ms)
                << ",\"gpu_median_ms\":" << number(measurement.gpu_median_ms)
                << ",\"cpu_per_candidate_ms\":" << number(measurement.cpu_median_ms / measurement.size)
                << ",\"gpu_per_candidate_ms\":" << number(measurement.gpu_median_ms / measurement.size)
                << ",\"cpu_candidates_per_second\":" << number(1000.0 * measurement.size / measurement.cpu_median_ms)
                << ",\"gpu_candidates_per_second\":" << number(1000.0 * measurement.size / measurement.gpu_median_ms)
                << ",\"speedup\":" << number(measurement.cpu_median_ms / measurement.gpu_median_ms)
                << ",\"cpu_rows\":" << measurement.size << ",\"gpu_rows\":" << measurement.size
                << ",\"representative_rows\":" << measurement.representative_rows
                << ",\"order_ok\":" << (measurement.order_ok ? "true" : "false")
                << ",\"numerical_ok\":" << (measurement.numerical_ok ? "true" : "false")
                << ",\"max_abs_delta\":" << number(measurement.max_abs_delta)
                << ",\"max_relative_delta\":" << number(measurement.max_relative_delta)
                << ",\"gpu_executed\":" << (report.gpu_executed ? "true" : "false")
                << ",\"backend\":" << quote(coilgun::simulation::cuda::to_string(report.backend))
                << ",\"requested_backend\":" << quote(coilgun::simulation::cuda::to_string(report.requested_backend))
                << ",\"solver\":" << quote(coilgun::simulation::cuda::to_string(report.solver))
                << ",\"precision\":" << quote(coilgun::simulation::cuda::to_string(report.precision))
                << ",\"requested_precision\":" << quote(coilgun::simulation::cuda::to_string(report.requested_precision))
                << ",\"fallback_events\":" << measurement.fallback_events
                << ",\"gpu_report_time_ms\":" << number(report.gpu_time_ms)
                << ",\"transfer_time_ms\":" << number(report.transfer_time_ms) << "}";
        }
        out << "],\"optimization\":{\"same_seed\":true,\"population_size\":8,\"max_generations\":4"
            << ",\"cpu\":{";
        out << "\"termination\":" << quote(to_string(cpu.termination.reason))
            << ",\"feasible\":" << (is_feasible(cpu_best.constraints) ? "true" : "false")
            << ",\"best_voltage\":" << number(cpu_best.variables.values.front())
            << ",\"best_objective\":" << number(cpu_best.objectives.front().value)
            << ",\"reference_objective\":" << number(cpu_reference_eval.objectives.front().value)
            << ",\"reference_error\":" << number(cpu_reference_error)
            << ",\"reference_tolerance\":" << number(cpu_reference_tolerance)
            << ",\"reference_valid\":" << (cpu_reference_valid ? "true" : "false")
            << ",\"statistics\":";
        write_statistics(out, cpu.statistics);
        out << "},\"gpu\":{";
        out << "\"termination\":" << quote(to_string(gpu.termination.reason))
            << ",\"feasible\":" << (is_feasible(gpu_best.constraints) ? "true" : "false")
            << ",\"best_voltage\":" << number(gpu_best.variables.values.front())
            << ",\"best_objective\":" << number(gpu_best.objectives.front().value)
            << ",\"reference_objective\":" << number(gpu_reference_eval.objectives.front().value)
            << ",\"reference_error\":" << number(gpu_reference_error)
            << ",\"reference_tolerance\":" << number(gpu_reference_tolerance)
            << ",\"reference_valid\":" << (gpu_reference_valid ? "true" : "false")
            << ",\"statistics\":";
        write_statistics(out, gpu.statistics);
        out << "}},\"decision\":{";
        // Compare per-candidate medians, not total batch latency.
        out << "\"batch32_gpu_faster\":" << (measurements[2].gpu_median_ms / 32.0 < measurements[2].cpu_median_ms / 32.0 ? "true" : "false")
            << ",\"batch128_gpu_faster\":" << (measurements[3].gpu_median_ms / 128.0 < measurements[3].cpu_median_ms / 128.0 ? "true" : "false")
            << ",\"success_criterion\":\"gpu per-candidate throughput exceeds CPU at batch 32 and 128\"}}\n";
        std::cout << "wrote " << output << '\n';
    } catch (const std::exception& error) {
        std::cerr << "benchmark failed: " << error.what() << '\n';
        return 2;
    }
    return 0;
}
