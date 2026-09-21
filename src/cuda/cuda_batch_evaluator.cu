#include "coilgun/optimization/cuda_batch_evaluator.hpp"
#include "coilgun/simulation/cuda/sim_batch.hpp"

#include "coilgun/physics/constants.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace coilgun::optimization {
namespace {

using simulation::CapacitorExcitation;
using simulation::CrowbarExcitation;
using simulation::Excitation;
using simulation::TriggerConfig;
using simulation::cuda::BackendMode;
using simulation::cuda::GpuBackend;
using simulation::cuda::SimBatch;

GpuBackend default_backend() {
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    return backend;
}

std::string binding_name(CoilgunParameter parameter) {
    switch (parameter) {
    case CoilgunParameter::CoilInnerRadius: return "CoilInnerRadius";
    case CoilgunParameter::CoilOuterRadius: return "CoilOuterRadius";
    case CoilgunParameter::CoilLength: return "CoilLength";
    case CoilgunParameter::CoilTurns: return "CoilTurns";
    case CoilgunParameter::CoilPosition: return "CoilPosition";
    case CoilgunParameter::ExcitationVoltage: return "ExcitationVoltage";
    case CoilgunParameter::ExcitationCapacitance: return "ExcitationCapacitance";
    case CoilgunParameter::TriggerValue: return "TriggerValue";
    case CoilgunParameter::ArmaturePosition: return "ArmaturePosition";
    case CoilgunParameter::ArmatureVelocity: return "ArmatureVelocity";
    case CoilgunParameter::ArmatureMass: return "ArmatureMass";
    }
    return "unknown";
}

bool supported_binding(CoilgunParameter parameter) {
    return parameter == CoilgunParameter::ExcitationVoltage ||
           parameter == CoilgunParameter::ExcitationCapacitance ||
           parameter == CoilgunParameter::TriggerValue;
}

std::string number(double value) {
    std::ostringstream stream;
    stream << std::setprecision(17) << value;
    return stream.str();
}

EvaluationResult invalid(const char* code, const std::string& message) {
    return EvaluationResult::invalid(code, message);
}

EvaluationResult failed(const char* code, const std::string& message) {
    return EvaluationResult::failed(code, message);
}

struct ValidRow {
    std::size_t original_index = 0;
    std::vector<CoilgunExcitationConfig> excitations;
    std::vector<TriggerConfig> triggers;
};

} // namespace

std::atomic<std::uint64_t> CudaBatchEvaluator::next_execution_identity_{0};

CudaBatchEvaluator::CudaBatchEvaluator(const CoilgunOptimizationProblem& problem)
    : CudaBatchEvaluator(problem, default_backend(), {}, {}) {}

CudaBatchEvaluator::CudaBatchEvaluator(const CoilgunOptimizationProblem& problem,
                                       GpuBackend backend)
    : CudaBatchEvaluator(problem, std::move(backend), {}, {}) {}

CudaBatchEvaluator::CudaBatchEvaluator(const CoilgunOptimizationProblem& problem,
                                       GpuBackend backend,
                                       CudaFallbackOptions options)
    : CudaBatchEvaluator(problem, std::move(backend), options, {}) {}

CudaBatchEvaluator::CudaBatchEvaluator(const CoilgunOptimizationProblem& problem,
                                       GpuBackend backend,
                                       CudaFallbackOptions options,
                                       ExecutionFunction execution)
    : problem_(std::make_shared<CoilgunOptimizationProblem>(problem.schema(), problem.config())),
      backend_(std::move(backend)), options_(options), execution_(std::move(execution)) {
    if (execution_) execution_identity_ = ++next_execution_identity_;
    const auto& config = problem_->config();
    const auto stage_count = config.coil_specs.empty() ? config.coils.size() : config.coil_specs.size();
    if (!config.armature)
        throw std::invalid_argument("CUDA batch evaluator requires an armature");
    if (stage_count == 0)
        throw std::invalid_argument("CUDA batch evaluator requires at least one coil");
    if (!config.coil_specs.empty() && !config.coils.empty() &&
        config.coil_specs.size() != config.coils.size())
        throw std::invalid_argument("CUDA batch evaluator requires matching coil geometry");
    if (config.excitations.size() != stage_count || config.triggers.size() + 1 != stage_count)
        throw std::invalid_argument("CUDA batch evaluator configuration has inconsistent stage counts");
    if (!std::isfinite(config.dt) || config.dt <= 0.0)
        throw std::invalid_argument("CUDA batch evaluator requires a finite positive dt");
    if (config.enable_thermal)
        throw std::invalid_argument("CUDA batch evaluator does not support thermal simulation");
    if (config.optimization_level != simulation::OptimizationLevel::Full)
        throw std::invalid_argument("CUDA batch evaluator requires OptimizationLevel::Full");
    for (const auto& excitation : config.excitations) {
        if (!std::isfinite(excitation.initial_voltage) || excitation.initial_voltage <= 0.0 ||
            !std::isfinite(excitation.capacitance) ||
            excitation.capacitance <= 0.0)
            throw std::invalid_argument(
                "CUDA batch evaluator excitation initial_voltage must be finite and positive; "
                "capacitance must be finite and positive");
    }
    for (const auto& trigger : config.triggers)
        simulation::validate_trigger_config(trigger);
    for (const auto& binding : config.bindings) {
        if (!supported_binding(binding.parameter))
            throw std::invalid_argument("CUDA batch evaluator does not support candidate binding: " +
                                        binding_name(binding.parameter));
        if (binding.variable_id.empty())
            throw std::invalid_argument("CUDA batch evaluator binding variable ID must not be empty");
        const auto variable = std::find_if(
            problem_->schema().variables().begin(), problem_->schema().variables().end(),
            [&](const auto& item) { return item.id == binding.variable_id; });
        if (variable == problem_->schema().variables().end())
            throw std::invalid_argument("CUDA batch evaluator binding references unknown variable: " +
                                        binding.variable_id);
        if (binding.parameter == CoilgunParameter::TriggerValue) {
            if (binding.index >= config.triggers.size())
                throw std::out_of_range("CUDA batch evaluator trigger binding index out of range");
        } else if (binding.index >= config.excitations.size()) {
            throw std::out_of_range("CUDA batch evaluator excitation binding index out of range");
        }
    }
    for (const auto& constraint : config.constraints) {
        if (constraint.metric == CoilgunMetric::MaximumTemperature)
            throw std::invalid_argument("CUDA batch evaluator does not support thermal metrics");
        if (constraint.metric == CoilgunMetric::PeakCurrent)
            throw std::invalid_argument("CUDA batch evaluator does not support PeakCurrent metrics");
    }
    switch (options_.policy) {
    case CudaFallbackPolicy::Strict:
    case CudaFallbackPolicy::PerCandidateCpu:
    case CudaFallbackPolicy::WholeBatchCpu:
        break;
    default:
        throw std::invalid_argument("CUDA fallback policy is invalid");
    }
    backend_.validate();
}

std::vector<EvaluationResult> CudaBatchEvaluator::evaluate_batch(
    const std::vector<CandidateVariables>& candidates, const EvaluationContext& context) {
    const auto start = std::chrono::steady_clock::now();
    if (candidates.empty()) return {};

    std::vector<EvaluationResult> results(candidates.size());
    std::vector<ValidRow> valid_rows;
    valid_rows.reserve(candidates.size());
    const auto& config = problem_->config();
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        const auto& candidate = candidates[index];
        if (candidate.values.size() != problem_->schema().size()) {
            results[index] = invalid("variable_count", "candidate variable count does not match schema");
            continue;
        }
        if (!std::all_of(candidate.values.begin(), candidate.values.end(),
                         [](double value) { return std::isfinite(value); })) {
            results[index] = invalid("non_finite_variable", "candidate variables must be finite");
            continue;
        }

        CandidateVariables decoded;
        try {
            decoded = problem_->schema().decode(candidate);
        } catch (const std::exception& error) {
            results[index] = invalid("variable_decode", error.what());
            continue;
        }
        ValidRow row;
        row.original_index = index;
        row.excitations = config.excitations;
        row.triggers = config.triggers;
        bool valid = true;
        std::string error_message;
        for (const auto& binding : config.bindings) {
            const auto variable = std::find_if(
                problem_->schema().variables().begin(), problem_->schema().variables().end(),
                [&](const auto& item) { return item.id == binding.variable_id; });
            if (variable == problem_->schema().variables().end()) {
                valid = false;
                error_message = "binding references unknown variable: " + binding.variable_id;
                break;
            }
            const auto variable_index = static_cast<std::size_t>(
                std::distance(problem_->schema().variables().begin(), variable));
            const double value = decoded.values[variable_index];
            if (binding.parameter == CoilgunParameter::ExcitationVoltage)
                row.excitations[binding.index].initial_voltage = value;
            else if (binding.parameter == CoilgunParameter::ExcitationCapacitance)
                row.excitations[binding.index].capacitance = value;
            else
                row.triggers[binding.index].value = value;
        }
        for (const auto& excitation : row.excitations) {
            if (!std::isfinite(excitation.initial_voltage) || excitation.initial_voltage <= 0.0 ||
                !std::isfinite(excitation.capacitance) ||
                excitation.capacitance <= 0.0) {
                valid = false;
                error_message = "excitation initial_voltage must be finite and positive; capacitance must be finite and positive";
                break;
            }
        }
        if (valid) {
            for (const auto& trigger : row.triggers) {
                try {
                    simulation::validate_trigger_config(trigger);
                } catch (const std::exception& error) {
                    valid = false;
                    error_message = error.what();
                    break;
                }
            }
        }
        if (!valid) {
            results[index] = invalid("invalid_configuration", error_message);
            continue;
        }
        valid_rows.push_back(std::move(row));
    }

    CudaExecutionSnapshot current;
    current.report.requested_backend = backend_.backend;
    current.report.requested_precision = simulation::cuda::PrecisionMode::Full;
    current.report.requested_thermal = simulation::cuda::ThermalMode::Disabled;
    EvaluationStatistics delta;
    delta.seed = context.seed;
    delta.gpu_requested_evaluations = valid_rows.size();

    auto add_diagnostic = [](EvaluationResult& result, std::string code,
                             std::string message, DiagnosticSeverity severity) {
        result.diagnostics.push_back({std::move(code), std::move(message), severity});
    };
    auto cpu_fallback = [&](const ValidRow& row, const EvaluationResult& gpu_result,
                            const std::string& reason) {
        auto fallback = problem_->evaluate_cpu(candidates[row.original_index]);
        std::string failure_code = "unknown";
        std::string failure_message = reason;
        if (!gpu_result.diagnostics.empty()) {
            failure_code = gpu_result.diagnostics.front().code;
            failure_message = gpu_result.diagnostics.front().message;
        }
        fallback.metadata["cpu_fallback"] = "true";
        fallback.metadata["gpu_failure_code"] = failure_code;
        fallback.metadata["gpu_failure_reason"] = failure_message;
        add_diagnostic(fallback, "gpu_row_failure", failure_message, DiagnosticSeverity::Warning);
        add_diagnostic(fallback, "cpu_fallback", reason, DiagnosticSeverity::Info);
        ++delta.cpu_fallback_evaluations;
        results[row.original_index] = std::move(fallback);
    };
    auto mark_batch_failed = [&](const char* code, const std::string& message,
                                 bool permit_whole_batch_fallback) {
        ++delta.gpu_failed_batches;
        if (permit_whole_batch_fallback) ++delta.fallbacks;
        const bool fallback_enabled = permit_whole_batch_fallback &&
            options_.policy == CudaFallbackPolicy::WholeBatchCpu;
        if (fallback_enabled) {
            ++delta.gpu_fallbacks;
            for (const auto& row : valid_rows) {
                auto gpu_failure = failed(code, message);
                cpu_fallback(row, gpu_failure, "whole CUDA batch fallback");
            }
        } else {
            if (permit_whole_batch_fallback) ++delta.gpu_fallbacks;
            for (const auto& row : valid_rows)
                results[row.original_index] = failed(code, message);
        }
    };

    if (valid_rows.empty()) {
        // No CUDA invocation was attempted. Do not attribute local validation
        // time to GPU elapsed time or the latest CUDA execution snapshot.
        current.host_time_ms = 0.0;
    } else {
        ++delta.gpu_batches;
        CudaExecutionResponse response;
        try {
            if (execution_) {
                std::vector<CandidateVariables> gpu_candidates;
                gpu_candidates.reserve(valid_rows.size());
                for (const auto& row : valid_rows) gpu_candidates.push_back(candidates[row.original_index]);
                response = execution_(gpu_candidates, context);
            } else {
                std::vector<components::DrivingCoil> fixed_coils = config.coils;
                if (fixed_coils.empty()) {
                    fixed_coils.reserve(config.coil_specs.size());
                    for (const auto& spec : config.coil_specs) fixed_coils.push_back(spec.make_coil());
                }
                SimBatch<simulation::EulerStepper> batch(
                    std::move(fixed_coils), *config.armature, static_cast<int>(valid_rows.size()), config.dt,
                    backend_);
                for (std::size_t row_index = 0; row_index < valid_rows.size(); ++row_index) {
                    const auto& row = valid_rows[row_index];
                    std::vector<std::unique_ptr<Excitation>> excitations;
                    excitations.reserve(row.excitations.size());
                    for (const auto& excitation : row.excitations) {
                        if (excitation.crowbar)
                            excitations.push_back(std::make_unique<CrowbarExcitation>(
                                excitation.initial_voltage, excitation.capacitance));
                        else
                            excitations.push_back(std::make_unique<CapacitorExcitation>(
                                excitation.initial_voltage, excitation.capacitance));
                    }
                    batch.set_excitations(static_cast<int>(row_index), std::move(excitations), row.triggers);
                }
                batch.run(config.termination);
                response.report = batch.execution_report();
                response.rows.reserve(valid_rows.size());
                for (std::size_t row_index = 0; row_index < valid_rows.size(); ++row_index) {
                    const auto& row = valid_rows[row_index];
                    response.rows.push_back({row_index,
                        problem_->result_from_simulation(batch.result(static_cast<int>(row_index)),
                                                        *config.armature, row.excitations)});
                }
            }
            current.report = response.report;
            current.host_time_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start).count();
            const bool gpu_resolved = current.report.gpu_executed &&
                current.report.backend != BackendMode::Fallback;
            if (!gpu_resolved) {
                std::string message = "CUDA batch resolved to a non-GPU backend";
                if (!current.report.fallback_reason.empty()) message += ": " + current.report.fallback_reason;
                mark_batch_failed("gpu_backend_fallback", message, true);
            } else {
                bool protocol_ok = response.rows.size() == valid_rows.size();
                if (!protocol_ok) {
                    mark_batch_failed("gpu_protocol_error", "CUDA batch returned an unexpected result count", false);
                } else {
                    for (std::size_t row_index = 0; row_index < response.rows.size(); ++row_index) {
                        if (response.rows[row_index].index != row_index) {
                            mark_batch_failed("gpu_protocol_error", "CUDA batch returned results out of order", false);
                            protocol_ok = false;
                            break;
                        }
                    }
                }
                if (protocol_ok) {
                    for (std::size_t row_index = 0; row_index < valid_rows.size(); ++row_index) {
                        const auto& row = valid_rows[row_index];
                        auto result = response.rows[row_index].result;
                        ++delta.gpu_executed_evaluations;
                        if (result.status == EvaluationStatus::Success) {
                            ++delta.gpu_successful_evaluations;
                            result.metadata["cuda_backend"] = simulation::cuda::to_string(current.report.backend);
                            result.metadata["cuda_solver"] = simulation::cuda::to_string(current.report.solver);
                            result.metadata["cuda_gpu_executed"] = "true";
                            result.metadata["cuda_transfer_time_ms"] = number(current.report.transfer_time_ms);
                            result.metadata["cuda_gpu_time_ms"] = number(current.report.gpu_time_ms);
                            result.metadata["cuda_host_time_ms"] = number(current.host_time_ms);
                            results[row.original_index] = std::move(result);
                        } else {
                            ++delta.gpu_failed_evaluations;
                            if (options_.policy == CudaFallbackPolicy::PerCandidateCpu)
                                cpu_fallback(row, result, "per-candidate CUDA row fallback");
                            else {
                                result.status = EvaluationStatus::Failed;
                                if (result.diagnostics.empty())
                                    add_diagnostic(result, "gpu_row_failure",
                                                   "CUDA result row failed before CPU fallback",
                                                   DiagnosticSeverity::Error);
                                results[row.original_index] = std::move(result);
                            }
                        }
                    }
                }
            }
        } catch (const std::exception& error) {
            current.host_time_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start).count();
            mark_batch_failed("gpu_execution_failed", error.what(), true);
        } catch (...) {
            current.host_time_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start).count();
            mark_batch_failed("gpu_execution_failed", "unknown CUDA exception", true);
        }
    }

    auto finite_nonnegative = [](double value) {
        return std::isfinite(value) && value >= 0.0 ? value : 0.0;
    };
    delta.gpu_transfer_seconds = finite_nonnegative(current.report.transfer_time_ms) / 1000.0;
    delta.gpu_kernel_seconds = finite_nonnegative(current.report.gpu_time_ms) / 1000.0;
    delta.gpu_elapsed_seconds = finite_nonnegative(current.host_time_ms) / 1000.0;
    {
        std::lock_guard lock(mutex_);
        snapshot_ = current;
        statistics_.seed = context.seed;
        // Retain the lifetime host timing for the compatibility snapshot;
        // run-local wall time remains optimizer-owned and GPU elapsed time is
        // carried separately through the collector.
        statistics_.elapsed_seconds += delta.gpu_elapsed_seconds;
        statistics_.fallbacks += delta.fallbacks;
        statistics_.gpu_requested_evaluations += delta.gpu_requested_evaluations;
        statistics_.gpu_executed_evaluations += delta.gpu_executed_evaluations;
        statistics_.gpu_successful_evaluations += delta.gpu_successful_evaluations;
        statistics_.gpu_failed_evaluations += delta.gpu_failed_evaluations;
        statistics_.cpu_fallback_evaluations += delta.cpu_fallback_evaluations;
        statistics_.gpu_batches += delta.gpu_batches;
        statistics_.gpu_failed_batches += delta.gpu_failed_batches;
        statistics_.gpu_fallbacks += delta.gpu_fallbacks;
        statistics_.gpu_transfer_seconds += delta.gpu_transfer_seconds;
        statistics_.gpu_kernel_seconds += delta.gpu_kernel_seconds;
        statistics_.gpu_elapsed_seconds += delta.gpu_elapsed_seconds;
    }
    if (context.statistics) context.statistics->add(delta);
    return results;
}

EvaluationCacheIdentity CudaBatchEvaluator::cache_identity() const {
    std::string namespace_id = "coilgun.optimization.cuda.batch:";
    namespace_id += problem_->cache_identity().namespace_id;
    namespace_id += ":backend=" + std::to_string(static_cast<int>(backend_.backend));
    namespace_id += ":device=" + std::to_string(backend_.device_id);
    namespace_id += ":threads=" + std::to_string(backend_.threads_per_block);
    namespace_id += ":batch=" + std::to_string(backend_.max_batch_sims);
    namespace_id += ":profiling=" + std::to_string(backend_.enable_profiling ? 1 : 0);
    namespace_id += ":persistent=" + std::to_string(backend_.use_persistent ? 1 : 0);
    namespace_id += ":fallback=" + std::to_string(static_cast<int>(options_.policy));
    namespace_id += ":execution=" + std::to_string(execution_identity_);
    return {std::move(namespace_id), "cuda-batch-result-v2"};
}

std::optional<EvaluationStatistics> CudaBatchEvaluator::statistics_snapshot() const {
    std::lock_guard lock(mutex_);
    return statistics_;
}

CudaExecutionSnapshot CudaBatchEvaluator::execution_snapshot() const {
    std::lock_guard lock(mutex_);
    return snapshot_;
}

} // namespace coilgun::optimization
