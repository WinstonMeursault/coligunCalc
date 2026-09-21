#pragma once

#include "coilgun/optimization/coilgun_problem.hpp"
#include "coilgun/optimization/evaluator.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include "coilgun/simulation/cuda/gpu_execution_report.hpp"

#include <cstddef>
#include <cstdint>
#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <vector>

namespace coilgun::optimization {

/** Snapshot of the most recent CUDA batch request and its host timing. */
struct CudaExecutionSnapshot {
    simulation::cuda::ExecutionReport report;
    double host_time_ms = 0.0;
};

/** Explicit policy for failures in the CUDA execution path.
 *
 * Strict is the default: a device or row failure remains Failed. The other
 * policies opt into CPU repair for, respectively, one failed row or the
 * whole CUDA-eligible batch. Local Invalid rows are never repaired.
 */
enum class CudaFallbackPolicy { Strict, PerCandidateCpu, WholeBatchCpu };

struct CudaFallbackOptions {
    CudaFallbackPolicy policy = CudaFallbackPolicy::Strict;
};

/** One ordered row used by the evaluator's private test-access envelope. */
struct CudaExecutionRow {
    /** Zero-based position in the compacted CUDA-eligible input batch. */
    std::size_t index = 0;
    EvaluationResult result;
};

/** Deterministic result envelope for one CUDA batch invocation's test access.
 * `rows` must contain one row for every submitted candidate, in submitted
 * order; violations are protocol errors and are never CPU-repaired.
 */
struct CudaExecutionResponse {
    simulation::cuda::ExecutionReport report;
    std::vector<CudaExecutionRow> rows;
};

/**
 * Batch evaluator for fixed-geometry CoilgunOptimizationProblem instances.
 *
 * Candidate bindings may change excitation voltage/capacitance and trigger
 * values. Every valid row is submitted to one SimBatch<EulerStepper> launch;
 * locally invalid rows remain Invalid in their original positions.
 */
class CudaBatchEvaluator final : public BatchEvaluator {
public:
    explicit CudaBatchEvaluator(const CoilgunOptimizationProblem& problem);
    CudaBatchEvaluator(const CoilgunOptimizationProblem& problem,
                       simulation::cuda::GpuBackend backend);
    CudaBatchEvaluator(const CoilgunOptimizationProblem& problem,
                       simulation::cuda::GpuBackend backend,
                       CudaFallbackOptions options);

    std::vector<EvaluationResult> evaluate_batch(
        const std::vector<CandidateVariables>& candidates,
        const EvaluationContext& context = {}) override;

    [[nodiscard]] EvaluationCacheIdentity cache_identity() const override;
    [[nodiscard]] std::optional<EvaluationStatistics> statistics_snapshot() const override;
    [[nodiscard]] CudaExecutionSnapshot execution_snapshot() const;
    [[nodiscard]] CudaFallbackOptions fallback_options() const noexcept { return options_; }

private:
    using ExecutionFunction = std::function<CudaExecutionResponse(
        const std::vector<CandidateVariables>&, const EvaluationContext&)>;

    CudaBatchEvaluator(const CoilgunOptimizationProblem& problem,
                       simulation::cuda::GpuBackend backend,
                       CudaFallbackOptions options,
                       ExecutionFunction execution);
    friend struct CudaBatchEvaluatorTestAccess;

    // The CUDA evaluator owns an immutable copy because snapshots can outlive
    // the caller's source problem, including during CPU fallback.
    std::shared_ptr<const CoilgunOptimizationProblem> problem_;
    simulation::cuda::GpuBackend backend_;
    CudaFallbackOptions options_;
    ExecutionFunction execution_;
    std::uint64_t execution_identity_ = 0;
    mutable std::mutex mutex_;
    CudaExecutionSnapshot snapshot_;
    EvaluationStatistics statistics_;
    static std::atomic<std::uint64_t> next_execution_identity_;
};

} // namespace coilgun::optimization
