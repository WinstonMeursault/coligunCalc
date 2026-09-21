#pragma once
#include <cstdint>
#include <cmath>
#include <mutex>
namespace coilgun::optimization {
struct EvaluationStatistics {
    // evaluations counts uncached delegate evaluations; cache_hits counts results served by the cache.
    std::uint64_t evaluations = 0, successful_evaluations = 0, failed_evaluations = 0;
    std::uint64_t cache_hits = 0, fallbacks = 0, seed = 0;
    double elapsed_seconds = 0.0;
    // CUDA-owned counters. These are additive and intentionally separate from
    // generic evaluation/success/failure counts owned by optimizer wrappers.
    // gpu_kernel_seconds is ExecutionReport::gpu_time_ms converted to seconds;
    // it is the report's measured CUDA physical-pipeline time, not a pure
    // device kernel-event timer.
    std::uint64_t gpu_requested_evaluations = 0;
    std::uint64_t gpu_executed_evaluations = 0;
    std::uint64_t gpu_successful_evaluations = 0;
    std::uint64_t gpu_failed_evaluations = 0;
    std::uint64_t cpu_fallback_evaluations = 0;
    std::uint64_t gpu_batches = 0;
    std::uint64_t gpu_failed_batches = 0;
    std::uint64_t gpu_fallbacks = 0;
    double gpu_transfer_seconds = 0.0;
    double gpu_kernel_seconds = 0.0;
    double gpu_elapsed_seconds = 0.0;
};

/** Thread-safe, additive statistics sink owned by one optimization run. */
class EvaluationStatisticsCollector final {
public:
    void add(const EvaluationStatistics& delta) noexcept {
        std::lock_guard lock(mutex_);
        statistics_.evaluations += delta.evaluations;
        statistics_.successful_evaluations += delta.successful_evaluations;
        statistics_.failed_evaluations += delta.failed_evaluations;
        statistics_.cache_hits += delta.cache_hits;
        statistics_.fallbacks += delta.fallbacks;
        statistics_.elapsed_seconds += finite_nonnegative(delta.elapsed_seconds);
        statistics_.gpu_requested_evaluations += delta.gpu_requested_evaluations;
        statistics_.gpu_executed_evaluations += delta.gpu_executed_evaluations;
        statistics_.gpu_successful_evaluations += delta.gpu_successful_evaluations;
        statistics_.gpu_failed_evaluations += delta.gpu_failed_evaluations;
        statistics_.cpu_fallback_evaluations += delta.cpu_fallback_evaluations;
        statistics_.gpu_batches += delta.gpu_batches;
        statistics_.gpu_failed_batches += delta.gpu_failed_batches;
        statistics_.gpu_fallbacks += delta.gpu_fallbacks;
        statistics_.gpu_transfer_seconds += finite_nonnegative(delta.gpu_transfer_seconds);
        statistics_.gpu_kernel_seconds += finite_nonnegative(delta.gpu_kernel_seconds);
        statistics_.gpu_elapsed_seconds += finite_nonnegative(delta.gpu_elapsed_seconds);
    }
    void set_seed(std::uint64_t value) noexcept {
        std::lock_guard lock(mutex_);
        statistics_.seed = value;
    }
    void add_evaluations(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::evaluations, value); }
    void record_evaluations(std::uint64_t value = 1) noexcept { add_evaluations(value); }
    void add_successful_evaluations(std::uint64_t value = 1) noexcept {
        add_member(&EvaluationStatistics::successful_evaluations, value);
    }
    void record_successful_evaluations(std::uint64_t value = 1) noexcept {
        add_successful_evaluations(value);
    }
    void add_failed_evaluations(std::uint64_t value = 1) noexcept {
        add_member(&EvaluationStatistics::failed_evaluations, value);
    }
    void record_failed_evaluations(std::uint64_t value = 1) noexcept {
        add_failed_evaluations(value);
    }
    void add_cache_hits(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::cache_hits, value); }
    void record_cache_hit() noexcept { add_cache_hits(); }
    void record_cache_hits(std::uint64_t value = 1) noexcept { add_cache_hits(value); }
    void add_fallbacks(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::fallbacks, value); }
    void record_fallback() noexcept { add_fallbacks(); }
    void record_fallbacks(std::uint64_t value = 1) noexcept { add_fallbacks(value); }
    void add_gpu_requested_evaluations(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::gpu_requested_evaluations, value); }
    void add_gpu_executed_evaluations(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::gpu_executed_evaluations, value); }
    void add_gpu_successful_evaluations(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::gpu_successful_evaluations, value); }
    void add_gpu_failed_evaluations(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::gpu_failed_evaluations, value); }
    void add_cpu_fallback_evaluations(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::cpu_fallback_evaluations, value); }
    void add_gpu_batches(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::gpu_batches, value); }
    void add_gpu_failed_batches(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::gpu_failed_batches, value); }
    void add_gpu_fallbacks(std::uint64_t value = 1) noexcept { add_member(&EvaluationStatistics::gpu_fallbacks, value); }
    void add_elapsed_seconds(double value) noexcept {
        std::lock_guard lock(mutex_);
        statistics_.elapsed_seconds += finite_nonnegative(value);
    }
    void set_elapsed_seconds(double value) noexcept {
        std::lock_guard lock(mutex_);
        statistics_.elapsed_seconds = finite_nonnegative(value);
    }
    void add_gpu_transfer_seconds(double value) noexcept { add_duration(&EvaluationStatistics::gpu_transfer_seconds, value); }
    void add_gpu_kernel_seconds(double value) noexcept { add_duration(&EvaluationStatistics::gpu_kernel_seconds, value); }
    void add_gpu_elapsed_seconds(double value) noexcept { add_duration(&EvaluationStatistics::gpu_elapsed_seconds, value); }
    [[nodiscard]] EvaluationStatistics snapshot() const noexcept {
        std::lock_guard lock(mutex_);
        return statistics_;
    }
    [[nodiscard]] EvaluationStatistics statistics() const noexcept { return snapshot(); }
    [[nodiscard]] EvaluationStatistics statistics_snapshot() const noexcept { return snapshot(); }

private:
    static double finite_nonnegative(double value) noexcept {
        return std::isfinite(value) && value >= 0.0 ? value : 0.0;
    }

    using Counter = std::uint64_t EvaluationStatistics::*;
    void add_member(Counter member, std::uint64_t value) noexcept {
        std::lock_guard lock(mutex_);
        statistics_.*member += value;
    }
    using Duration = double EvaluationStatistics::*;
    void add_duration(Duration member, double value) noexcept {
        std::lock_guard lock(mutex_);
        if (value >= 0.0 && std::isfinite(value)) statistics_.*member += value;
    }
    mutable std::mutex mutex_;
    EvaluationStatistics statistics_;
};

using EvaluationStatisticsSink = EvaluationStatisticsCollector;
using RunStatisticsCollector = EvaluationStatisticsCollector;

}
