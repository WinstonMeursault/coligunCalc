#include "coilgun/optimization/evaluator.hpp"
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <unordered_map>
#include <utility>
namespace coilgun::optimization {
namespace {
void copy_cuda_statistics(EvaluationStatistics& destination,
                          const EvaluationStatistics& source) {
    destination.gpu_requested_evaluations = source.gpu_requested_evaluations;
    destination.gpu_executed_evaluations = source.gpu_executed_evaluations;
    destination.gpu_successful_evaluations = source.gpu_successful_evaluations;
    destination.gpu_failed_evaluations = source.gpu_failed_evaluations;
    destination.cpu_fallback_evaluations = source.cpu_fallback_evaluations;
    destination.gpu_batches = source.gpu_batches;
    destination.gpu_failed_batches = source.gpu_failed_batches;
    destination.gpu_fallbacks = source.gpu_fallbacks;
    destination.gpu_transfer_seconds = source.gpu_transfer_seconds;
    destination.gpu_kernel_seconds = source.gpu_kernel_seconds;
    destination.gpu_elapsed_seconds = source.gpu_elapsed_seconds;
}

EvaluationResult exception_result(const std::exception& error) {
    return EvaluationResult::failed("evaluation_exception", error.what());
}

EvaluationResult unknown_exception_result() {
    return EvaluationResult::failed("evaluation_exception", "unknown exception");
}

EvaluationResult normalize_result(EvaluationResult result) {
    if (result.status == EvaluationStatus::Unevaluated) {
        return EvaluationResult::failed("evaluation_unevaluated", "evaluator returned an unevaluated result");
    }
    return result;
}

std::vector<EvaluationResult> normalize_batch(std::vector<EvaluationResult> results,
                                               std::size_t candidate_count) {
    if (results.size() < candidate_count) {
        const auto missing = candidate_count - results.size();
        for (std::size_t i = 0; i < missing; ++i) {
            results.push_back(EvaluationResult::failed(
                "evaluation_batch_output", "batch evaluator returned too few results"));
        }
    }
    if (results.size() > candidate_count) results.resize(candidate_count);
    for (auto& result : results) result = normalize_result(std::move(result));
    return results;
}

std::vector<EvaluationResult> safe_batch(BatchEvaluator& evaluator,
                                         const std::vector<CandidateVariables>& candidates,
                                         const EvaluationContext& context) {
    if (candidates.empty()) return {};
    try {
        return normalize_batch(evaluator.evaluate_batch(candidates, context), candidates.size());
    } catch (const std::exception& error) {
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for (std::size_t i = 0; i < candidates.size(); ++i) results.push_back(exception_result(error));
        return results;
    } catch (...) {
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for (std::size_t i = 0; i < candidates.size(); ++i) results.push_back(unknown_exception_result());
        return results;
    }
}

std::vector<EvaluationResult> safe_batch_with_isolation(
    BatchEvaluator& evaluator, const std::vector<CandidateVariables>& candidates,
    const EvaluationContext& context) {
    if (candidates.empty()) return {};
    try {
        return normalize_batch(evaluator.evaluate_batch(candidates, context), candidates.size());
    } catch (...) {
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for (const auto& candidate : candidates) {
            const auto single = safe_batch(evaluator, {candidate}, context);
            results.push_back(single.front());
        }
        return results;
    }
}

std::vector<EvaluationResult> safe_snapshot_with_isolation(
    const BatchEvaluationSnapshot& snapshot,
    const std::vector<CandidateVariables>& candidates,
    const EvaluationContext& context) {
    if (candidates.empty()) return {};
    try {
        return normalize_batch(snapshot.evaluate(candidates, context), candidates.size());
    } catch (...) {
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for (const auto& candidate : candidates) {
            try {
                results.push_back(normalize_batch(snapshot.evaluate({candidate}, context), 1).front());
            } catch (const std::exception& error) {
                results.push_back(exception_result(error));
            } catch (...) {
                results.push_back(unknown_exception_result());
            }
        }
        return results;
    }
}

EvaluationResult one(const Evaluator& evaluator, const CandidateVariables& variables,
                     const EvaluationContext& context) {
    try {
        return normalize_result(evaluator.evaluate(variables, context));
    } catch (const std::exception& error) {
        return exception_result(error);
    } catch (...) {
        return unknown_exception_result();
    }
}

void record_status(EvaluationStatistics& statistics, const EvaluationResult& result) {
    if (result.status == EvaluationStatus::Success) {
        ++statistics.successful_evaluations;
    } else if (result.status == EvaluationStatus::Failed || result.status == EvaluationStatus::Invalid ||
               result.status == EvaluationStatus::Unevaluated) {
        ++statistics.failed_evaluations;
    }
}

void record_status(const std::shared_ptr<EvaluationStatisticsCollector>& collector,
                   const EvaluationResult& result) {
    if (!collector) return;
    if (result.status == EvaluationStatus::Success) {
        collector->add_successful_evaluations();
    } else if (result.status == EvaluationStatus::Failed || result.status == EvaluationStatus::Invalid ||
               result.status == EvaluationStatus::Unevaluated) {
        collector->add_failed_evaluations();
    }
}

class LegacyFallbackScope final {
public:
    explicit LegacyFallbackScope(const EvaluationContext& context)
        : context_(context), active_(context.fallback && context.statistics) {
        if (active_) {
            const auto nested = std::find(fallback_contexts_.rbegin(), fallback_contexts_.rend(), &context_);
            owner_ = nested == fallback_contexts_.rend();
            fallback_contexts_.push_back(&context_);
        }
    }
    ~LegacyFallbackScope() {
        if (active_) fallback_contexts_.pop_back();
    }
    [[nodiscard]] bool owner() const noexcept { return owner_; }

private:
    const EvaluationContext& context_;
    bool active_ = false;
    bool owner_ = false;
    static thread_local std::vector<const EvaluationContext*> fallback_contexts_;
};

thread_local std::vector<const EvaluationContext*> LegacyFallbackScope::fallback_contexts_;
}

BatchEvaluator::BatchEvaluator() = default;

BatchEvaluationSnapshot BatchEvaluator::evaluation_snapshot() {
    std::shared_ptr<BatchEvaluator> owner;
    try {
        owner = shared_from_this();
    } catch (const std::bad_weak_ptr&) {
        throw std::logic_error(
            "BatchEvaluator::evaluation_snapshot requires shared ownership");
    }
    return make_owned_snapshot(make_evaluation_snapshot(), std::move(owner));
}

BatchEvaluationSnapshot BatchEvaluator::make_evaluation_snapshot() {
    auto identity = cache_identity();
    return {std::move(identity), [this](const auto& candidates, const auto& context) {
                return evaluate_batch(candidates, context);
            }};
}

BatchEvaluationSnapshot BatchEvaluator::make_owned_snapshot(
    BatchEvaluationSnapshot snapshot, std::shared_ptr<BatchEvaluator> owner) {
    auto delegate = std::move(snapshot.evaluate);
    snapshot.evaluate = [owner = std::move(owner), delegate = std::move(delegate)](
                            const std::vector<CandidateVariables>& candidates,
                            const EvaluationContext& context) {
        (void)owner;
        return delegate(candidates, context);
    };
    return snapshot;
}

SerialBatchEvaluator::SerialBatchEvaluator(std::shared_ptr<const Evaluator> e) : evaluator_(std::move(e)) { if (!evaluator_) throw std::invalid_argument("evaluator must not be null"); }
std::vector<EvaluationResult> SerialBatchEvaluator::evaluate_batch(const std::vector<CandidateVariables>& c, const EvaluationContext& x) { std::vector<EvaluationResult> r; r.reserve(c.size()); for (const auto& v:c) r.push_back(one(*evaluator_,v,x)); return r; }
StatisticsBatchEvaluator::StatisticsBatchEvaluator(std::shared_ptr<BatchEvaluator> e) : evaluator_(std::move(e)) { if (!evaluator_) throw std::invalid_argument("evaluator must not be null"); }
EvaluationCacheIdentity StatisticsBatchEvaluator::cache_identity() const { return evaluator_->cache_identity(); }
EvaluationStatistics StatisticsBatchEvaluator::statistics() const {
    std::lock_guard lock(statistics_mutex_);
    return statistics_;
}
std::optional<EvaluationStatistics> StatisticsBatchEvaluator::statistics_snapshot() const {
    auto snapshot = statistics();
    if (const auto nested = evaluator_->statistics_snapshot()) {
        snapshot.fallbacks = nested->fallbacks;
        copy_cuda_statistics(snapshot, *nested);
    }
    return snapshot;
}
std::vector<EvaluationResult> StatisticsBatchEvaluator::evaluate_batch(const std::vector<CandidateVariables>& c, const EvaluationContext& x) {
    const auto start = std::chrono::steady_clock::now();
    {
        std::lock_guard lock(statistics_mutex_);
        statistics_.seed = x.seed;
        if (x.fallback) ++statistics_.fallbacks;
    }
    LegacyFallbackScope fallback_scope{x};
    auto results = safe_batch_with_isolation(*evaluator_, c, x);
    {
        std::lock_guard lock(statistics_mutex_);
        statistics_.evaluations += c.size();
        for (const auto& result : results) record_status(statistics_, result);
        statistics_.elapsed_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
    if (x.statistics) {
        x.statistics->add_evaluations(c.size());
        for (const auto& result : results) record_status(x.statistics, result);
        if (fallback_scope.owner()) x.statistics->add_fallbacks();
    }
    return results;
}
BatchEvaluationSnapshot StatisticsBatchEvaluator::make_evaluation_snapshot() {
    auto nested = evaluator_->evaluation_snapshot();
    return {nested.identity, [this, nested = std::move(nested)](
                                const std::vector<CandidateVariables>& candidates,
                                const EvaluationContext& context) {
        const auto start = std::chrono::steady_clock::now();
        {
            std::lock_guard lock(statistics_mutex_);
            statistics_.seed = context.seed;
            if (context.fallback) ++statistics_.fallbacks;
        }
        LegacyFallbackScope fallback_scope{context};
        auto results = safe_snapshot_with_isolation(nested, candidates, context);
        {
            std::lock_guard lock(statistics_mutex_);
            statistics_.evaluations += candidates.size();
            for (const auto& result : results) record_status(statistics_, result);
            statistics_.elapsed_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
        }
        if (context.statistics) {
            context.statistics->add_evaluations(candidates.size());
            for (const auto& result : results) record_status(context.statistics, result);
            if (fallback_scope.owner()) context.statistics->add_fallbacks();
        }
        return results;
    }};
}
CachedBatchEvaluator::CachedBatchEvaluator(std::shared_ptr<BatchEvaluator> e,std::shared_ptr<EvaluationCache> c):evaluator_(std::move(e)),cache_(std::move(c)){if(!evaluator_||!cache_)throw std::invalid_argument("evaluator and cache must not be null");}
EvaluationCacheIdentity CachedBatchEvaluator::cache_identity() const { return evaluator_->cache_identity(); }
EvaluationStatistics CachedBatchEvaluator::statistics() const {
    std::lock_guard lock(statistics_mutex_);
    return statistics_;
}
std::optional<EvaluationStatistics> CachedBatchEvaluator::statistics_snapshot() const {
    auto snapshot = statistics();
    if (const auto nested = evaluator_->statistics_snapshot()) {
        snapshot.fallbacks = nested->fallbacks;
        copy_cuda_statistics(snapshot, *nested);
    }
    return snapshot;
}
BatchEvaluationSnapshot CachedBatchEvaluator::make_evaluation_snapshot() {
    auto nested = evaluator_->evaluation_snapshot();
    return {nested.identity, [this, nested = std::move(nested)](
                                const std::vector<CandidateVariables>& candidates,
                                const EvaluationContext& context) {
        return evaluate_batch_with_snapshot(nested, candidates, context);
    }};
}
std::vector<EvaluationResult> CachedBatchEvaluator::evaluate_batch(
    const std::vector<CandidateVariables>& c, const EvaluationContext& x) {
    return evaluate_batch_with_snapshot(evaluator_->evaluation_snapshot(), c, x);
}
std::vector<EvaluationResult> CachedBatchEvaluator::evaluate_batch_with_snapshot(
    const BatchEvaluationSnapshot& snapshot,
    const std::vector<CandidateVariables>& c,
    const EvaluationContext& x) {
    const auto start = std::chrono::steady_clock::now();
    std::vector<EvaluationResult> results;
    results.reserve(c.size());
    {
        std::lock_guard lock(statistics_mutex_);
        statistics_.seed = x.seed;
        if (x.fallback) ++statistics_.fallbacks;
    }
    std::vector<CandidateVariables> misses;
    std::vector<std::size_t> miss_candidate_indices;
    std::vector<std::size_t> miss_result_indices;
    std::vector<std::string> miss_keys;
    std::unordered_map<std::string, std::size_t> pending;
    misses.reserve(c.size());
    miss_candidate_indices.reserve(c.size());
    miss_result_indices.reserve(c.size());
    miss_keys.reserve(c.size());
    results.resize(c.size());
    const auto& identity = snapshot.identity;
    for (std::size_t i = 0; i < c.size(); ++i) {
        const auto& variables = c[i];
        const auto key = make_cache_key(identity, variables, x);
        if (auto hit = cache_->get(key)) {
            results[i] = normalize_result(*hit);
            {
                std::lock_guard lock(statistics_mutex_);
                ++statistics_.cache_hits;
            }
            if (x.statistics) x.statistics->add_cache_hits();
            continue;
        }
        const auto [pending_it, inserted] = pending.emplace(key, misses.size());
        if (inserted) {
            misses.push_back(variables);
            miss_keys.push_back(key);
        }
        miss_candidate_indices.push_back(i);
        miss_result_indices.push_back(pending_it->second);
    }
    const auto fresh = safe_snapshot_with_isolation(snapshot, misses, x);
    for (std::size_t i = 0; i < fresh.size(); ++i) {
        const auto result = normalize_result(fresh[i]);
        cache_->put(miss_keys[i], result);
        {
            std::lock_guard lock(statistics_mutex_);
            record_status(statistics_, result);
        }
    }
    for (std::size_t i = 0; i < miss_candidate_indices.size(); ++i) {
        results[miss_candidate_indices[i]] = fresh[miss_result_indices[i]];
    }
    {
        std::lock_guard lock(statistics_mutex_);
        statistics_.evaluations += misses.size();
        statistics_.elapsed_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
    return results;
}
}
