#include "coilgun/optimization/evaluator.hpp"
#include <chrono>
#include <stdexcept>
#include <unordered_map>
#include <utility>
namespace coilgun::optimization {
namespace {
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
}

SerialBatchEvaluator::SerialBatchEvaluator(std::shared_ptr<const Evaluator> e) : evaluator_(std::move(e)) { if (!evaluator_) throw std::invalid_argument("evaluator must not be null"); }
std::vector<EvaluationResult> SerialBatchEvaluator::evaluate_batch(const std::vector<CandidateVariables>& c, const EvaluationContext& x) { std::vector<EvaluationResult> r; r.reserve(c.size()); for (const auto& v:c) r.push_back(one(*evaluator_,v,x)); return r; }
StatisticsBatchEvaluator::StatisticsBatchEvaluator(std::shared_ptr<BatchEvaluator> e) : evaluator_(std::move(e)) { if (!evaluator_) throw std::invalid_argument("evaluator must not be null"); }
std::vector<EvaluationResult> StatisticsBatchEvaluator::evaluate_batch(const std::vector<CandidateVariables>& c, const EvaluationContext& x) {
    const auto start = std::chrono::steady_clock::now();
    statistics_.seed = x.seed;
    if (x.fallback) ++statistics_.fallbacks;
    auto results = safe_batch(*evaluator_, c, x);
    statistics_.evaluations += c.size();
    for (const auto& result : results) record_status(statistics_, result);
    statistics_.elapsed_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return results;
}
CachedBatchEvaluator::CachedBatchEvaluator(std::shared_ptr<BatchEvaluator> e,std::shared_ptr<EvaluationCache> c):evaluator_(std::move(e)),cache_(std::move(c)){if(!evaluator_||!cache_)throw std::invalid_argument("evaluator and cache must not be null");}
std::vector<EvaluationResult> CachedBatchEvaluator::evaluate_batch(const std::vector<CandidateVariables>& c,const EvaluationContext& x){
    const auto start = std::chrono::steady_clock::now();
    std::vector<EvaluationResult> results;
    results.reserve(c.size());
    statistics_.seed = x.seed;
    if (x.fallback) ++statistics_.fallbacks;
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
    for (std::size_t i = 0; i < c.size(); ++i) {
        const auto& variables = c[i];
        const auto key = make_cache_key(variables, x);
        if (auto hit = cache_->get(key)) {
            results[i] = normalize_result(*hit);
            ++statistics_.cache_hits;
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
    const auto fresh = safe_batch_with_isolation(*evaluator_, misses, x);
    for (std::size_t i = 0; i < fresh.size(); ++i) {
        const auto result = normalize_result(fresh[i]);
        cache_->put(miss_keys[i], result);
        record_status(statistics_, result);
    }
    for (std::size_t i = 0; i < miss_candidate_indices.size(); ++i) {
        results[miss_candidate_indices[i]] = fresh[miss_result_indices[i]];
    }
    statistics_.evaluations += misses.size();
    statistics_.elapsed_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return results;
}
}
