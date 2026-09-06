#include "coilgun/optimization/evaluator.hpp"
#include <chrono>
#include <stdexcept>
namespace coilgun::optimization {
namespace {
EvaluationResult exception_result(const std::exception& error) {
    return EvaluationResult::failed("evaluation_exception", error.what());
}

EvaluationResult unknown_exception_result() {
    return EvaluationResult::failed("evaluation_exception", "unknown exception");
}

std::vector<EvaluationResult> malformed_batch(std::size_t count, const char* message) {
    std::vector<EvaluationResult> results;
    results.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        results.push_back(EvaluationResult::failed("evaluation_batch_output", message));
    }
    return results;
}

std::vector<EvaluationResult> safe_batch(BatchEvaluator& evaluator,
                                         const std::vector<CandidateVariables>& candidates,
                                         const EvaluationContext& context) {
    if (candidates.empty()) return {};
    try {
        auto results = evaluator.evaluate_batch(candidates, context);
        if (results.size() != candidates.size()) {
            return malformed_batch(candidates.size(), "batch evaluator returned an unexpected number of results");
        }
        return results;
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

EvaluationResult one(const Evaluator& evaluator, const CandidateVariables& variables,
                     const EvaluationContext& context) {
    try {
        return evaluator.evaluate(variables, context);
    } catch (const std::exception& error) {
        return exception_result(error);
    } catch (...) {
        return unknown_exception_result();
    }
}

void record_status(EvaluationStatistics& statistics, const EvaluationResult& result) {
    if (result.status == EvaluationStatus::Success) {
        ++statistics.successful_evaluations;
    } else if (result.status == EvaluationStatus::Failed || result.status == EvaluationStatus::Invalid) {
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
    std::uint64_t misses = 0;
    statistics_.seed = x.seed;
    if (x.fallback) ++statistics_.fallbacks;
    for (const auto& variables : c) {
        const auto key = make_cache_key(variables, x);
        if (auto hit = cache_->get(key)) {
            results.push_back(*hit);
            ++statistics_.cache_hits;
            continue;
        }
        const auto fresh = safe_batch(*evaluator_, {variables}, x);
        const auto& result = fresh.front();
        results.push_back(result);
        cache_->put(key, result);
        ++misses;
        record_status(statistics_, result);
    }
    statistics_.evaluations += misses;
    statistics_.elapsed_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return results;
}
}
