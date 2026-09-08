#pragma once
#include "coilgun/optimization/cache.hpp"
#include <memory>
#include <optional>
#include <vector>
namespace coilgun::optimization {
class Evaluator {
public: virtual ~Evaluator() = default;
    virtual EvaluationResult evaluate(const CandidateVariables&, const EvaluationContext&) const = 0;
};
class BatchEvaluator {
public: virtual ~BatchEvaluator() = default;
    virtual std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&, const EvaluationContext&) = 0;
    [[nodiscard]] virtual std::optional<EvaluationStatistics> statistics_snapshot() const { return std::nullopt; }
};
class SerialBatchEvaluator final : public BatchEvaluator {
public: explicit SerialBatchEvaluator(std::shared_ptr<const Evaluator>);
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&, const EvaluationContext&) override;
private: std::shared_ptr<const Evaluator> evaluator_;
};
class StatisticsBatchEvaluator final : public BatchEvaluator {
public: explicit StatisticsBatchEvaluator(std::shared_ptr<BatchEvaluator>);
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&, const EvaluationContext&) override;
    const EvaluationStatistics& statistics() const { return statistics_; }
    [[nodiscard]] std::optional<EvaluationStatistics> statistics_snapshot() const override;
private: std::shared_ptr<BatchEvaluator> evaluator_; EvaluationStatistics statistics_;
};
class CachedBatchEvaluator final : public BatchEvaluator {
public: CachedBatchEvaluator(std::shared_ptr<BatchEvaluator>, std::shared_ptr<EvaluationCache>);
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&, const EvaluationContext&) override;
    const EvaluationStatistics& statistics() const { return statistics_; }
    [[nodiscard]] std::optional<EvaluationStatistics> statistics_snapshot() const override;
private: std::shared_ptr<BatchEvaluator> evaluator_; std::shared_ptr<EvaluationCache> cache_; EvaluationStatistics statistics_;
};
}
