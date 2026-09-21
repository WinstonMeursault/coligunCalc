#pragma once
#include "coilgun/optimization/cache.hpp"
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <utility>
#include <vector>
namespace coilgun::optimization {
class Evaluator {
public: virtual ~Evaluator() = default;
    virtual EvaluationResult evaluate(const CandidateVariables&, const EvaluationContext&) const = 0;
};
struct BatchEvaluationSnapshot {
    EvaluationCacheIdentity identity;
    std::function<std::vector<EvaluationResult>(const std::vector<CandidateVariables>&,
                                                const EvaluationContext&)> evaluate;
};
class BatchEvaluator : public std::enable_shared_from_this<BatchEvaluator> {
public:
    BatchEvaluator();
    virtual ~BatchEvaluator() = default;
    virtual std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&, const EvaluationContext&) = 0;
    /** Optional safe batch hook for callers that only have a const evaluator. */
    virtual std::optional<std::vector<EvaluationResult>> evaluate_batch_const(
        const std::vector<CandidateVariables>&, const EvaluationContext&) const {
        return std::nullopt;
    }
    /// Returns the stable cache namespace and result-schema version owned by this evaluator.
    /// The legacy identity preserves source compatibility for existing evaluator subclasses.
    [[nodiscard]] virtual EvaluationCacheIdentity cache_identity() const {
        return {"coilgun.optimization.legacy", "1"};
    }
    /**
     * Final lifecycle boundary for an escaping identity/callable snapshot.
     * The evaluator must be managed by std::shared_ptr; the returned callable
     * retains that owner. Legacy subclasses retain the old API through the
     * default protected hook below; stateful evaluators should override the
     * hook to bind identity and state together.
     */
    virtual BatchEvaluationSnapshot evaluation_snapshot() final;
    [[nodiscard]] virtual std::optional<EvaluationStatistics> statistics_snapshot() const { return std::nullopt; }
protected:
    /** Builds one identity/callable pair for evaluation_snapshot(). */
    virtual BatchEvaluationSnapshot make_evaluation_snapshot();
private:
    BatchEvaluationSnapshot make_owned_snapshot(BatchEvaluationSnapshot,
                                                 std::shared_ptr<BatchEvaluator>);
};
class SerialBatchEvaluator final : public BatchEvaluator {
public: explicit SerialBatchEvaluator(std::shared_ptr<const Evaluator>);
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&, const EvaluationContext&) override;
private: std::shared_ptr<const Evaluator> evaluator_;
};
class StatisticsBatchEvaluator final : public BatchEvaluator {
public: explicit StatisticsBatchEvaluator(std::shared_ptr<BatchEvaluator>);
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&, const EvaluationContext&) override;
    [[nodiscard]] EvaluationCacheIdentity cache_identity() const override;
    [[nodiscard]] EvaluationStatistics statistics() const;
    [[nodiscard]] std::optional<EvaluationStatistics> statistics_snapshot() const override;
protected:
    BatchEvaluationSnapshot make_evaluation_snapshot() override;
private:
    std::shared_ptr<BatchEvaluator> evaluator_;
    mutable std::mutex statistics_mutex_;
    EvaluationStatistics statistics_;
};
class CachedBatchEvaluator final : public BatchEvaluator {
public: CachedBatchEvaluator(std::shared_ptr<BatchEvaluator>, std::shared_ptr<EvaluationCache>);
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&, const EvaluationContext&) override;
    [[nodiscard]] EvaluationCacheIdentity cache_identity() const override;
    [[nodiscard]] EvaluationStatistics statistics() const;
    [[nodiscard]] std::optional<EvaluationStatistics> statistics_snapshot() const override;
protected:
    BatchEvaluationSnapshot make_evaluation_snapshot() override;
private:
    std::vector<EvaluationResult> evaluate_batch_with_snapshot(
        const BatchEvaluationSnapshot&, const std::vector<CandidateVariables>&,
        const EvaluationContext&);
    std::shared_ptr<BatchEvaluator> evaluator_;
    std::shared_ptr<EvaluationCache> cache_;
    mutable std::mutex statistics_mutex_;
    EvaluationStatistics statistics_;
};
}
