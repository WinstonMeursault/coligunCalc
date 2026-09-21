#pragma once
#include "coilgun/optimization/statistics.hpp"
#include "coilgun/optimization/types.hpp"
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
namespace coilgun::optimization {
struct EvaluationContext {
    std::uint64_t seed = 0;
    bool fallback = false;
    // Optional run-local sink. Existing {seed, fallback} construction remains valid.
    std::shared_ptr<EvaluationStatisticsCollector> statistics;

    [[nodiscard]] const std::shared_ptr<EvaluationStatisticsCollector>& statistics_collector() const noexcept {
        return statistics;
    }
};

/// Stable, caller-supplied identity for an evaluator's cached result schema.
/// Namespace and version must be non-empty; increment the version when result semantics change.
class EvaluationCacheIdentity {
public:
    EvaluationCacheIdentity(std::string namespace_id, std::string version);

    // Public const fields preserve the original read-only access shape while
    // preventing mutation after construction.
    const std::string namespace_id;
    const std::string version;
};

class EvaluationCache {
public: virtual ~EvaluationCache() = default;
    virtual std::optional<EvaluationResult> get(const std::string&) const = 0;
    virtual void put(std::string, EvaluationResult) = 0;
};
class InMemoryEvaluationCache final : public EvaluationCache {
public: std::optional<EvaluationResult> get(const std::string&) const override;
    void put(std::string, EvaluationResult) override;
private:
    mutable std::mutex mutex_;
    std::unordered_map<std::string, EvaluationResult> values_;
};
std::string make_cache_key(const EvaluationCacheIdentity&, const CandidateVariables&,
                           const EvaluationContext&);
/// Builds a key using the legacy `coilgun.optimization.legacy` namespace at version `1`.
std::string make_cache_key(const CandidateVariables&, const EvaluationContext&);
}
