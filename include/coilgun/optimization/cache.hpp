#pragma once
#include "coilgun/optimization/statistics.hpp"
#include "coilgun/optimization/types.hpp"
#include <optional>
#include <string>
#include <unordered_map>
namespace coilgun::optimization {
struct EvaluationContext { std::uint64_t seed = 0; bool fallback = false; };
class EvaluationCache {
public: virtual ~EvaluationCache() = default;
    virtual std::optional<EvaluationResult> get(const std::string&) const = 0;
    virtual void put(std::string, EvaluationResult) = 0;
};
class InMemoryEvaluationCache final : public EvaluationCache {
public: std::optional<EvaluationResult> get(const std::string&) const override;
    void put(std::string, EvaluationResult) override;
private: std::unordered_map<std::string, EvaluationResult> values_;
};
std::string make_cache_key(const CandidateVariables&, const EvaluationContext&);
}
