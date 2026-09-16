#include "coilgun/optimization/cache.hpp"
#include <cstring>
#include <iomanip>
#include <sstream>
namespace coilgun::optimization {
std::optional<EvaluationResult> InMemoryEvaluationCache::get(const std::string& key) const {
    auto it = values_.find(key); return it == values_.end() ? std::nullopt : std::optional<EvaluationResult>(it->second);
}
void InMemoryEvaluationCache::put(std::string key, EvaluationResult result) { values_[std::move(key)] = std::move(result); }
std::string make_cache_key(const CandidateVariables& variables, const EvaluationContext& context) {
    std::ostringstream out; out << context.seed << ':' << context.fallback << ':';
    for (double value : variables.values) { std::uint64_t bits; std::memcpy(&bits, &value, sizeof bits); out << std::hex << bits << ','; }
    return out.str();
}
}
