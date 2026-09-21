#include "coilgun/optimization/cache.hpp"
#include <cstring>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>
namespace coilgun::optimization {
EvaluationCacheIdentity::EvaluationCacheIdentity(std::string namespace_value,
                                                 std::string version_value)
    : namespace_id(std::move(namespace_value)), version(std::move(version_value)) {
    if (namespace_id.empty()) throw std::invalid_argument("cache identity namespace must not be empty");
    if (version.empty()) throw std::invalid_argument("cache identity version must not be empty");
}

std::optional<EvaluationResult> InMemoryEvaluationCache::get(const std::string& key) const {
    std::lock_guard lock(mutex_);
    auto it = values_.find(key); return it == values_.end() ? std::nullopt : std::optional<EvaluationResult>(it->second);
}
void InMemoryEvaluationCache::put(std::string key, EvaluationResult result) {
    std::lock_guard lock(mutex_);
    values_[std::move(key)] = std::move(result);
}
std::string make_cache_key(const EvaluationCacheIdentity& identity,
                           const CandidateVariables& variables,
                           const EvaluationContext& context) {
    std::ostringstream out;
    out << identity.namespace_id.size() << ':' << identity.namespace_id << ':'
        << identity.version.size() << ':' << identity.version << ':' << context.seed << ':'
        << context.fallback << ':';
    for (double value : variables.values) { std::uint64_t bits; std::memcpy(&bits, &value, sizeof bits); out << std::hex << bits << ','; }
    return out.str();
}

std::string make_cache_key(const CandidateVariables& variables, const EvaluationContext& context) {
    static const EvaluationCacheIdentity legacy_identity{"coilgun.optimization.legacy", "1"};
    return make_cache_key(legacy_identity, variables, context);
}
}
