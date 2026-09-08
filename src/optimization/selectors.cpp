#include "coilgun/optimization/selectors.hpp"

#include "coilgun/optimization/constraint.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace coilgun::optimization {
namespace {

const ObjectiveValue& objective(const Candidate& candidate, const std::string& id) {
    const auto it = std::find_if(candidate.objectives.begin(), candidate.objectives.end(),
                                 [&](const ObjectiveValue& value) { return value.id == id; });
    if (it == candidate.objectives.end()) throw std::invalid_argument("objective not found: " + id);
    if (!std::isfinite(it->value)) throw std::invalid_argument("objective value is not finite: " + id);
    return *it;
}

double oriented(const ObjectiveValue& value) { return value.maximize ? value.value : -value.value; }

double normalized_violation(const Candidate& candidate) {
    if (candidate.evaluation_status != EvaluationStatus::Success)
        return std::numeric_limits<double>::infinity();
    const double violation =
        aggregate_normalized_violation(candidate.constraints, ConstraintKind::Hard) +
        aggregate_normalized_violation(candidate.constraints, ConstraintKind::Soft);
    return std::isfinite(violation) && violation >= 0.0
        ? violation
        : std::numeric_limits<double>::infinity();
}

template <typename Better>
std::optional<Candidate> choose(const OptimizationResult& result, Better better) {
    if (result.pareto_front.empty()) return std::nullopt;
    std::size_t best = 0;
    for (std::size_t i = 1; i < result.pareto_front.size(); ++i)
        if (better(result.pareto_front[i], result.pareto_front[best])) best = i;
    return result.pareto_front[best];
}

std::vector<std::vector<double>> normalized_objectives(const OptimizationResult& result) {
    if (result.pareto_front.empty()) return {};
    const auto count = result.pareto_front.front().objectives.size();
    if (count == 0) throw std::invalid_argument("Pareto candidates have no objectives");
    std::vector<double> minimum(count, std::numeric_limits<double>::infinity());
    std::vector<double> maximum(count, -std::numeric_limits<double>::infinity());
    for (const auto& candidate : result.pareto_front) {
        if (candidate.objectives.size() != count)
            throw std::invalid_argument("Pareto candidates have inconsistent objective counts");
        for (std::size_t i = 0; i < count; ++i) {
            const double value = oriented(candidate.objectives[i]);
            if (!std::isfinite(value)) throw std::invalid_argument("objective value is not finite");
            minimum[i] = std::min(minimum[i], value);
            maximum[i] = std::max(maximum[i], value);
        }
    }
    std::vector<std::vector<double>> normalized(result.pareto_front.size(), std::vector<double>(count));
    for (std::size_t row = 0; row < result.pareto_front.size(); ++row) {
        for (std::size_t i = 0; i < count; ++i) {
            const double range = maximum[i] - minimum[i];
            normalized[row][i] = range > 0.0
                ? (oriented(result.pareto_front[row].objectives[i]) - minimum[i]) / range
                : 0.5;
        }
    }
    return normalized;
}

}

std::optional<Candidate> MaxObjective::select(const OptimizationResult& result) const {
    return choose(result, [&](const Candidate& lhs, const Candidate& rhs) {
        return oriented(objective(lhs, objective_id_)) > oriented(objective(rhs, objective_id_));
    });
}

std::optional<Candidate> MinConstraintViolationMargin::select(const OptimizationResult& result) const {
    return choose(result, [](const Candidate& lhs, const Candidate& rhs) {
        return normalized_violation(lhs) < normalized_violation(rhs);
    });
}

std::optional<Candidate> IdealPointDistance::select(const OptimizationResult& result) const {
    if (result.pareto_front.empty()) return std::nullopt;
    if (!std::isfinite(distance_power_) || distance_power_ <= 0.0)
        throw std::invalid_argument("distance power must be finite and positive");
    const auto values = normalized_objectives(result);
    return choose(result, [&](const Candidate& lhs, const Candidate& rhs) {
        const auto index = static_cast<std::size_t>(&lhs - result.pareto_front.data());
        const auto other = static_cast<std::size_t>(&rhs - result.pareto_front.data());
        double left = 0.0, right = 0.0;
        for (const double value : values[index]) left += std::pow(1.0 - value, distance_power_);
        for (const double value : values[other]) right += std::pow(1.0 - value, distance_power_);
        return left < right;
    });
}

std::optional<Candidate> WeightedScore::select(const OptimizationResult& result) const {
    if (result.pareto_front.empty()) return std::nullopt;
    if (weights_.empty() || std::any_of(weights_.begin(), weights_.end(),
                                        [](double weight) { return !std::isfinite(weight) || weight < 0.0; }))
        throw std::invalid_argument("weights must be finite, non-negative, and non-empty");
    const auto values = normalized_objectives(result);
    if (!values.empty() && weights_.size() != values.front().size())
        throw std::invalid_argument("weight count must match objective count");
    return choose(result, [&](const Candidate& lhs, const Candidate& rhs) {
        const auto left_index = static_cast<std::size_t>(&lhs - result.pareto_front.data());
        const auto right_index = static_cast<std::size_t>(&rhs - result.pareto_front.data());
        double left = 0.0, right = 0.0;
        for (std::size_t i = 0; i < weights_.size(); ++i) {
            left += weights_[i] * values[left_index][i];
            right += weights_[i] * values[right_index][i];
        }
        return left > right;
    });
}

std::optional<Candidate> LexicographicObjectives::select(const OptimizationResult& result) const {
    if (result.pareto_front.empty()) return std::nullopt;
    if (objective_ids_.empty()) throw std::invalid_argument("lexicographic objective list cannot be empty");
    return choose(result, [&](const Candidate& lhs, const Candidate& rhs) {
        for (const auto& id : objective_ids_) {
            const double left = oriented(objective(lhs, id));
            const double right = oriented(objective(rhs, id));
            if (left != right) return left > right;
        }
        return false;
    });
}

std::optional<Candidate> CallbackSelector::select(const OptimizationResult& result) const {
    if (result.pareto_front.empty()) return std::nullopt;
    if (!callback_) throw std::invalid_argument("custom selector callback is empty");
    return callback_(result);
}

}
