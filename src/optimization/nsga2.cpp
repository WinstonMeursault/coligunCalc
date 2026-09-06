#include "coilgun/optimization/nsga2.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace coilgun::optimization {
namespace {

constexpr double inf = std::numeric_limits<double>::infinity();

struct ObjectiveView {
    std::size_t count = 0;
    const std::vector<ObjectiveDefinition>* definitions = nullptr;

    double oriented(const Candidate& candidate, std::size_t objective) const {
        if (candidate.evaluation_status != EvaluationStatus::Success ||
            objective >= candidate.objectives.size())
            return inf;
        const auto& value = candidate.objectives[objective];
        if (definitions != nullptr && !definitions->empty()) {
            try {
                return (*definitions)[objective].oriented(value.value);
            } catch (const std::invalid_argument&) {
                return inf;
            }
        }
        if (!std::isfinite(value.value)) return inf;
        return value.maximize ? -value.value : value.value;
    }
};

void validate_input(const std::vector<Candidate>& candidates,
                    const std::vector<ObjectiveDefinition>& definitions) {
    std::size_t count = definitions.size();
    bool has_success = false;
    if (!definitions.empty()) {
        for (const auto& definition : definitions) definition.validate();
    }

    // Failed and invalid evaluations are represented by empty objective lists.
    // Infer the schema only from successful candidates and leave an entirely
    // unsuccessful population objective-free (its candidates tie as infeasible).
    if (count == 0) {
        for (const auto& candidate : candidates) {
            if (candidate.evaluation_status == EvaluationStatus::Success) {
                has_success = true;
                count = candidate.objectives.size();
                break;
            }
        }
    }
    if (count != 0 && count < 2) throw std::invalid_argument("NSGA-II requires at least two objectives");
    if (has_success && count == 0) throw std::invalid_argument("NSGA-II requires at least two objectives");
    for (const auto& candidate : candidates) {
        if (candidate.evaluation_status != EvaluationStatus::Success) continue;
        if (candidate.objectives.size() != count)
            throw std::invalid_argument("all candidates must have the fixed objective count");
    }
}

std::size_t objective_count(const std::vector<Candidate>& candidates,
                            const std::vector<ObjectiveDefinition>& definitions) {
    if (!definitions.empty()) return definitions.size();
    for (const auto& candidate : candidates) {
        if (candidate.evaluation_status == EvaluationStatus::Success)
            return candidate.objectives.size();
    }
    return 0;
}

double hard_violation(const Candidate& candidate) {
    if (candidate.evaluation_status != EvaluationStatus::Success) return inf;
    const double violation = aggregate_normalized_violation(candidate.constraints, ConstraintKind::Hard);
    return std::isfinite(violation) && violation >= 0.0 ? violation : inf;
}

bool dominates(const Candidate& lhs, const Candidate& rhs, const ObjectiveView& view) {
    const double lhs_violation = hard_violation(lhs);
    const double rhs_violation = hard_violation(rhs);
    const bool lhs_feasible = lhs_violation == 0.0;
    const bool rhs_feasible = rhs_violation == 0.0;
    if (lhs_feasible != rhs_feasible) return lhs_feasible;
    if (!lhs_feasible) return lhs_violation < rhs_violation;

    bool strictly_better = false;
    for (std::size_t i = 0; i < view.count; ++i) {
        const double left = view.oriented(lhs, i);
        const double right = view.oriented(rhs, i);
        if (left > right) return false;
        if (left < right) strictly_better = true;
    }
    return strictly_better;
}

std::vector<double> crowding_for_front(const std::vector<Candidate>& candidates,
                                       const std::vector<std::size_t>& front,
                                       const ObjectiveView& view) {
    std::vector<double> distances(candidates.size(), 0.0);
    if (front.empty()) return distances;
    if (front.size() <= 2) {
        for (const auto index : front) distances[index] = inf;
        return distances;
    }

    for (std::size_t objective = 0; objective < view.count; ++objective) {
        std::vector<std::size_t> order = front;
        std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
            return view.oriented(candidates[lhs], objective) < view.oriented(candidates[rhs], objective);
        });
        const double minimum = view.oriented(candidates[order.front()], objective);
        const double maximum = view.oriented(candidates[order.back()], objective);
        if (!std::isfinite(minimum) || !std::isfinite(maximum) || maximum == minimum) {
            for (const auto index : order) distances[index] = inf;
            continue;
        }

        // Mark every tied boundary point. This avoids making duplicate points
        // dependent on an arbitrary sort position and keeps ties stable.
        for (const auto index : order) {
            const double value = view.oriented(candidates[index], objective);
            if (value == minimum || value == maximum) distances[index] = inf;
        }
        for (std::size_t position = 1; position + 1 < order.size(); ++position) {
            const auto index = order[position];
            if (std::isinf(distances[index])) continue;
            const double previous = view.oriented(candidates[order[position - 1]], objective);
            const double next = view.oriented(candidates[order[position + 1]], objective);
            distances[index] += (next - previous) / (maximum - minimum);
        }
    }
    return distances;
}

} // namespace

Nsga2Ranking nsga2_rank(const std::vector<Candidate>& candidates,
                        const std::vector<ObjectiveDefinition>& definitions) {
    validate_input(candidates, definitions);
    Nsga2Ranking result;
    result.ranks.assign(candidates.size(), 0);
    result.crowding_distances.assign(candidates.size(), 0.0);
    if (candidates.empty()) return result;

    const ObjectiveView view{objective_count(candidates, definitions),
                             definitions.empty() ? nullptr : &definitions};
    std::vector<std::vector<std::size_t>> dominated(candidates.size());
    std::vector<std::size_t> domination_count(candidates.size(), 0);
    std::vector<std::size_t> first;
    for (std::size_t lhs = 0; lhs < candidates.size(); ++lhs) {
        for (std::size_t rhs = lhs + 1; rhs < candidates.size(); ++rhs) {
            if (dominates(candidates[lhs], candidates[rhs], view)) {
                dominated[lhs].push_back(rhs);
                ++domination_count[rhs];
            } else if (dominates(candidates[rhs], candidates[lhs], view)) {
                dominated[rhs].push_back(lhs);
                ++domination_count[lhs];
            }
        }
        if (domination_count[lhs] == 0) first.push_back(lhs);
    }

    result.fronts.push_back(std::move(first));
    for (std::size_t rank = 0; rank < result.fronts.size(); ++rank) {
        for (const auto index : result.fronts[rank]) result.ranks[index] = rank;
        std::vector<std::size_t> next;
        for (const auto index : result.fronts[rank]) {
            for (const auto dominated_index : dominated[index]) {
                if (--domination_count[dominated_index] == 0) next.push_back(dominated_index);
            }
        }
        if (!next.empty()) result.fronts.push_back(std::move(next));
    }

    for (const auto& front : result.fronts) {
        const auto distances = crowding_for_front(candidates, front, view);
        for (const auto index : front) result.crowding_distances[index] = distances[index];
    }
    return result;
}

std::vector<std::vector<std::size_t>> non_dominated_sort(
    const std::vector<Candidate>& candidates,
    const std::vector<ObjectiveDefinition>& definitions) {
    return nsga2_rank(candidates, definitions).fronts;
}

std::vector<double> crowding_distances(
    const std::vector<Candidate>& candidates,
    const std::vector<std::size_t>& front,
    const std::vector<ObjectiveDefinition>& definitions) {
    validate_input(candidates, definitions);
    if (candidates.empty()) return {};
    const ObjectiveView view{objective_count(candidates, definitions),
                             definitions.empty() ? nullptr : &definitions};
    for (const auto index : front) {
        if (index >= candidates.size()) throw std::out_of_range("front index out of range");
    }
    return crowding_for_front(candidates, front, view);
}

std::vector<Candidate> nsga2_select(
    const std::vector<Candidate>& parents,
    const std::vector<Candidate>& offspring,
    std::size_t target_size,
    const std::vector<ObjectiveDefinition>& definitions) {
    std::vector<Candidate> merged;
    merged.reserve(parents.size() + offspring.size());
    merged.insert(merged.end(), parents.begin(), parents.end());
    merged.insert(merged.end(), offspring.begin(), offspring.end());
    if (merged.empty() || target_size == 0) return {};
    const auto ranking = nsga2_rank(merged, definitions);
    const std::size_t count = std::min(target_size, merged.size());
    std::vector<std::size_t> selected;
    selected.reserve(count);
    for (const auto& front : ranking.fronts) {
        if (selected.size() + front.size() <= count) {
            selected.insert(selected.end(), front.begin(), front.end());
            continue;
        }
        std::vector<std::size_t> remainder = front;
        std::stable_sort(remainder.begin(), remainder.end(), [&](std::size_t lhs, std::size_t rhs) {
            return ranking.crowding_distances[lhs] > ranking.crowding_distances[rhs];
        });
        remainder.resize(count - selected.size());
        selected.insert(selected.end(), remainder.begin(), remainder.end());
        break;
    }
    std::vector<Candidate> result;
    result.reserve(selected.size());
    for (const auto index : selected) result.push_back(merged[index]);
    return result;
}

Population nsga2_select(
    const Population& parents,
    const Population& offspring,
    std::size_t target_size,
    const std::vector<ObjectiveDefinition>& definitions) {
    std::vector<Candidate> parent_values(parents.begin(), parents.end());
    std::vector<Candidate> offspring_values(offspring.begin(), offspring.end());
    const auto selected = nsga2_select(parent_values, offspring_values, target_size, definitions);
    Population result;
    for (auto candidate : selected) result.push_back(std::move(candidate));
    return result;
}

} // namespace coilgun::optimization
