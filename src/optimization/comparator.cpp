#include "coilgun/optimization/comparator.hpp"
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace coilgun::optimization {
namespace {
double checked_aggregate_normalized_violation(const Candidate& candidate, ConstraintKind kind) {
    for (const auto& constraint : candidate.constraints) {
        if (!std::isfinite(constraint.normalized_violation) || constraint.normalized_violation < 0.0)
            throw std::invalid_argument(
                "normalized constraint violation must be finite and non-negative");
    }
    return aggregate_normalized_violation(candidate.constraints, kind);
}

int compare_constraint_priorities(const std::vector<ConstraintReport>& lhs,
                                  const std::vector<ConstraintReport>& rhs,
                                  ConstraintKind kind) {
    std::map<int, std::pair<double, double>> levels;
    for (const auto& report : lhs) {
        if (report.kind == kind) levels[report.priority].first += report.normalized_violation;
    }
    for (const auto& report : rhs) {
        if (report.kind == kind) levels[report.priority].second += report.normalized_violation;
    }
    for (const auto& [priority, violations] : levels) {
        (void)priority;
        if (violations.first != violations.second)
            return violations.first < violations.second ? -1 : 1;
    }
    return 0;
}

int status_rank(EvaluationStatus status) {
    switch (status) {
    case EvaluationStatus::Success: return 0;
    case EvaluationStatus::Invalid: return 1;
    case EvaluationStatus::Failed: return 2;
    case EvaluationStatus::Unevaluated: return 3;
    }
    return 3;
}
}

FeasibilityComparator::FeasibilityComparator(FeasibilityStrategy strategy, double penalty_weight)
    : strategy_(strategy), penalty_weight_(penalty_weight) {
    if (!std::isfinite(penalty_weight) || penalty_weight < 0.0) throw std::invalid_argument("penalty weight must be finite and non-negative");
}
FeasibilityComparator::FeasibilityComparator(FeasibilityStrategy strategy, double penalty_weight,
                                               ObjectiveDefinition objective_definition)
    : FeasibilityComparator(strategy, penalty_weight) {
    objective_definition.validate();
    objective_definition_ = std::move(objective_definition);
}

int FeasibilityComparator::compare(const Candidate& lhs, const Candidate& rhs) const {
    return compare_impl(lhs, rhs, objective_definition_ ? &*objective_definition_ : nullptr);
}

int FeasibilityComparator::compare(const Candidate& lhs, const Candidate& rhs,
                                   const ObjectiveDefinition& objective_definition) const {
    objective_definition.validate();
    return compare_impl(lhs, rhs, &objective_definition);
}

int FeasibilityComparator::compare_impl(const Candidate& lhs, const Candidate& rhs,
                                        const ObjectiveDefinition* objective_definition) const {
    const double lh = checked_aggregate_normalized_violation(lhs, ConstraintKind::Hard);
    const double rh = checked_aggregate_normalized_violation(rhs, ConstraintKind::Hard);
    const double ls = checked_aggregate_normalized_violation(lhs, ConstraintKind::Soft);
    const double rs = checked_aggregate_normalized_violation(rhs, ConstraintKind::Soft);

    const int lhs_rank = status_rank(lhs.evaluation_status);
    const int rhs_rank = status_rank(rhs.evaluation_status);
    if (lhs_rank != rhs_rank) return lhs_rank < rhs_rank ? -1 : 1;
    const bool successful = lhs.evaluation_status == EvaluationStatus::Success;
    if (!successful) {
        if (strategy_ == FeasibilityStrategy::Lexicographic) {
            const int hard_comparison = compare_constraint_priorities(lhs.constraints, rhs.constraints,
                                                                       ConstraintKind::Hard);
            if (hard_comparison != 0) return hard_comparison;
            const int soft_comparison = compare_constraint_priorities(lhs.constraints, rhs.constraints,
                                                                       ConstraintKind::Soft);
            if (soft_comparison != 0) return soft_comparison;
        }
    } else {
        const bool lf = lh == 0.0;
        const bool rf = rh == 0.0;
        if (lf != rf) return lf ? -1 : 1;
        if (!lf && lh != rh) return lh < rh ? -1 : 1;

        if (strategy_ == FeasibilityStrategy::Lexicographic) {
            if (!lf) {
                const int hard_comparison = compare_constraint_priorities(lhs.constraints, rhs.constraints,
                                                                           ConstraintKind::Hard);
                if (hard_comparison != 0) return hard_comparison;
            }
            const int soft_comparison = compare_constraint_priorities(lhs.constraints, rhs.constraints,
                                                                        ConstraintKind::Soft);
            if (soft_comparison != 0) return soft_comparison;
        }
    }

    // A cleared objective array denotes a constraint-only comparison. In that
    // mode even Penalty must not turn soft constraints into a prefilter.
    if (!lhs.objectives.empty() && !rhs.objectives.empty()) {
        const auto oriented = [&](const Candidate& candidate) {
            const auto& objective = candidate.objectives.front();
            if (objective_definition != nullptr) return objective_definition->oriented(objective.value);
            ObjectiveDefinition implicit{objective.id, objective.maximize, 1.0};
            return implicit.oriented(objective.value);
        };
        const auto penalized = [&](const Candidate& candidate, double soft_violation) {
            const double base = oriented(candidate);
            return successful && strategy_ == FeasibilityStrategy::Penalty
                ? base + penalty_weight_ * soft_violation
                : base;
        };
        const double lhs_score = penalized(lhs, ls);
        const double rhs_score = penalized(rhs, rs);
        if (lhs_score != rhs_score) return lhs_score < rhs_score ? -1 : 1;
    }
    return 0;
}
}
