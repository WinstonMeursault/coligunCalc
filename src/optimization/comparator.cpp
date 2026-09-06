#include "coilgun/optimization/comparator.hpp"
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace coilgun::optimization {
namespace {
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
int FeasibilityComparator::compare(const Candidate& lhs, const Candidate& rhs) const {
    const bool lhs_success = lhs.evaluation_status == EvaluationStatus::Success;
    const bool rhs_success = rhs.evaluation_status == EvaluationStatus::Success;
    const double lh = aggregate_normalized_violation(lhs.constraints, ConstraintKind::Hard);
    const double rh = aggregate_normalized_violation(rhs.constraints, ConstraintKind::Hard);
    const bool lf = lhs_success && lh == 0.0, rf = rhs_success && rh == 0.0;
    if (strategy_ == FeasibilityStrategy::FeasibilityFirst || strategy_ == FeasibilityStrategy::Lexicographic) {
        if (lf != rf) return lf ? -1 : 1;
        if (!lf) {
            if (strategy_ == FeasibilityStrategy::Lexicographic) {
                if (lhs_success != rhs_success) return lhs_success ? -1 : 1;
                const int hard_comparison = compare_constraint_priorities(lhs.constraints, rhs.constraints,
                                                                           ConstraintKind::Hard);
                if (hard_comparison != 0) return hard_comparison;
            } else if (lh != rh) {
                const double lhs_violation = lhs_success ? lh : std::numeric_limits<double>::infinity();
                const double rhs_violation = rhs_success ? rh : std::numeric_limits<double>::infinity();
                if (lhs_violation != rhs_violation) return lhs_violation < rhs_violation ? -1 : 1;
            }
        }
        if (strategy_ == FeasibilityStrategy::Lexicographic) {
            const int soft_comparison = compare_constraint_priorities(lhs.constraints, rhs.constraints,
                                                                        ConstraintKind::Soft);
            if (soft_comparison != 0) return soft_comparison;
        }
    } else {
        const double lp = lhs_success
            ? lh + penalty_weight_ * aggregate_normalized_violation(lhs.constraints, ConstraintKind::Soft)
            : std::numeric_limits<double>::infinity();
        const double rp = rhs_success
            ? rh + penalty_weight_ * aggregate_normalized_violation(rhs.constraints, ConstraintKind::Soft)
            : std::numeric_limits<double>::infinity();
        if (lp != rp) return lp < rp ? -1 : 1;
    }
    if (lhs.evaluation_status != rhs.evaluation_status) {
        const int lhs_rank = status_rank(lhs.evaluation_status);
        const int rhs_rank = status_rank(rhs.evaluation_status);
        if (lhs_rank != rhs_rank) return lhs_rank < rhs_rank ? -1 : 1;
    }
    if (!lhs.objectives.empty() && !rhs.objectives.empty()) {
        ObjectiveDefinition d{lhs.objectives.front().id, lhs.objectives.front().maximize, 1.0};
        const double lo = d.oriented(lhs.objectives.front().value), ro = d.oriented(rhs.objectives.front().value);
        if (lo != ro) return lo < ro ? -1 : 1;
    }
    return 0;
}
}
