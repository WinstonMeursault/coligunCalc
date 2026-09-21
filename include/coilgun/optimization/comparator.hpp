#pragma once

#include "coilgun/optimization/objective.hpp"
#include "coilgun/optimization/constraint.hpp"

#include <optional>

namespace coilgun::optimization {

enum class FeasibilityStrategy { FeasibilityFirst, Penalty, Lexicographic };

class FeasibilityComparator {
public:
    explicit FeasibilityComparator(FeasibilityStrategy strategy = FeasibilityStrategy::FeasibilityFirst,
                                   double penalty_weight = 1.0);
    FeasibilityComparator(FeasibilityStrategy strategy, double penalty_weight,
                          ObjectiveDefinition objective_definition);
    int compare(const Candidate& lhs, const Candidate& rhs) const;
    int compare(const Candidate& lhs, const Candidate& rhs,
                const ObjectiveDefinition& objective_definition) const;
    bool better(const Candidate& lhs, const Candidate& rhs) const { return compare(lhs, rhs) < 0; }
    bool better(const Candidate& lhs, const Candidate& rhs,
                const ObjectiveDefinition& objective_definition) const {
        return compare(lhs, rhs, objective_definition) < 0;
    }

    FeasibilityStrategy strategy() const { return strategy_; }
    double penalty_weight() const { return penalty_weight_; }

private:
    int compare_impl(const Candidate& lhs, const Candidate& rhs,
                     const ObjectiveDefinition* objective_definition) const;

    FeasibilityStrategy strategy_;
    double penalty_weight_;
    std::optional<ObjectiveDefinition> objective_definition_;
};

}
