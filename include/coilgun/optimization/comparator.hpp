#pragma once

#include "coilgun/optimization/objective.hpp"
#include "coilgun/optimization/constraint.hpp"

namespace coilgun::optimization {

enum class FeasibilityStrategy { FeasibilityFirst, Penalty, Lexicographic };

class FeasibilityComparator {
public:
    explicit FeasibilityComparator(FeasibilityStrategy strategy = FeasibilityStrategy::FeasibilityFirst,
                                   double penalty_weight = 1.0);
    int compare(const Candidate& lhs, const Candidate& rhs) const;
    bool better(const Candidate& lhs, const Candidate& rhs) const { return compare(lhs, rhs) < 0; }

private:
    FeasibilityStrategy strategy_;
    double penalty_weight_;
};

}
