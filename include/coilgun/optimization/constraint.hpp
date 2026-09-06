#pragma once

#include "coilgun/optimization/types.hpp"
#include <string>

namespace coilgun::optimization {

struct ConstraintDefinition {
    std::string id;
    ConstraintKind kind = ConstraintKind::Hard;
    ConstraintRelation relation = ConstraintRelation::LessEqual;
    double lower_bound = 0.0;
    double upper_bound = 0.0;
    double scale = 1.0;
    // Lower values represent higher precedence in Lexicographic comparisons.
    int priority = 0;

    void validate() const;
    ConstraintReport evaluate(double value) const;
};

double aggregate_normalized_violation(const std::vector<ConstraintReport>& constraints,
                                      ConstraintKind kind = ConstraintKind::Hard);
bool is_feasible(const std::vector<ConstraintReport>& constraints);

}
