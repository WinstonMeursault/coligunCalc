#include "coilgun/optimization/constraint.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace coilgun::optimization {
void ConstraintDefinition::validate() const {
    if (id.empty()) throw std::invalid_argument("constraint id must not be empty");
    if (!std::isfinite(scale) || scale <= 0.0) throw std::invalid_argument("constraint scale must be finite and positive");
    if (!std::isfinite(lower_bound) || !std::isfinite(upper_bound)) throw std::invalid_argument("constraint bounds must be finite");
    if (relation == ConstraintRelation::InRange && lower_bound > upper_bound)
        throw std::invalid_argument("range lower bound cannot exceed upper bound");
}
ConstraintReport ConstraintDefinition::evaluate(double value) const {
    validate();
    if (!std::isfinite(value)) throw std::invalid_argument("constraint value must be finite");
    double violation = 0.0;
    switch (relation) {
    case ConstraintRelation::Equal: violation = std::abs(value - lower_bound); break;
    case ConstraintRelation::LessEqual: violation = std::max(0.0, value - upper_bound); break;
    case ConstraintRelation::GreaterEqual: violation = std::max(0.0, lower_bound - value); break;
    case ConstraintRelation::InRange:
        violation = value < lower_bound ? lower_bound - value : (value > upper_bound ? value - upper_bound : 0.0); break;
    }
    ConstraintReport report{id, kind, relation, value, lower_bound, upper_bound, violation,
                           violation / scale, violation == 0.0, priority};
    return report;
}
double aggregate_normalized_violation(const std::vector<ConstraintReport>& constraints, ConstraintKind kind) {
    double total = 0.0;
    for (const auto& c : constraints) if (c.kind == kind) total += c.normalized_violation;
    return total;
}
bool is_feasible(const std::vector<ConstraintReport>& constraints) {
    return aggregate_normalized_violation(constraints, ConstraintKind::Hard) == 0.0;
}
}
