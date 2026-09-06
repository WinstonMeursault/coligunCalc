#include "coilgun/optimization/objective.hpp"
#include <cmath>
#include <stdexcept>

namespace coilgun::optimization {
void ObjectiveDefinition::validate() const {
    if (id.empty()) throw std::invalid_argument("objective id must not be empty");
    if (!std::isfinite(scale) || scale <= 0.0) throw std::invalid_argument("objective scale must be finite and positive");
}
double ObjectiveDefinition::normalize(double value) const {
    if (!std::isfinite(value)) throw std::invalid_argument("objective value must be finite");
    validate();
    return value / scale;
}
double ObjectiveDefinition::oriented(double value) const { return maximize ? -normalize(value) : normalize(value); }
}
