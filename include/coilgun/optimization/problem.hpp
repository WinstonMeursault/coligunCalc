#pragma once

#include "coilgun/optimization/types.hpp"

namespace coilgun::optimization {

class OptimizationProblem {
public:
    virtual ~OptimizationProblem() = default;
    virtual EvaluationResult evaluate(const CandidateVariables& variables) const = 0;
};

} // namespace coilgun::optimization
