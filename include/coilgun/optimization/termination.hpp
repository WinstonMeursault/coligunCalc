#pragma once

#include "coilgun/optimization/types.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>

namespace coilgun::optimization {

// Compatibility name retained for callers of the original optimizer API.
using GeneticTerminationReason = TerminationReason;

struct TerminationConfig {
    std::size_t max_generations = 0; // zero uses OptimizationConfig::max_generations
    std::size_t max_evaluations = 0; // zero means unlimited
    std::size_t max_no_improvement_generations = 0; // zero disables this check
    std::optional<double> target_value;
    double improvement_tolerance = 0.0;

    void validate() const {
        if (target_value && !std::isfinite(*target_value))
            throw std::invalid_argument("target_value must be finite");
        if (!std::isfinite(improvement_tolerance) || improvement_tolerance < 0.0)
            throw std::invalid_argument("improvement_tolerance must be finite and non-negative");
    }
};

using TerminationPolicy = TerminationConfig;

inline const char* to_string(TerminationReason reason) {
    switch (reason) {
    case TerminationReason::None: return "none";
    case TerminationReason::MaxGenerations: return "maximum generations";
    case TerminationReason::TargetReached: return "target reached";
    case TerminationReason::Converged: return "no improvement";
    case TerminationReason::Cancelled: return "cancelled";
    case TerminationReason::ConfigurationError: return "configuration error";
    case TerminationReason::EvaluationFailure: return "evaluation failure";
    case TerminationReason::MaxEvaluations: return "evaluation budget";
    }
    return "unknown";
}

inline GeneticTerminationReason genetic_termination_reason(const OptimizationTermination& termination) {
    return termination.reason;
}

} // namespace coilgun::optimization
