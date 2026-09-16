#pragma once

#include "coilgun/optimization/types.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>

namespace coilgun::optimization {

// The domain result predates the evaluation-budget criterion.  This local
// reason preserves that public type while allowing the optimizer to report the
// complete set of single-objective stop conditions in its own API.
enum class GeneticTerminationReason {
    None,
    MaxGenerations,
    MaxEvaluations,
    TargetReached,
    Converged,
    ConfigurationError,
    EvaluationFailure,
};

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

inline const char* to_string(GeneticTerminationReason reason) {
    switch (reason) {
    case GeneticTerminationReason::None: return "none";
    case GeneticTerminationReason::MaxGenerations: return "maximum generations";
    case GeneticTerminationReason::MaxEvaluations: return "evaluation budget";
    case GeneticTerminationReason::TargetReached: return "target reached";
    case GeneticTerminationReason::Converged: return "no improvement";
    case GeneticTerminationReason::ConfigurationError: return "configuration error";
    case GeneticTerminationReason::EvaluationFailure: return "evaluation failure";
    }
    return "unknown";
}

inline GeneticTerminationReason genetic_termination_reason(const OptimizationTermination& termination) {
    if (termination.reason == TerminationReason::MaxGenerations &&
        termination.message.find("evaluation budget") != std::string::npos)
        return GeneticTerminationReason::MaxEvaluations;
    switch (termination.reason) {
    case TerminationReason::MaxGenerations: return GeneticTerminationReason::MaxGenerations;
    case TerminationReason::TargetReached: return GeneticTerminationReason::TargetReached;
    case TerminationReason::Converged: return GeneticTerminationReason::Converged;
    case TerminationReason::ConfigurationError: return GeneticTerminationReason::ConfigurationError;
    case TerminationReason::EvaluationFailure: return GeneticTerminationReason::EvaluationFailure;
    case TerminationReason::Cancelled: return GeneticTerminationReason::None;
    case TerminationReason::None: return GeneticTerminationReason::None;
    }
    return GeneticTerminationReason::None;
}

} // namespace coilgun::optimization
