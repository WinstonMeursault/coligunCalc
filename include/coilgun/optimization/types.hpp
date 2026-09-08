#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace coilgun::optimization {

class RepresentativeSelector;

using CandidateId = std::uint64_t;

struct CandidateVariables {
    std::vector<double> values;
    CandidateVariables() = default;
    explicit CandidateVariables(std::vector<double> values_) : values(std::move(values_)) {}
};

enum class EvaluationStatus { Unevaluated, Success, Invalid, Failed };
enum class DiagnosticSeverity { Info, Warning, Error };

struct Diagnostic {
    std::string code;
    std::string message;
    DiagnosticSeverity severity = DiagnosticSeverity::Error;
};

struct ObjectiveValue {
    std::string id;
    double value = 0.0;
    bool maximize = true;
};

enum class ConstraintKind { Hard, Soft };
enum class ConstraintRelation { Equal, LessEqual, GreaterEqual, InRange };

struct ConstraintReport {
    std::string id;
    ConstraintKind kind = ConstraintKind::Hard;
    ConstraintRelation relation = ConstraintRelation::LessEqual;
    double value = 0.0;
    double lower_bound = 0.0;
    double upper_bound = 0.0;
    double violation = 0.0;
    double normalized_violation = 0.0;
    bool satisfied = true;
    // Lower values represent higher precedence in Lexicographic comparisons.
    int priority = 0;
};

struct EvaluationResult {
    EvaluationStatus status = EvaluationStatus::Unevaluated;
    std::vector<ObjectiveValue> objectives;
    std::vector<ConstraintReport> constraints;
    std::vector<Diagnostic> diagnostics;
    std::unordered_map<std::string, std::string> metadata;

    static EvaluationResult success() {
        EvaluationResult result;
        result.status = EvaluationStatus::Success;
        return result;
    }
    static EvaluationResult invalid(std::string code, std::string message) {
        EvaluationResult result;
        result.status = EvaluationStatus::Invalid;
        result.diagnostics.push_back({std::move(code), std::move(message), DiagnosticSeverity::Error});
        return result;
    }
    static EvaluationResult failed(std::string code, std::string message) {
        EvaluationResult result;
        result.status = EvaluationStatus::Failed;
        result.diagnostics.push_back({std::move(code), std::move(message), DiagnosticSeverity::Error});
        return result;
    }
};

struct Candidate {
    CandidateId id = 0;
    CandidateVariables variables;
    std::vector<ObjectiveValue> objectives;
    std::vector<ConstraintReport> constraints;
    EvaluationStatus evaluation_status = EvaluationStatus::Unevaluated;
    std::vector<Diagnostic> diagnostics;
    std::unordered_map<std::string, std::string> metadata;
};

struct OptimizationStatistics {
    std::uint64_t evaluations = 0;
    std::uint64_t successful_evaluations = 0;
    std::uint64_t failed_evaluations = 0;
    std::uint64_t cache_hits = 0;
    std::uint64_t gpu_fallbacks = 0;
    std::uint64_t generations = 0;
    double elapsed_seconds = 0.0;
};

enum class TerminationReason { None, MaxGenerations, TargetReached, Converged, Cancelled,
                               ConfigurationError, EvaluationFailure };

struct OptimizationTermination {
    TerminationReason reason = TerminationReason::None;
    std::string message;
    std::uint64_t generation = 0;
};

struct OptimizationResult {
    std::vector<Candidate> pareto_front;
    std::unordered_map<std::string, Candidate> best_by_objective;
    OptimizationStatistics statistics;
    OptimizationTermination termination;

    std::optional<Candidate> select_representative(const RepresentativeSelector& selector) const;
};

} // namespace coilgun::optimization
