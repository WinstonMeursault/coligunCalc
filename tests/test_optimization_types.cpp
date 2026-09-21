#include <doctest/doctest.h>

#include "coilgun/optimization/config.hpp"
#include "coilgun/optimization/problem.hpp"
#include "coilgun/optimization/termination.hpp"
#include "coilgun/optimization/types.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

using namespace coilgun::optimization;

static_assert(std::is_same_v<GeneticTerminationReason, TerminationReason>);

TEST_CASE("optimization domain types have useful empty defaults") {
    CandidateVariables variables;
    Candidate candidate;
    EvaluationResult result;
    OptimizationStatistics statistics;
    OptimizationTermination termination;

    CHECK(variables.values.empty());
    CHECK(candidate.id == 0);
    CHECK(candidate.evaluation_status == EvaluationStatus::Unevaluated);
    CHECK(result.status == EvaluationStatus::Unevaluated);
    CHECK(statistics.evaluations == 0);
    CHECK(termination.reason == TerminationReason::None);
}

TEST_CASE("candidate domain values are copyable and retain stable identity") {
    Candidate original{42, CandidateVariables{{1.0, 2.0}}};
    original.objectives.push_back(ObjectiveValue{"velocity", 12.5, true});
    original.constraints.push_back(ConstraintReport{"mass", ConstraintKind::Hard,
        ConstraintRelation::GreaterEqual, 0.2, 0.1, 0.0, 0.0, true});
    Candidate copy = original;

    CHECK(copy.id == 42);
    CHECK(copy.variables.values == std::vector<double>{1.0, 2.0});
    CHECK(copy.objectives.front().id == "velocity");
    CHECK(copy.constraints.front().satisfied);
}

TEST_CASE("optimization problem provides an abstract evaluation boundary") {
    struct Problem final : OptimizationProblem {
        EvaluationResult evaluate(const CandidateVariables& variables) const override {
            EvaluationResult result;
            result.status = EvaluationStatus::Success;
            result.objectives.push_back(ObjectiveValue{"sum", variables.values[0] + variables.values[1], false});
            return result;
        }
    } problem;

    const auto result = problem.evaluate(CandidateVariables{{2.0, 3.0}});
    REQUIRE(result.status == EvaluationStatus::Success);
    CHECK(result.objectives.front().value == doctest::Approx(5.0));
}

TEST_CASE("optimization config has deterministic defaults and validates required settings") {
    const auto defaults = OptimizationConfig::defaults();
    CHECK(defaults.population_size > 0);
    CHECK(defaults.max_generations > 0);
    CHECK_NOTHROW(defaults.validate());

    auto invalid = defaults;
    invalid.population_size = 0;
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);
    invalid = defaults;
    invalid.crossover_rate = 1.5;
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);
}

TEST_CASE("termination reasons have complete structured names") {
    const std::vector<std::pair<TerminationReason, std::string>> reasons = {
        {TerminationReason::None, "none"},
        {TerminationReason::MaxGenerations, "maximum generations"},
        {TerminationReason::TargetReached, "target reached"},
        {TerminationReason::Converged, "no improvement"},
        {TerminationReason::Cancelled, "cancelled"},
        {TerminationReason::ConfigurationError, "configuration error"},
        {TerminationReason::EvaluationFailure, "evaluation failure"},
        {TerminationReason::MaxEvaluations, "evaluation budget"},
    };
    for (const auto& [reason, name] : reasons)
        CHECK(std::string(to_string(reason)) == name);
}

TEST_CASE("termination compatibility helper ignores diagnostic message text") {
    const OptimizationTermination misleading{
        TerminationReason::MaxGenerations, "evaluation budget exhausted", 3};
    CHECK(genetic_termination_reason(misleading) == TerminationReason::MaxGenerations);

    const OptimizationTermination localized{
        TerminationReason::MaxEvaluations, "预算已耗尽", 3};
    CHECK(genetic_termination_reason(localized) == TerminationReason::MaxEvaluations);

    const OptimizationTermination cancelled{
        TerminationReason::Cancelled, "cancelled by operator", 3};
    CHECK(genetic_termination_reason(cancelled) == TerminationReason::Cancelled);
}
