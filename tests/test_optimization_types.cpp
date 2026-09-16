#include <doctest/doctest.h>

#include "coilgun/optimization/config.hpp"
#include "coilgun/optimization/problem.hpp"
#include "coilgun/optimization/types.hpp"

#include <stdexcept>

using namespace coilgun::optimization;

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
