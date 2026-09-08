#include <doctest/doctest.h>

#include "coilgun/coilgun.hpp"

#include <memory>

namespace {
class SmokeProblem final : public coilgun::optimization::OptimizationProblem {
public:
    coilgun::optimization::EvaluationResult evaluate(
        const coilgun::optimization::CandidateVariables&) const override {
        auto result = coilgun::optimization::EvaluationResult::success();
        result.objectives.push_back({"score", 1.0, true});
        return result;
    }
};
}

TEST_CASE("umbrella header exposes optimization problem and result APIs") {
    using namespace coilgun::optimization;

    SmokeProblem problem;
    VariableSchema schema({VariableSpec::continuous("x", 0.0, 1.0)});
    auto config = OptimizationConfig::defaults();
    config.population_size = 4;
    config.max_generations = 1;
    config.random_seed = 7;

    const auto result = GeneticOptimizer(schema, problem, config).run();
    REQUIRE(result.termination.reason == TerminationReason::MaxGenerations);
    REQUIRE_FALSE(result.pareto_front.empty());

    const auto selected = result.select_representative(MaxObjective{"score"});
    REQUIRE(selected.has_value());
    CHECK(selected->evaluation_status == EvaluationStatus::Success);
}
