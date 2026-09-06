#include <doctest/doctest.h>

#include "coilgun/optimization/genetic_optimizer.hpp"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

using namespace coilgun::optimization;

namespace {
VariableSchema one_variable_schema() {
    return VariableSchema({VariableSpec::continuous("x", 0.0, 10.0)});
}

class ScoreEvaluator final : public BatchEvaluator {
public:
    explicit ScoreEvaluator(bool maximize = true, bool constant = false)
        : maximize_(maximize), constant_(constant) {}

    std::size_t calls = 0;
    std::size_t candidates = 0;

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        ++calls;
        candidates += values.size();
        std::vector<EvaluationResult> results;
        results.reserve(values.size());
        for (const auto& value : values) {
            auto result = EvaluationResult::success();
            const double x = value.values.front();
            result.objectives.push_back({"score", constant_ ? 1.0 : x, maximize_});
            results.push_back(std::move(result));
        }
        return results;
    }

private:
    bool maximize_;
    bool constant_;
};

class TwoObjectiveEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        std::vector<EvaluationResult> results;
        for (const auto& value : values) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"a", value.values.front(), true});
            result.objectives.push_back({"b", value.values.front(), false});
            results.push_back(std::move(result));
        }
        return results;
    }
};

class MixedDirectionEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        std::vector<EvaluationResult> results;
        for (std::size_t i = 0; i < values.size(); ++i) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"score", values[i].values.front(), i == 0});
            results.push_back(std::move(result));
        }
        return results;
    }
};

class FailedEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        return std::vector<EvaluationResult>(values.size(),
            EvaluationResult::failed("failed", "synthetic failure"));
    }
};

class TinyImprovementEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        const double score = calls++ == 0 ? 1.0 : 1.01;
        std::vector<EvaluationResult> results;
        results.reserve(values.size());
        for (const auto& value : values) {
            (void)value;
            auto result = EvaluationResult::success();
            result.objectives.push_back({"score", score, true});
            results.push_back(std::move(result));
        }
        return results;
    }

    std::size_t calls = 0;
};

class EliteDropEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        const bool first_generation = calls++ == 0;
        std::vector<EvaluationResult> results;
        results.reserve(values.size());
        for (const auto& value : values) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"score", first_generation ? 1.0 + value.values.front() : 0.0, true});
            results.push_back(std::move(result));
        }
        return results;
    }

    std::size_t calls = 0;
};

class ElitePropagationEvaluator final : public BatchEvaluator {
public:
    std::vector<double> first_generation_values;
    double elite_value = 0.0;
    std::size_t elite_occurrences_in_second_generation = 0;
    std::size_t calls = 0;

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        if (calls++ == 0) {
            first_generation_values.reserve(values.size());
            for (const auto& value : values) first_generation_values.push_back(value.values.front());
            elite_value = *std::max_element(first_generation_values.begin(), first_generation_values.end());
        } else {
            elite_occurrences_in_second_generation = static_cast<std::size_t>(std::count_if(
                values.begin(), values.end(), [&](const CandidateVariables& value) {
                    return value.values.front() == elite_value;
                }));
        }

        std::vector<EvaluationResult> results;
        results.reserve(values.size());
        for (const auto& value : values) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"score", value.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    }
};

class FeasibilityFirstEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        const bool first_generation = calls++ == 0;
        std::vector<EvaluationResult> results;
        results.reserve(values.size());
        const ConstraintDefinition limit{"limit", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                        0.0, 0.0, 1.0, 0};
        for (const auto& value : values) {
            (void)value;
            auto result = EvaluationResult::success();
            result.objectives.push_back({"score", first_generation ? 1.0 : 100.0, true});
            result.constraints.push_back(limit.evaluate(first_generation ? 0.0 : 1.0));
            results.push_back(std::move(result));
        }
        return results;
    }

    std::size_t calls = 0;
};

OptimizationConfig test_config() {
    OptimizationConfig config;
    config.population_size = 12;
    config.max_generations = 4;
    config.crossover_rate = 0.8;
    config.mutation_rate = 0.2;
    config.elite_count = 2;
    config.random_seed = 19;
    return config;
}
}

TEST_CASE("single objective optimizer honors maximize and minimize directions") {
    auto max_evaluator = std::make_shared<ScoreEvaluator>(true);
    auto max_config = test_config();
    max_config.max_generations = 2;
    const auto maximum = GeneticOptimizer(one_variable_schema(), max_evaluator, max_config).optimize();
    REQUIRE(maximum.termination.reason == TerminationReason::MaxGenerations);
    REQUIRE(maximum.pareto_front.size() == 1);
    REQUIRE(maximum.best_by_objective.count("score") == 1);
    const auto max_value = maximum.best_by_objective.at("score").objectives.front().value;

    auto min_evaluator = std::make_shared<ScoreEvaluator>(false);
    const auto minimum = GeneticOptimizer(one_variable_schema(), min_evaluator, max_config).optimize();
    REQUIRE(minimum.pareto_front.size() == 1);
    REQUIRE(minimum.best_by_objective.count("score") == 1);
    const auto min_value = minimum.best_by_objective.at("score").objectives.front().value;
    CHECK(max_value >= min_value);
    CHECK(maximum.best_by_objective.at("score").objectives.front().maximize);
    CHECK_FALSE(minimum.best_by_objective.at("score").objectives.front().maximize);
}

TEST_CASE("single objective optimizer is deterministic and preserves elites") {
    auto first_evaluator = std::make_shared<ScoreEvaluator>();
    auto first_config = test_config();
    first_config.max_generations = 3;
    const auto first = GeneticOptimizer(one_variable_schema(), first_evaluator, first_config).optimize();
    auto second_evaluator = std::make_shared<ScoreEvaluator>();
    const auto second = GeneticOptimizer(one_variable_schema(), second_evaluator, first_config).optimize();
    REQUIRE(first.best_by_objective.count("score") == 1);
    REQUIRE(second.best_by_objective.count("score") == 1);
    CHECK(first.best_by_objective.at("score").variables.values ==
          second.best_by_objective.at("score").variables.values);
    CHECK(first.statistics.evaluations == second.statistics.evaluations);
    CHECK(first.statistics.generations == second.statistics.generations);
    CHECK(first_evaluator->candidates == second_evaluator->candidates);
}

TEST_CASE("single objective optimizer supports target, evaluation, and no-improvement termination") {
    auto target_evaluator = std::make_shared<ScoreEvaluator>(true, true);
    auto config = test_config();
    config.max_generations = 20;
    TerminationConfig target;
    target.target_value = 1.0;
    const auto reached = GeneticOptimizer(one_variable_schema(), target_evaluator, config, target).optimize();
    CHECK(reached.termination.reason == TerminationReason::TargetReached);
    CHECK(reached.statistics.generations == 1);

    auto budget_evaluator = std::make_shared<ScoreEvaluator>();
    TerminationConfig budget;
    budget.max_evaluations = config.population_size;
    const auto budget_result = GeneticOptimizer(one_variable_schema(), budget_evaluator, config, budget).optimize();
    CHECK(budget_result.termination.reason == TerminationReason::MaxGenerations);
    CHECK(genetic_termination_reason(budget_result.termination) == GeneticTerminationReason::MaxEvaluations);
    CHECK(budget_result.statistics.evaluations <= budget.max_evaluations);

    auto stagnant_evaluator = std::make_shared<ScoreEvaluator>(true, true);
    TerminationConfig stagnant;
    stagnant.max_no_improvement_generations = 1;
    const auto converged = GeneticOptimizer(one_variable_schema(), stagnant_evaluator, config, stagnant).optimize();
    CHECK(converged.termination.reason == TerminationReason::Converged);
    CHECK(converged.statistics.generations == 2);
}

TEST_CASE("single objective optimizer rejects multiple objectives and all failed evaluations") {
    auto multiple = std::make_shared<TwoObjectiveEvaluator>();
    const auto invalid = GeneticOptimizer(one_variable_schema(), multiple, test_config()).optimize();
    CHECK(invalid.termination.reason == TerminationReason::ConfigurationError);
    CHECK(invalid.pareto_front.empty());

    auto mixed = std::make_shared<MixedDirectionEvaluator>();
    const auto mixed_result = GeneticOptimizer(one_variable_schema(), mixed, test_config()).optimize();
    CHECK(mixed_result.termination.reason == TerminationReason::ConfigurationError);

    auto failed = std::make_shared<FailedEvaluator>();
    const auto failure = GeneticOptimizer(one_variable_schema(), failed, test_config()).optimize();
    CHECK(failure.termination.reason == TerminationReason::EvaluationFailure);
    CHECK(failure.pareto_front.empty());
    CHECK(failure.statistics.failed_evaluations == test_config().population_size);
}

TEST_CASE("single objective optimizer converges when improvement is below tolerance") {
    auto evaluator = std::make_shared<TinyImprovementEvaluator>();
    auto config = test_config();
    config.max_generations = 5;
    TerminationConfig termination;
    termination.max_no_improvement_generations = 1;
    termination.improvement_tolerance = 0.1;

    const auto result = GeneticOptimizer(one_variable_schema(), evaluator, config, termination).optimize();

    CHECK(result.termination.reason == TerminationReason::Converged);
    CHECK(result.statistics.generations == 2);
    REQUIRE(result.best_by_objective.count("score") == 1);
    CHECK(result.best_by_objective.at("score").objectives.front().value == doctest::Approx(1.0));
}

TEST_CASE("single objective optimizer preserves the actual elite across generations") {
    auto evaluator = std::make_shared<EliteDropEvaluator>();
    auto config = test_config();
    config.max_generations = 2;
    config.crossover_rate = 0.0;
    config.mutation_rate = 0.0;

    const auto result = GeneticOptimizer(one_variable_schema(), evaluator, config).optimize();

    REQUIRE(result.best_by_objective.count("score") == 1);
    CHECK(result.best_by_objective.at("score").objectives.front().value > 1.0);
    CHECK(result.termination.reason == TerminationReason::MaxGenerations);
}

TEST_CASE("single objective elite propagation inserts the elite into the next generation") {
    auto evaluator = std::make_shared<ElitePropagationEvaluator>();
    auto config = test_config();
    config.population_size = 8;
    config.max_generations = 2;
    config.elite_count = 1;
    config.crossover_rate = 0.0;
    config.mutation_rate = 1.0;
    config.random_seed = 31;

    const auto result = GeneticOptimizer(one_variable_schema(), evaluator, config).optimize();

    REQUIRE(result.termination.reason == TerminationReason::MaxGenerations);
    REQUIRE(evaluator->first_generation_values.size() == config.population_size);
    CHECK(evaluator->elite_occurrences_in_second_generation == config.elite_count);
}

TEST_CASE("single objective optimizer keeps feasible incumbent ahead of infeasible objective gain") {
    auto evaluator = std::make_shared<FeasibilityFirstEvaluator>();
    auto config = test_config();
    config.max_generations = 3;
    config.crossover_rate = 0.0;
    config.mutation_rate = 0.0;
    TerminationConfig termination;
    termination.max_no_improvement_generations = 1;

    const auto result = GeneticOptimizer(one_variable_schema(), evaluator, config, termination).optimize();

    CHECK(result.termination.reason == TerminationReason::Converged);
    CHECK(result.statistics.generations == 2);
    REQUIRE(result.best_by_objective.count("score") == 1);
    const auto& best = result.best_by_objective.at("score");
    CHECK(best.objectives.front().value == doctest::Approx(1.0));
    CHECK(is_feasible(best.constraints));
}
