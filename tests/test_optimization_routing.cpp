#include <doctest/doctest.h>

#include "coilgun/optimization/genetic_optimizer.hpp"

#include <memory>
#include <limits>
#include <vector>

using namespace coilgun::optimization;

namespace {

VariableSchema schema() {
    return VariableSchema({VariableSpec::continuous("x", 0.0, 1.0)});
}

OptimizationConfig config(SelectionStrategy strategy = SelectionStrategy::Auto) {
    OptimizationConfig result;
    result.population_size = 12;
    result.max_generations = 3;
    result.crossover_rate = 0.0;
    result.mutation_rate = 0.0;
    result.elite_count = 2;
    result.random_seed = 7;
    result.strategy = strategy;
    return result;
}

class OneObjectiveEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
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

class TwoObjectiveEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        std::vector<EvaluationResult> results;
        results.reserve(values.size());
        for (const auto& value : values) {
            const double x = value.values.front();
            auto result = EvaluationResult::success();
            result.objectives.push_back({"left", x, true});
            result.objectives.push_back({"right", 1.0 - x, true});
            results.push_back(std::move(result));
        }
        return results;
    }
};

class ChangingObjectiveEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        const bool first_generation = calls++ == 0;
        std::vector<EvaluationResult> results;
        results.reserve(values.size());
        for (const auto& value : values) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"left", value.values.front(), true});
            if (first_generation) result.objectives.push_back({"right", 1.0 - value.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    }

    std::size_t calls = 0;
};

class SchemaMutationEvaluator final : public BatchEvaluator {
public:
    enum class Mutation { Id, Direction, NonFinite };

    explicit SchemaMutationEvaluator(Mutation mutation) : mutation_(mutation) {}

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        const bool second_generation = calls++ > 0;
        std::vector<EvaluationResult> results;
        results.reserve(values.size());
        for (const auto& value : values) {
            const double x = value.values.front();
            auto result = EvaluationResult::success();
            const bool mutate_id = second_generation && mutation_ == Mutation::Id;
            const bool mutate_direction = second_generation && mutation_ == Mutation::Direction;
            const bool mutate_value = second_generation && mutation_ == Mutation::NonFinite;
            result.objectives.push_back({mutate_id ? "changed" : "left", x, !mutate_direction});
            result.objectives.push_back({"right", mutate_value
                ? std::numeric_limits<double>::quiet_NaN() : 1.0 - x, true});
            results.push_back(std::move(result));
        }
        return results;
    }

private:
    Mutation mutation_;
    std::size_t calls = 0;
};

} // namespace

TEST_CASE("Auto routing is equivalent to explicit single-objective optimization") {
    auto auto_evaluator = std::make_shared<OneObjectiveEvaluator>();
    auto explicit_evaluator = std::make_shared<OneObjectiveEvaluator>();

    const auto automatic = GeneticOptimizer(schema(), auto_evaluator, config()).optimize();
    const auto explicit_single = GeneticOptimizer(
        schema(), explicit_evaluator, config(SelectionStrategy::SingleObjective)).optimize();

    REQUIRE(automatic.termination.reason == TerminationReason::MaxGenerations);
    REQUIRE(explicit_single.termination.reason == TerminationReason::MaxGenerations);
    REQUIRE(automatic.best_by_objective.count("score") == 1);
    REQUIRE(explicit_single.best_by_objective.count("score") == 1);
    CHECK(automatic.best_by_objective.at("score").variables.values ==
          explicit_single.best_by_objective.at("score").variables.values);
    CHECK(automatic.statistics.evaluations == explicit_single.statistics.evaluations);
}

TEST_CASE("Auto routing uses NSGA-II for multiple objectives") {
    auto evaluator = std::make_shared<TwoObjectiveEvaluator>();
    const auto result = GeneticOptimizer(schema(), evaluator, config()).optimize();

    REQUIRE(result.termination.reason == TerminationReason::MaxGenerations);
    CHECK(result.pareto_front.size() >= 2);
    CHECK(result.best_by_objective.empty());
    for (const auto& candidate : result.pareto_front)
        CHECK(candidate.objectives.size() == 2);
}

TEST_CASE("Explicit strategies reject incompatible objective counts") {
    auto multi_evaluator = std::make_shared<TwoObjectiveEvaluator>();
    const auto single_result = GeneticOptimizer(
        schema(), multi_evaluator, config(SelectionStrategy::SingleObjective)).optimize();
    CHECK(single_result.termination.reason == TerminationReason::ConfigurationError);

    auto single_evaluator = std::make_shared<OneObjectiveEvaluator>();
    const auto nsga_result = GeneticOptimizer(
        schema(), single_evaluator, config(SelectionStrategy::NSGA2)).optimize();
    CHECK(nsga_result.termination.reason == TerminationReason::ConfigurationError);
}

TEST_CASE("Auto routing freezes objective count for the whole run") {
    auto evaluator = std::make_shared<ChangingObjectiveEvaluator>();
    const auto result = GeneticOptimizer(schema(), evaluator, config()).optimize();

    CHECK(result.termination.reason == TerminationReason::ConfigurationError);
    CHECK(result.termination.message.find("objective") != std::string::npos);
}

TEST_CASE("Auto routing freezes objective IDs, directions, and finite values") {
    for (const auto mutation : {SchemaMutationEvaluator::Mutation::Id,
                                SchemaMutationEvaluator::Mutation::Direction,
                                SchemaMutationEvaluator::Mutation::NonFinite}) {
        auto evaluator = std::make_shared<SchemaMutationEvaluator>(mutation);
        const auto result = GeneticOptimizer(schema(), evaluator, config()).optimize();
        CHECK(result.termination.reason == TerminationReason::ConfigurationError);
        CHECK(result.termination.message.find("objective") != std::string::npos);
    }
}
