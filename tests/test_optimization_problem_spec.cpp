#include <doctest/doctest.h>

#include "coilgun/optimization/coilgun_problem.hpp"
#include "coilgun/optimization/genetic_optimizer.hpp"
#include "coilgun/optimization/problem.hpp"

#include <concepts>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace coilgun::optimization;

namespace {

template <typename T>
concept HasConstBatchEvaluate = requires(
    const T& evaluator,
    const std::vector<CandidateVariables>& candidates,
    const EvaluationContext& context) {
    { evaluator.evaluate_batch(candidates, context) } -> std::same_as<std::vector<EvaluationResult>>;
};

static_assert(!HasConstBatchEvaluate<BatchEvaluator>);
static_assert(HasConstBatchEvaluate<CoilgunOptimizationProblem>);

VariableSchema schema() {
    return VariableSchema({VariableSpec::continuous("x", 0.0, 1.0),
                           VariableSpec::integer("n", 1, 3)});
}

ObjectiveDefinition score() { return {"score", true, 2.0}; }

ConstraintDefinition limit() {
    return {"limit", ConstraintKind::Hard, ConstraintRelation::LessEqual,
            0.0, 1.0, 1.0, 0};
}

OptimizationConfig config() {
    OptimizationConfig result;
    result.population_size = 4;
    result.max_generations = 1;
    result.crossover_rate = 0.0;
    result.mutation_rate = 0.0;
    result.random_seed = 11;
    return result;
}

class LegacyProblem final : public OptimizationProblem {
public:
    mutable std::size_t evaluations = 0;
    EvaluationResult evaluate(const CandidateVariables&) const override {
        ++evaluations;
        auto result = EvaluationResult::success();
        result.objectives.push_back({"score", 1.0, true});
        return result;
    }
};

class SpecProblem final : public OptimizationProblem {
public:
    explicit SpecProblem(ProblemSpec spec) : OptimizationProblem(std::move(spec)) {}
    EvaluationResult evaluate(const CandidateVariables& value) const override {
        auto result = EvaluationResult::success();
        result.objectives.push_back({"score", value.values.front(), true});
        result.constraints.push_back(limit().evaluate(value.values.front()));
        return result;
    }
};

class BatchSpecProblem final : public OptimizationProblem, public BatchEvaluator {
public:
    explicit BatchSpecProblem(ProblemSpec spec) : OptimizationProblem(std::move(spec)) {}
    mutable std::size_t serial_evaluations = 0;
    std::size_t batch_calls = 0;
    EvaluationResult evaluate(const CandidateVariables&) const override {
        ++serial_evaluations;
        return EvaluationResult::failed("unexpected_serial", "serial path used");
    }
    std::vector<EvaluationResult> evaluate_batch(
        const std::vector<CandidateVariables>& values, const EvaluationContext&) override {
        ++batch_calls;
        std::vector<EvaluationResult> results;
        for (const auto& value : values) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"score", value.values.front(), true});
            result.constraints.push_back(limit().evaluate(value.values.front()));
            results.push_back(std::move(result));
        }
        return results;
    }
};

class ConstDualSpecProblem final : public OptimizationProblem, public BatchEvaluator {
public:
    explicit ConstDualSpecProblem(ProblemSpec spec) : OptimizationProblem(std::move(spec)) {}
    mutable std::size_t serial_evaluations = 0;
    mutable std::size_t legacy_batch_calls = 0;

    EvaluationResult evaluate(const CandidateVariables& value) const override {
        ++serial_evaluations;
        auto result = EvaluationResult::success();
        result.objectives.push_back({"score", value.values.front(), true});
        result.constraints.push_back(limit().evaluate(value.values.front()));
        return result;
    }

    std::vector<EvaluationResult> evaluate_batch(
        const std::vector<CandidateVariables>& values, const EvaluationContext&) override {
        ++legacy_batch_calls;
        return std::vector<EvaluationResult>(values.size(),
            EvaluationResult::failed("unexpected_legacy_batch", "legacy non-const batch path used"));
    }
};

} // namespace

TEST_CASE("ProblemSpec validates and exposes only read-only definitions") {
    ProblemSpec spec(schema(), {score()}, {limit()});
    CHECK(spec.schema().size() == 2);
    CHECK(spec.objectives().front().id == "score");
    CHECK(spec.constraints().front().id == "limit");
    CHECK(spec.repair(CandidateVariables{{2.0, 2.8}}).values ==
          std::vector<double>{1.0, 3.0});
    CHECK_THROWS_AS(ProblemSpec(schema(), {ObjectiveDefinition{"", true, 1.0}}),
                    std::invalid_argument);
    CHECK_THROWS_AS(ProblemSpec(schema(), {score(), score()}), std::invalid_argument);
    CHECK_THROWS_AS(ProblemSpec(schema(), {score()}, {limit(), limit()}),
                    std::invalid_argument);
    CHECK_THROWS_AS(ProblemSpec(schema(), {score()}, {}, RepairPolicy{}),
                    std::invalid_argument);
}

TEST_CASE("ProblemSpec custom repair is invoked and schema remains canonical authority") {
    std::size_t calls = 0;
    ProblemSpec spec(schema(), {score()}, {},
        [&](const VariableSchema&, const CandidateVariables&) {
            ++calls;
            return CandidateVariables{{-100.0, 100.25}};
        });
    CHECK(spec.repair(CandidateVariables{{0.4, 2.0}}).values ==
          std::vector<double>{0.0, 3.0});
    CHECK(calls == 1);
}

TEST_CASE("spec-aware optimizer applies custom repair to initial candidates and offspring") {
    std::size_t calls = 0;
    ProblemSpec spec(schema(), {score()}, {},
        [&](const VariableSchema& owner, const CandidateVariables& candidate) {
            ++calls;
            return owner.repair(candidate);
        });
    class Evaluator final : public BatchEvaluator {
    public:
        std::vector<EvaluationResult> evaluate_batch(
            const std::vector<CandidateVariables>& values, const EvaluationContext&) override {
            std::vector<EvaluationResult> results;
            for (const auto& value : values) {
                auto result = EvaluationResult::success();
                result.objectives.push_back({"score", value.values.front(), true});
                results.push_back(std::move(result));
            }
            return results;
        }
    } evaluator;
    auto optimizer_config = config();
    optimizer_config.max_generations = 2;
    GeneticOptimizer(spec, evaluator, optimizer_config).optimize();
    CHECK(calls == optimizer_config.population_size +
          (optimizer_config.population_size - optimizer_config.elite_count));
}

TEST_CASE("spec-aware problem optimizer needs no duplicate schema and rejects legacy problem") {
    SpecProblem problem(ProblemSpec(schema(), {score()}, {limit()}));
    const auto result = GeneticOptimizer(problem, config()).optimize();
    CHECK(result.termination.reason == TerminationReason::MaxGenerations);

    LegacyProblem legacy;
    CHECK_THROWS_WITH_AS(GeneticOptimizer(legacy, config()),
                         "optimization problem does not provide a ProblemSpec",
                         std::invalid_argument);
    CHECK(legacy.evaluations == 0);
}

TEST_CASE("spec-aware batch problem preserves batch dispatch") {
    BatchSpecProblem problem(ProblemSpec(schema(), {score()}, {limit()}));
    const auto result = GeneticOptimizer(problem, config()).optimize();
    CHECK(result.termination.reason == TerminationReason::MaxGenerations);
    CHECK(problem.batch_calls == 1);
    CHECK(problem.serial_evaluations == 0);
}

TEST_CASE("const problem-only optimization never dispatches legacy non-const batch") {
    const ConstDualSpecProblem problem(ProblemSpec(schema(), {score()}, {limit()}));
    const auto result = GeneticOptimizer(problem, config()).optimize();

    CHECK(result.termination.reason == TerminationReason::MaxGenerations);
    CHECK(problem.legacy_batch_calls == 0);
    CHECK(problem.serial_evaluations == config().population_size);
}

TEST_CASE("spec Auto NSGA-II stagnation rejects before repair or evaluator work") {
    std::size_t repair_calls = 0;
    ProblemSpec spec(schema(), {score(), {"other", true, 1.0}}, {},
        [&](const VariableSchema& owner, const CandidateVariables& candidate) {
            ++repair_calls;
            return owner.repair(candidate);
        });
    class Evaluator final : public BatchEvaluator {
    public:
        std::size_t calls = 0;
        std::vector<EvaluationResult> evaluate_batch(
            const std::vector<CandidateVariables>& values, const EvaluationContext&) override {
            ++calls;
            return std::vector<EvaluationResult>(values.size(), EvaluationResult::failed("unexpected", "unexpected"));
        }
    } evaluator;
    auto termination = TerminationConfig{};
    termination.max_no_improvement_generations = 2;

    const auto result = GeneticOptimizer(spec, evaluator, config(), termination).optimize();

    CHECK(result.termination.reason == TerminationReason::ConfigurationError);
    CHECK(result.statistics.generations == 0);
    CHECK(repair_calls == 0);
    CHECK(evaluator.calls == 0);
}

TEST_CASE("spec-aware optimizer structurally validates declared output schema") {
    class Mismatch final : public BatchEvaluator {
    public:
        std::vector<EvaluationResult> evaluate_batch(
            const std::vector<CandidateVariables>& values, const EvaluationContext&) override {
            std::vector<EvaluationResult> results;
            for (std::size_t i = 0; i < values.size(); ++i) {
                auto result = EvaluationResult::success();
                result.objectives.push_back({"other", 1.0, true});
                results.push_back(std::move(result));
            }
            return results;
        }
    } evaluator;
    ProblemSpec spec(schema(), {score()});
    const auto result = GeneticOptimizer(spec, evaluator, config()).optimize();
    CHECK(result.termination.reason == TerminationReason::ConfigurationError);
}
