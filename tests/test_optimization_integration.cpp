#include <doctest/doctest.h>

#include "coilgun/optimization/coilgun_problem.hpp"
#include "coilgun/optimization/genetic_optimizer.hpp"
#include "coilgun/physics/constants.hpp"

#include <cmath>
#include <memory>
#include <stdexcept>

using namespace coilgun::optimization;
using coilgun::components::Armature;
using coilgun::physics::ALUMINUM;
using coilgun::physics::COPPER;
using coilgun::simulation::TriggerConfig;
using coilgun::simulation::TriggerMode;

namespace {
// Fixed-seed baseline after the post-RNG stream implementation on this branch.
constexpr double kCurrentFullVelocityBaseline = 0.0096453039804834419;
constexpr double kFullVelocityRegressionAbsoluteTolerance = 5e-8;
constexpr double kFullVelocityRegressionRelativeTolerance = 1e-6;

CoilgunOptimizationProblem::Config workflow_config() {
    CoilgunOptimizationProblem::Config config;
    config.coils.emplace_back(0.005, 0.010, 0.010, 12,
                              COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    config.armature = Armature(0.002, 0.008, 0.010,
                               ALUMINUM.resistivity_ref, ALUMINUM.density,
                               0.0, 0.005, 1, 1, 0.015);
    config.excitations = {{500.0, 500e-6, true}};
    config.triggers.clear();
    config.dt = 1e-6;
    config.termination.max_steps = 8;
    config.termination.enable_velocity_check = false;
    config.objective_id = "muzzle_velocity";
    config.constraints.push_back({"velocity_floor", CoilgunMetric::TerminalVelocity,
        ConstraintDefinition{"velocity_floor", ConstraintKind::Hard,
            ConstraintRelation::GreaterEqual, 0.0095, 0.0, 0.001}});
    return config;
}

VariableSchema workflow_schema() {
    return VariableSchema({VariableSpec::continuous("voltage", 450.0, 550.0)});
}

OptimizationConfig workflow_optimization() {
    OptimizationConfig config;
    config.population_size = 4;
    config.max_generations = 2;
    config.elite_count = 1;
    config.crossover_rate = 0.8;
    config.mutation_rate = 0.2;
    config.random_seed = 20260908;
    return config;
}
}

TEST_CASE("optimization workflow is reproducible, constrained, and reference-checkable") {
    const auto schema = workflow_schema();
    auto config = workflow_config();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    CoilgunOptimizationProblem first_problem(schema, config);
    CoilgunOptimizationProblem second_problem(schema, config);

    const auto optimization = workflow_optimization();
    const auto first = GeneticOptimizer(schema, first_problem, optimization).run();
    const auto second = GeneticOptimizer(schema, second_problem, optimization).run();

    REQUIRE(first.termination.reason == TerminationReason::MaxGenerations);
    REQUIRE(first.pareto_front.size() == 1);
    REQUIRE(first.best_by_objective.count("muzzle_velocity") == 1);
    const auto& best = first.best_by_objective.at("muzzle_velocity");
    REQUIRE(best.evaluation_status == EvaluationStatus::Success);
    REQUIRE(is_feasible(best.constraints));
    REQUIRE(best.constraints.front().value >= 0.0095);
    CHECK(best.constraints.front().violation == doctest::Approx(0.0));
    REQUIRE(std::isfinite(best.objectives.front().value));
    REQUIRE(best.variables.values.size() == 1);
    const double full_velocity_tolerance = kFullVelocityRegressionAbsoluteTolerance +
        kFullVelocityRegressionRelativeTolerance * std::abs(kCurrentFullVelocityBaseline);
    INFO("current Full velocity baseline = ", kCurrentFullVelocityBaseline,
         ", tolerance = ", full_velocity_tolerance,
         ", measured = ", best.objectives.front().value);
    CHECK(std::abs(best.objectives.front().value - kCurrentFullVelocityBaseline) <=
          full_velocity_tolerance);
    REQUIRE(second.best_by_objective.count("muzzle_velocity") == 1);
    CHECK(best.variables.values == second.best_by_objective.at("muzzle_velocity").variables.values);
    CHECK(best.objectives.front().value ==
          second.best_by_objective.at("muzzle_velocity").objectives.front().value);
    CHECK(first.statistics.evaluations == second.statistics.evaluations);
    CHECK(first.statistics.seed == 20260908);
    CHECK(first.statistics.successful_evaluations == first.statistics.evaluations);
    CHECK(first.statistics.failed_evaluations == 0);
    CHECK(first.statistics.cache_hits == 0);
    CHECK(first.statistics.gpu_fallbacks == 0);
    CHECK(first.statistics.gpu_requested_evaluations == 0);
    CHECK(first.statistics.gpu_executed_evaluations == 0);
    CHECK(first.statistics.gpu_successful_evaluations == 0);
    CHECK(first.statistics.gpu_failed_evaluations == 0);
    CHECK(first.statistics.cpu_fallback_evaluations == 0);
    CHECK(first.statistics.gpu_batches == 0);
    CHECK(first.statistics.gpu_failed_batches == 0);
    CHECK(first.statistics.gpu_transfer_seconds == doctest::Approx(0.0));
    CHECK(first.statistics.gpu_kernel_seconds == doctest::Approx(0.0));
    CHECK(first.statistics.gpu_elapsed_seconds == doctest::Approx(0.0));
    CHECK(first.statistics.skipped_due_to_budget == 0);
    CHECK(first.statistics.generations == 2);
    CHECK(std::isfinite(first.statistics.elapsed_seconds));
    CHECK(first.statistics.elapsed_seconds >= 0.0);

    auto reference_config = config;
    reference_config.optimization_level = coilgun::simulation::OptimizationLevel::Reference;
    CoilgunOptimizationProblem reference_problem(schema, std::move(reference_config));
    const auto reference = reference_problem.evaluate(best.variables);
    REQUIRE(reference.status == EvaluationStatus::Success);
    REQUIRE(reference.objectives.size() == 1);
    CHECK(std::isfinite(reference.objectives.front().value));
    const double reference_error = std::abs(reference.objectives.front().value -
                                            best.objectives.front().value);
    REQUIRE_MESSAGE(reference_error > 1e-12,
                    "validation workload must exercise distinct Full and Reference paths");
    const double reference_tolerance = kFullVelocityRegressionAbsoluteTolerance +
        kFullVelocityRegressionRelativeTolerance * std::abs(reference.objectives.front().value);
    CHECK(reference_error <= reference_tolerance);
}

TEST_CASE("optimization workflow isolates failed batch candidates and records adapter fallback boundary") {
    const auto schema = workflow_schema();
    auto config = workflow_config();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    auto problem = std::make_shared<CoilgunOptimizationProblem>(schema, std::move(config));
    std::size_t callback_calls = 0;
    problem->set_gpu_batch_evaluator([&callback_calls](const std::vector<CandidateVariables>& candidates,
                                                      const EvaluationContext&) {
        ++callback_calls;
        std::vector<EvaluationResult> output;
        output.reserve(candidates.size());
        for (const auto& candidate : candidates) {
            auto result = EvaluationResult::failed("candidate_failure", "synthetic candidate failure");
            if (candidate.values.front() > 500.0) {
                result = EvaluationResult::success();
                result.objectives.push_back({"muzzle_velocity", 0.0, true});
                result.constraints.push_back({"velocity_floor", ConstraintKind::Hard,
                    ConstraintRelation::GreaterEqual, 0.01, 0.0095, 0.0,
                    0.0, 0.0, true, 0});
            }
            output.push_back(std::move(result));
        }
        return output;
    });

    std::shared_ptr<BatchEvaluator> problem_evaluator = problem;
    StatisticsBatchEvaluator tracked(problem_evaluator);
    const auto batch = tracked.evaluate_batch({CandidateVariables{{490.0}}, CandidateVariables{{510.0}}},
                                              EvaluationContext{77, true});
    REQUIRE(batch.size() == 2);
    CHECK(batch[0].status == EvaluationStatus::Failed);
    CHECK(batch[1].status == EvaluationStatus::Success);
    CHECK_FALSE(problem->last_batch_used_fallback());
    CHECK(callback_calls >= 1);
    CHECK(tracked.statistics().failed_evaluations == 1);
    CHECK(tracked.statistics().successful_evaluations == 1);

    problem->set_gpu_batch_evaluator([&callback_calls](const std::vector<CandidateVariables>&,
                                                      const EvaluationContext&)
                                        -> std::vector<EvaluationResult> {
        ++callback_calls;
        throw std::runtime_error("synthetic GPU unavailable");
    });
    const auto fallback = problem->evaluate_batch({CandidateVariables{{490.0}}, CandidateVariables{{510.0}}},
                                                 EvaluationContext{77, true});
    REQUIRE(fallback.size() == 2);
    CHECK(fallback[0].status == EvaluationStatus::Success);
    CHECK(fallback[1].status == EvaluationStatus::Success);
    CHECK(problem->last_batch_used_fallback());

    auto cache = std::make_shared<InMemoryEvaluationCache>();
    CachedBatchEvaluator cached(problem_evaluator, cache);
    const auto cached_batch = cached.evaluate_batch({CandidateVariables{{510.0}}, CandidateVariables{{510.0}}},
                                                    EvaluationContext{77, false});
    REQUIRE(cached_batch.size() == 2);
    CHECK(cached.statistics().evaluations == 1);
    CHECK(cached.statistics().cache_hits == 0);
    CHECK(problem->last_batch_used_fallback());
    const auto second = cached.evaluate_batch({CandidateVariables{{510.0}}}, EvaluationContext{77, false});
    REQUIRE(second.size() == 1);
    CHECK(cached.statistics().cache_hits == 1);
    CHECK(cached.statistics().fallbacks == 0);
}
