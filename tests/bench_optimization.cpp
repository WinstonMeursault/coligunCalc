#include "coilgun/optimization/coilgun_problem.hpp"
#include "coilgun/optimization/genetic_optimizer.hpp"
#include "coilgun/physics/constants.hpp"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace coilgun::optimization;
using coilgun::components::Armature;
using coilgun::physics::ALUMINUM;
using coilgun::physics::COPPER;

namespace {
CoilgunOptimizationProblem::Config workload_config() {
    CoilgunOptimizationProblem::Config config;
    config.coils.emplace_back(0.005, 0.010, 0.010, 12,
                              COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    config.armature = Armature(0.002, 0.008, 0.010,
                               ALUMINUM.resistivity_ref, ALUMINUM.density,
                               0.0, 0.005, 1, 1, 0.015);
    config.excitations = {{500.0, 500e-6, true}};
    config.dt = 1e-6;
    config.termination.max_steps = 8;
    config.termination.enable_velocity_check = false;
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    config.constraints.push_back({"velocity_floor", CoilgunMetric::TerminalVelocity,
        ConstraintDefinition{"velocity_floor", ConstraintKind::Hard,
            ConstraintRelation::GreaterEqual, 0.0095, 0.0, 0.001}});
    return config;
}
}

int main() {
    using clock = std::chrono::steady_clock;
    const auto setup_start = clock::now();
    const VariableSchema schema({VariableSpec::continuous("voltage", 450.0, 550.0)});
    CoilgunOptimizationProblem problem(schema, workload_config());
    OptimizationConfig config;
    config.population_size = 4;
    config.max_generations = 2;
    config.elite_count = 1;
    config.random_seed = 20260908;
    const auto setup_seconds = std::chrono::duration<double>(clock::now() - setup_start).count();

    const std::vector<CandidateVariables> batch{
        CandidateVariables{{500.0}}, CandidateVariables{{510.0}},
        CandidateVariables{{520.0}}, CandidateVariables{{530.0}}};
    const auto first_start = clock::now();
    const auto first = problem.evaluate_batch(batch, EvaluationContext{config.random_seed, false});
    const auto first_seconds = std::chrono::duration<double>(clock::now() - first_start).count();

    constexpr int warmup_runs = 2;
    for (int i = 0; i < warmup_runs; ++i)
        (void)problem.evaluate_batch(batch, EvaluationContext{config.random_seed, false});
    const auto steady_start = clock::now();
    constexpr int steady_runs = 5;
    for (int i = 0; i < steady_runs; ++i)
        (void)problem.evaluate_batch(batch, EvaluationContext{config.random_seed, false});
    const auto steady_seconds = std::chrono::duration<double>(clock::now() - steady_start).count();

    auto cache = std::make_shared<InMemoryEvaluationCache>();
    auto evaluator = std::shared_ptr<BatchEvaluator>(&problem, [](BatchEvaluator*) {});
    auto cached = std::make_shared<CachedBatchEvaluator>(evaluator, cache);
    const auto result = GeneticOptimizer(schema, cached, config).run();
    const auto& best = result.best_by_objective.at("muzzle_velocity");
    auto reference_config = problem.config();
    reference_config.optimization_level = coilgun::simulation::OptimizationLevel::Reference;
    CoilgunOptimizationProblem reference_problem(schema, std::move(reference_config));
    const auto reference = reference_problem.evaluate(best.variables);
    const double reference_error = std::abs(reference.objectives.front().value - best.objectives.front().value);

    std::size_t callback_calls = 0;
    problem.set_gpu_batch_evaluator([&callback_calls](const std::vector<CandidateVariables>& candidates,
                                                      const EvaluationContext&) {
        ++callback_calls;
        std::vector<EvaluationResult> output;
        output.reserve(candidates.size());
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            if (i == 0) {
                output.push_back(EvaluationResult::failed("candidate_failure", "synthetic failure"));
            } else {
                auto success = EvaluationResult::success();
                success.objectives.push_back({"muzzle_velocity", 0.0, true});
                output.push_back(std::move(success));
            }
        }
        return output;
    });
    const auto isolated = problem.evaluate_batch(batch, EvaluationContext{config.random_seed, false});
    const auto isolated_failures = std::count_if(isolated.begin(), isolated.end(), [](const auto& value) {
        return value.status != EvaluationStatus::Success;
    });

    problem.set_gpu_batch_evaluator([&callback_calls](const std::vector<CandidateVariables>&,
                                                      const EvaluationContext&)
                                        -> std::vector<EvaluationResult> {
        ++callback_calls;
        throw std::runtime_error("synthetic GPU unavailable");
    });
    const auto fallback_start = clock::now();
    const auto fallback_batch = problem.evaluate_batch(batch, EvaluationContext{config.random_seed, true});
    const auto fallback_seconds = std::chrono::duration<double>(clock::now() - fallback_start).count();
    const bool fallback_observed = problem.last_batch_used_fallback();

    std::cout << std::setprecision(10)
              << "source_revision=" << OPTIMIZATION_BENCH_SOURCE_REVISION << '\n'
              << "worktree_state=" << OPTIMIZATION_BENCH_WORKTREE_STATE << '\n'
              << "preset=" << OPTIMIZATION_BENCH_PRESET << '\n'
              << "seed=" << config.random_seed << '\n'
              << "setup_seconds=" << setup_seconds << '\n'
              << "first_batch_seconds=" << first_seconds << '\n'
              << "warmup_runs=" << warmup_runs << '\n'
              << "steady_state_runs=" << steady_runs << '\n'
              << "steady_state_seconds=" << steady_seconds << '\n'
              << "steady_state_per_batch_seconds=" << steady_seconds / steady_runs << '\n'
              << "evaluations=" << result.statistics.evaluations << '\n'
              << "cache_hits=" << result.statistics.cache_hits << '\n'
              << "optimizer_failed_evaluations=" << result.statistics.failed_evaluations << '\n'
              << "failed_evaluations=" << isolated_failures << '\n'
              << "fallback_count=" << static_cast<unsigned>(fallback_observed) << '\n'
              << "terminal_velocity=" << best.objectives.front().value << '\n'
              << "feasible=" << is_feasible(best.constraints) << '\n'
              << "reference_terminal_velocity=" << reference.objectives.front().value << '\n'
              << "reference_recheck_error=" << reference_error << '\n'
              << "gpu_batch_backend=unavailable (no concrete CUDA optimizer backend)" << '\n'
              << "gpu_callback_calls=" << callback_calls << '\n'
              << "gpu_callback_isolated_failures=" << isolated_failures << '\n'
              << "gpu_callback_isolated_successes=" << isolated.size() - isolated_failures << '\n'
              << "gpu_callback_fallback_observed=" << fallback_observed << '\n'
              << "gpu_callback_fallback_seconds=" << fallback_seconds << '\n'
              << "gpu_callback_fallback_successes=" << std::count_if(
                     fallback_batch.begin(), fallback_batch.end(), [](const auto& value) {
                         return value.status == EvaluationStatus::Success;
                     }) << '\n'
              << "first_batch_successes=" << std::count_if(first.begin(), first.end(), [](const auto& value) {
                     return value.status == EvaluationStatus::Success;
                 }) << '\n';
    return 0;
}
