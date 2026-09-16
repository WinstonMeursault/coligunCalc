#include "coilgun/optimization/genetic_optimizer.hpp"

#include "coilgun/optimization/nsga2.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace coilgun::optimization {
namespace {

class ProblemBatchEvaluator final : public BatchEvaluator {
public:
    explicit ProblemBatchEvaluator(const OptimizationProblem& problem) : problem_(problem) {}

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& variables,
                                                 const EvaluationContext&) override {
        std::vector<EvaluationResult> results;
        results.reserve(variables.size());
        for (const auto& value : variables) {
            try {
                results.push_back(problem_.evaluate(value));
            } catch (const std::exception& error) {
                results.push_back(EvaluationResult::failed("evaluation_exception", error.what()));
            } catch (...) {
                results.push_back(EvaluationResult::failed("evaluation_exception", "unknown exception"));
            }
        }
        return results;
    }

private:
    const OptimizationProblem& problem_;
};

void mark_failure(Candidate& candidate, EvaluationStatus status, std::string code, std::string message) {
    candidate.evaluation_status = status;
    candidate.objectives.clear();
    candidate.constraints.clear();
    candidate.diagnostics.clear();
    candidate.diagnostics.push_back({std::move(code), std::move(message), DiagnosticSeverity::Error});
}

bool objective_is_better(const Candidate& lhs, const Candidate& rhs, double tolerance) {
    if (lhs.evaluation_status != EvaluationStatus::Success || rhs.evaluation_status != EvaluationStatus::Success ||
        lhs.objectives.empty() || rhs.objectives.empty()) return false;
    const auto& l = lhs.objectives.front();
    const auto& r = rhs.objectives.front();
    if (l.maximize != r.maximize) return false;
    return l.maximize ? l.value > r.value + tolerance : l.value < r.value - tolerance;
}

bool non_objective_is_better(const Candidate& lhs, const Candidate& rhs,
                             const FeasibilityComparator& comparator) {
    Candidate lhs_without_objective = lhs;
    Candidate rhs_without_objective = rhs;
    lhs_without_objective.objectives.clear();
    rhs_without_objective.objectives.clear();
    return comparator.better(lhs_without_objective, rhs_without_objective);
}

bool non_objective_is_tied(const Candidate& lhs, const Candidate& rhs,
                           const FeasibilityComparator& comparator) {
    return !non_objective_is_better(lhs, rhs, comparator) &&
           !non_objective_is_better(rhs, lhs, comparator);
}

struct EvaluatedBatch {
    std::size_t successful = 0;
    std::size_t failed = 0;
    bool objective_schema_set = false;
    bool schema_error = false;
    std::vector<std::pair<std::string, bool>> objective_schema;
};

EvaluatedBatch assign_results(Population& population, const std::vector<EvaluationResult>& results,
                              std::optional<std::vector<std::pair<std::string, bool>>> objective_schema,
                              std::size_t evaluated_count) {
    EvaluatedBatch summary;
    auto active_schema = objective_schema;
    for (std::size_t i = 0; i < population.size(); ++i) {
        Candidate& candidate = population[i];
        if (i >= evaluated_count) {
            candidate.evaluation_status = EvaluationStatus::Unevaluated;
            candidate.objectives.clear();
            candidate.constraints.clear();
            candidate.diagnostics.clear();
            candidate.metadata.clear();
            continue;
        }
        const EvaluationResult* evaluation = i < results.size() ? &results[i] : nullptr;
        if (!evaluation) {
            mark_failure(candidate, EvaluationStatus::Failed, "evaluation_batch_output",
                         "batch evaluator returned too few results");
            ++summary.failed;
            continue;
        }
        candidate.evaluation_status = evaluation->status;
        candidate.objectives = evaluation->objectives;
        candidate.constraints = evaluation->constraints;
        candidate.diagnostics = evaluation->diagnostics;
        candidate.metadata = evaluation->metadata;
        if (evaluation->status != EvaluationStatus::Success) {
            ++summary.failed;
            continue;
        }
        if (evaluation->objectives.empty()) {
            mark_failure(candidate, EvaluationStatus::Invalid, "single_objective_required",
                         "optimization requires at least one objective");
            summary.schema_error = true;
            ++summary.failed;
            continue;
        }
        if (!active_schema) {
            std::vector<std::pair<std::string, bool>> discovered;
            discovered.reserve(evaluation->objectives.size());
            for (const auto& objective : evaluation->objectives)
                discovered.emplace_back(objective.id, objective.maximize);
            active_schema = discovered;
            summary.objective_schema_set = true;
            summary.objective_schema = discovered;
        }
        bool valid_objectives = evaluation->objectives.size() == active_schema->size();
        if (valid_objectives) {
            for (std::size_t objective_index = 0; objective_index < evaluation->objectives.size(); ++objective_index) {
                const auto& objective = evaluation->objectives[objective_index];
                if (objective.id.empty() || !std::isfinite(objective.value) ||
                    objective.id != (*active_schema)[objective_index].first ||
                    objective.maximize != (*active_schema)[objective_index].second) {
                    valid_objectives = false;
                    break;
                }
            }
        }
        if (!valid_objectives) {
            mark_failure(candidate, EvaluationStatus::Invalid, "objective_schema_mismatch",
                         "objective count, ids, and directions must remain fixed during a run");
            summary.schema_error = true;
            ++summary.failed;
            continue;
        }
        ++summary.successful;
    }
    return summary;
}

Candidate best_candidate(const Population& population, const FeasibilityComparator& comparator) {
    if (population.empty()) return {};
    return *std::min_element(population.begin(), population.end(),
        [&](const Candidate& lhs, const Candidate& rhs) { return comparator.better(lhs, rhs); });
}

std::vector<CandidateVariables> variables_for(const Population& population) {
    std::vector<CandidateVariables> values;
    values.reserve(population.size());
    for (const auto& candidate : population) values.push_back(candidate.variables);
    return values;
}

void fill_result(OptimizationResult& result, const Candidate& best) {
    if (best.evaluation_status != EvaluationStatus::Success || best.objectives.size() != 1) return;
    result.pareto_front = {best};
    result.best_by_objective[best.objectives.front().id] = best;
}

std::vector<ObjectiveDefinition> objective_definitions(
    const std::vector<std::pair<std::string, bool>>& schema) {
    std::vector<ObjectiveDefinition> definitions;
    definitions.reserve(schema.size());
    for (const auto& [id, maximize] : schema) definitions.push_back({id, maximize, 1.0});
    return definitions;
}

void fill_multi_result(OptimizationResult& result, const Population& population,
                       const std::vector<ObjectiveDefinition>& definitions,
                       const FeasibilityComparator& comparator) {
    std::vector<Candidate> candidates(population.begin(), population.end());
    if (candidates.empty()) return;
    const auto ranking = nsga2_rank(candidates, definitions, comparator);
    if (ranking.fronts.empty()) return;
    result.pareto_front.clear();
    for (const auto index : ranking.fronts.front()) result.pareto_front.push_back(candidates[index]);
}

EvaluationStatistics subtract_statistics(const EvaluationStatistics& after, const EvaluationStatistics& before) {
    EvaluationStatistics delta;
    delta.seed = after.seed;
    delta.evaluations = after.evaluations - before.evaluations;
    delta.successful_evaluations = after.successful_evaluations - before.successful_evaluations;
    delta.failed_evaluations = after.failed_evaluations - before.failed_evaluations;
    delta.cache_hits = after.cache_hits - before.cache_hits;
    delta.fallbacks = after.fallbacks - before.fallbacks;
    delta.elapsed_seconds = after.elapsed_seconds - before.elapsed_seconds;
    return delta;
}

std::vector<EvaluationResult> retry_singletons(BatchEvaluator& evaluator,
                                                const std::vector<CandidateVariables>& variables,
                                                const EvaluationContext& context) {
    std::vector<EvaluationResult> results;
    results.reserve(variables.size());
    for (const auto& variable : variables) {
        try {
            auto single = evaluator.evaluate_batch({variable}, context);
            if (single.size() == 1) results.push_back(std::move(single.front()));
            else results.push_back(EvaluationResult::failed("evaluation_batch_output", "singleton retry returned wrong result count"));
        } catch (const std::exception& single_error) {
            results.push_back(EvaluationResult::failed("evaluation_exception", single_error.what()));
        } catch (...) {
            results.push_back(EvaluationResult::failed("evaluation_exception", "unknown exception"));
        }
    }
    return results;
}

std::vector<EvaluationResult> evaluate_safely(BatchEvaluator& evaluator,
                                               const std::vector<CandidateVariables>& variables,
                                               const EvaluationContext& context) {
    try {
        return evaluator.evaluate_batch(variables, context);
    } catch (const std::exception&) {
        return retry_singletons(evaluator, variables, context);
    } catch (...) {
        return retry_singletons(evaluator, variables, context);
    }
}

} // namespace

GeneticOptimizer::GeneticOptimizer(VariableSchema schema, std::shared_ptr<BatchEvaluator> evaluator,
                                   OptimizationConfig config, TerminationConfig termination,
                                   FeasibilityComparator comparator)
    : schema_(std::move(schema)), owned_evaluator_(std::move(evaluator)), evaluator_(owned_evaluator_.get()),
      config_(config), termination_(std::move(termination)), comparator_(std::move(comparator)) {
    if (!evaluator_) throw std::invalid_argument("evaluator must not be null");
}

GeneticOptimizer::GeneticOptimizer(VariableSchema schema, BatchEvaluator& evaluator,
                                   OptimizationConfig config, TerminationConfig termination,
                                   FeasibilityComparator comparator)
    : schema_(std::move(schema)), evaluator_(&evaluator), config_(config),
      termination_(std::move(termination)), comparator_(std::move(comparator)) {}

GeneticOptimizer::GeneticOptimizer(VariableSchema schema, const OptimizationProblem& problem,
                                   OptimizationConfig config, TerminationConfig termination,
                                   FeasibilityComparator comparator)
    : schema_(std::move(schema)), owned_evaluator_(std::make_shared<ProblemBatchEvaluator>(problem)),
      evaluator_(owned_evaluator_.get()), config_(config), termination_(std::move(termination)),
      comparator_(std::move(comparator)) {}

OptimizationResult GeneticOptimizer::optimize() {
    OptimizationResult result;
    result.statistics.seed = config_.random_seed;
    try {
        config_.validate();
        termination_.validate();
    } catch (const std::exception& error) {
        result.termination.reason = TerminationReason::ConfigurationError;
        result.termination.message = error.what();
        return result;
    }

    RandomContext rng(config_.random_seed);
    Population population = Population::initialize(schema_, config_.population_size, rng);
    const std::size_t max_generations = termination_.max_generations == 0
        ? config_.max_generations : termination_.max_generations;
    std::optional<std::vector<std::pair<std::string, bool>>> objective_schema;
    std::optional<SelectionStrategy> resolved_strategy;
    std::optional<Population> pending_parents;
    std::optional<Candidate> best;
    std::size_t no_improvement = 0;
    std::size_t evaluations = 0;
    std::uint64_t next_id = static_cast<std::uint64_t>(config_.population_size);
    const auto evaluator_before = evaluator_->statistics_snapshot();

    for (std::size_t generation = 0; generation < max_generations; ++generation) {
        const std::size_t remaining = termination_.max_evaluations == 0
            ? population.size() : (evaluations >= termination_.max_evaluations
                ? 0 : termination_.max_evaluations - evaluations);
        std::vector<EvaluationResult> evaluated;
        std::size_t evaluated_count = 0;
        if (remaining >= population.size()) {
            evaluated = evaluate_safely(*evaluator_, variables_for(population),
                                         EvaluationContext{config_.random_seed, false});
            evaluations += population.size();
            evaluated_count = population.size();
        } else if (remaining > 0) {
            std::vector<CandidateVariables> prefix;
            prefix.reserve(remaining);
            for (std::size_t i = 0; i < remaining; ++i) prefix.push_back(population[i].variables);
            evaluated = evaluate_safely(*evaluator_, prefix, EvaluationContext{config_.random_seed, false});
            evaluations += remaining;
            evaluated_count = remaining;
            evaluated.resize(population.size());
            for (std::size_t i = remaining; i < population.size(); ++i)
                evaluated[i] = EvaluationResult::failed("evaluation_budget", "evaluation budget exhausted");
        } else {
            evaluated.resize(population.size(), EvaluationResult::failed(
                "evaluation_budget", "evaluation budget exhausted"));
        }

        const auto summary = assign_results(population, evaluated, objective_schema, evaluated_count);
        if (!objective_schema && summary.objective_schema_set)
            objective_schema = summary.objective_schema;
        result.statistics.evaluations = evaluations;
        result.statistics.successful_evaluations += summary.successful;
        result.statistics.failed_evaluations += summary.failed;
        result.statistics.skipped_due_to_budget += population.size() - evaluated_count;
        result.statistics.generations = generation + 1;
        if (const auto evaluator_after = evaluator_->statistics_snapshot()) {
            const auto baseline = evaluator_before.value_or(EvaluationStatistics{});
            const auto delta = subtract_statistics(*evaluator_after, baseline);
            result.statistics.cache_hits = delta.cache_hits;
            result.statistics.gpu_fallbacks = delta.fallbacks;
            result.statistics.elapsed_seconds = delta.elapsed_seconds;
        }

        if (summary.schema_error) {
            result.termination = {TerminationReason::ConfigurationError,
                                  "objective count, ids, and directions must remain fixed during a run", generation};
            return result;
        }
        if (summary.successful == 0) {
            result.termination = {TerminationReason::EvaluationFailure,
                                  "all candidates failed evaluation", generation};
            return result;
        }

        if (!resolved_strategy && objective_schema) {
            const std::size_t objective_count = objective_schema->size();
            if (config_.strategy == SelectionStrategy::Auto) {
                resolved_strategy = objective_count == 1 ? SelectionStrategy::SingleObjective
                                                         : SelectionStrategy::NSGA2;
            } else if (config_.strategy == SelectionStrategy::SingleObjective && objective_count != 1) {
                result.termination = {TerminationReason::ConfigurationError,
                                      "SingleObjective strategy requires exactly one objective", generation};
                return result;
            } else if (config_.strategy == SelectionStrategy::NSGA2 && objective_count < 2) {
                result.termination = {TerminationReason::ConfigurationError,
                                      "NSGA2 strategy requires at least two objectives", generation};
                return result;
            } else {
                resolved_strategy = config_.strategy;
            }
        }

        if (resolved_strategy == SelectionStrategy::NSGA2 && pending_parents) {
            population = nsga2_select(*pending_parents, population, config_.population_size,
                                      objective_definitions(*objective_schema), comparator_);
            pending_parents.reset();
        }

        if (resolved_strategy == SelectionStrategy::NSGA2) {
            fill_multi_result(result, population, objective_definitions(*objective_schema), comparator_);
        }

        const Candidate current_best = best_candidate(population, comparator_);
        if (resolved_strategy == SelectionStrategy::SingleObjective &&
            current_best.evaluation_status == EvaluationStatus::Success) {
            const bool improved = !best ||
                                  non_objective_is_better(current_best, *best, comparator_) ||
                                  (non_objective_is_tied(current_best, *best, comparator_) &&
                                   objective_is_better(current_best, *best, termination_.improvement_tolerance));
            if (improved) {
                best = current_best;
                no_improvement = 0;
            } else {
                ++no_improvement;
            }
            fill_result(result, *best);

            if (termination_.target_value && is_feasible(current_best.constraints)) {
                const double value = current_best.objectives.front().value;
                const bool reached = current_best.objectives.front().maximize
                    ? value >= *termination_.target_value : value <= *termination_.target_value;
                if (reached) {
                    result.termination = {TerminationReason::TargetReached, "target reached", generation};
                    return result;
                }
            }
        }

        if (evaluations >= termination_.max_evaluations && termination_.max_evaluations != 0) {
            // TerminationReason predates the evaluation-budget criterion; the
            // message carries the more specific reason without changing T1 API.
            result.termination = {TerminationReason::MaxGenerations, "evaluation budget exhausted", generation};
            return result;
        }
        if (termination_.max_no_improvement_generations != 0 &&
            no_improvement >= termination_.max_no_improvement_generations) {
            result.termination = {TerminationReason::Converged, "no objective improvement", generation};
            return result;
        }
        if (generation + 1 >= max_generations) {
            result.termination = {TerminationReason::MaxGenerations, "maximum generations reached", generation};
            return result;
        }

        Population next;
        if (resolved_strategy == SelectionStrategy::NSGA2) {
            pending_parents = population;
        } else {
            const auto elites = population.elites(config_.elite_count, comparator_);
            for (auto elite : elites) next.push_back(std::move(elite));
        }
        std::optional<Nsga2Ranking> mating_ranking;
        if (resolved_strategy == SelectionStrategy::NSGA2) {
            const std::vector<Candidate> candidates(population.begin(), population.end());
            mating_ranking = nsga2_rank(candidates, objective_definitions(*objective_schema), comparator_);
        }
        while (next.size() < config_.population_size) {
            const auto parent_a = resolved_strategy == SelectionStrategy::NSGA2
                ? nsga2_tournament_select(population, *mating_ranking, rng)
                : tournament_select(population, comparator_, rng);
            const auto parent_b = resolved_strategy == SelectionStrategy::NSGA2
                ? nsga2_tournament_select(population, *mating_ranking, rng)
                : tournament_select(population, comparator_, rng);
            Candidate child;
            child.id = next_id++;
            child.variables = sbx_crossover(parent_a.variables, parent_b.variables, schema_, rng, config_.crossover_rate);
            polynomial_mutation(child.variables, schema_, rng, config_.mutation_rate);
            next.push_back(std::move(child));
        }
        population = std::move(next);
    }
    result.termination = {TerminationReason::MaxGenerations, "maximum generations reached", max_generations};
    return result;
}

OptimizationResult optimize_single_objective(const VariableSchema& schema, BatchEvaluator& evaluator,
                                              const OptimizationConfig& config,
                                              const TerminationConfig& termination,
                                              const FeasibilityComparator& comparator) {
    auto single_config = config;
    single_config.strategy = SelectionStrategy::SingleObjective;
    return GeneticOptimizer(schema, evaluator, single_config, termination, comparator).optimize();
}

} // namespace coilgun::optimization
