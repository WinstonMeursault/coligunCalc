#include "coilgun/optimization/genetic_optimizer.hpp"

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

struct EvaluatedBatch {
    std::size_t successful = 0;
    std::size_t failed = 0;
    bool objective_schema_set = false;
    bool schema_error = false;
    std::string objective_id;
    bool maximize = true;
};

EvaluatedBatch assign_results(Population& population, const std::vector<EvaluationResult>& results,
                              std::optional<std::pair<std::string, bool>> objective_schema) {
    EvaluatedBatch summary;
    auto active_schema = objective_schema;
    for (std::size_t i = 0; i < population.size(); ++i) {
        Candidate& candidate = population[i];
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
        if (evaluation->objectives.size() != 1) {
            mark_failure(candidate, EvaluationStatus::Invalid, "single_objective_required",
                         "single-objective optimization requires exactly one objective");
            summary.schema_error = true;
            ++summary.failed;
            continue;
        }
        const auto& objective = evaluation->objectives.front();
        if (objective.id.empty() || !std::isfinite(objective.value)) {
            mark_failure(candidate, EvaluationStatus::Invalid, "invalid_objective",
                         "objective id must be non-empty and objective value finite");
            summary.schema_error = true;
            ++summary.failed;
            continue;
        }
        if (!active_schema) {
            summary.objective_schema_set = true;
            summary.objective_id = objective.id;
            summary.maximize = objective.maximize;
            active_schema = std::make_pair(objective.id, objective.maximize);
        } else if (objective.id != active_schema->first || objective.maximize != active_schema->second) {
            mark_failure(candidate, EvaluationStatus::Invalid, "objective_schema_mismatch",
                         "objective id and direction must remain fixed during a run");
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

void copy_evaluator_statistics(OptimizationResult& result, const BatchEvaluator& evaluator) {
    const EvaluationStatistics* statistics = nullptr;
    if (const auto* tracked = dynamic_cast<const StatisticsBatchEvaluator*>(&evaluator))
        statistics = &tracked->statistics();
    else if (const auto* cached = dynamic_cast<const CachedBatchEvaluator*>(&evaluator))
        statistics = &cached->statistics();
    if (!statistics) return;
    result.statistics.evaluations = statistics->evaluations;
    result.statistics.successful_evaluations = statistics->successful_evaluations;
    result.statistics.failed_evaluations = statistics->failed_evaluations;
    result.statistics.cache_hits = statistics->cache_hits;
    result.statistics.gpu_fallbacks = statistics->fallbacks;
    result.statistics.elapsed_seconds = statistics->elapsed_seconds;
}

std::vector<EvaluationResult> evaluate_safely(BatchEvaluator& evaluator,
                                               const std::vector<CandidateVariables>& variables,
                                               const EvaluationContext& context) {
    try {
        return evaluator.evaluate_batch(variables, context);
    } catch (const std::exception& error) {
        std::vector<EvaluationResult> results;
        results.reserve(variables.size());
        for (std::size_t i = 0; i < variables.size(); ++i)
            results.push_back(EvaluationResult::failed("evaluation_exception", error.what()));
        return results;
    } catch (...) {
        std::vector<EvaluationResult> results;
        results.reserve(variables.size());
        for (std::size_t i = 0; i < variables.size(); ++i)
            results.push_back(EvaluationResult::failed("evaluation_exception", "unknown exception"));
        return results;
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
    try {
        config_.validate();
        termination_.validate();
        if (config_.strategy == SelectionStrategy::NSGA2)
            throw std::invalid_argument("NSGA2 strategy is not valid for single-objective optimization");
    } catch (const std::exception& error) {
        result.termination.reason = TerminationReason::ConfigurationError;
        result.termination.message = error.what();
        return result;
    }

    RandomContext rng(config_.random_seed);
    Population population = Population::initialize(schema_, config_.population_size, rng);
    const std::size_t max_generations = termination_.max_generations == 0
        ? config_.max_generations : termination_.max_generations;
    std::optional<std::pair<std::string, bool>> objective_schema;
    std::optional<Candidate> best;
    std::size_t no_improvement = 0;
    std::size_t evaluations = 0;
    std::uint64_t next_id = static_cast<std::uint64_t>(config_.population_size);

    for (std::size_t generation = 0; generation < max_generations; ++generation) {
        const std::size_t remaining = termination_.max_evaluations == 0
            ? population.size() : (evaluations >= termination_.max_evaluations
                ? 0 : termination_.max_evaluations - evaluations);
        std::vector<EvaluationResult> evaluated;
        if (remaining >= population.size()) {
            evaluated = evaluate_safely(*evaluator_, variables_for(population),
                                         EvaluationContext{config_.random_seed, false});
            evaluations += population.size();
        } else if (remaining > 0) {
            std::vector<CandidateVariables> prefix;
            prefix.reserve(remaining);
            for (std::size_t i = 0; i < remaining; ++i) prefix.push_back(population[i].variables);
            evaluated = evaluate_safely(*evaluator_, prefix, EvaluationContext{config_.random_seed, false});
            evaluations += remaining;
            evaluated.resize(population.size());
            for (std::size_t i = remaining; i < population.size(); ++i)
                evaluated[i] = EvaluationResult::failed("evaluation_budget", "evaluation budget exhausted");
        } else {
            evaluated.resize(population.size(), EvaluationResult::failed(
                "evaluation_budget", "evaluation budget exhausted"));
        }

        const auto summary = assign_results(population, evaluated, objective_schema);
        if (!objective_schema && summary.objective_schema_set)
            objective_schema = std::make_pair(summary.objective_id, summary.maximize);
        result.statistics.evaluations = evaluations;
        result.statistics.successful_evaluations += summary.successful;
        result.statistics.failed_evaluations += summary.failed;
        result.statistics.generations = generation + 1;
        copy_evaluator_statistics(result, *evaluator_);

        if (summary.schema_error) {
            result.termination = {TerminationReason::ConfigurationError,
                                  "single-objective objective schema is invalid", generation};
            return result;
        }
        if (summary.successful == 0) {
            result.termination = {TerminationReason::EvaluationFailure,
                                  "all candidates failed evaluation", generation};
            return result;
        }

        const Candidate current_best = best_candidate(population, comparator_);
        if (current_best.evaluation_status == EvaluationStatus::Success) {
            const bool improved = !best || objective_is_better(current_best, *best, termination_.improvement_tolerance) ||
                                  comparator_.better(current_best, *best);
            if (improved) {
                best = current_best;
                no_improvement = 0;
            } else {
                ++no_improvement;
            }
            fill_result(result, *best);

            if (termination_.target_value) {
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
        const auto elites = population.elites(config_.elite_count, comparator_);
        for (auto elite : elites) next.push_back(std::move(elite));
        while (next.size() < config_.population_size) {
            const auto parent_a = tournament_select(population, comparator_, rng);
            const auto parent_b = tournament_select(population, comparator_, rng);
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
    return GeneticOptimizer(schema, evaluator, config, termination, comparator).optimize();
}

} // namespace coilgun::optimization
