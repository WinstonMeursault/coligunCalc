#pragma once

#include "coilgun/optimization/config.hpp"
#include "coilgun/optimization/evaluator.hpp"
#include "coilgun/optimization/genetic_operators.hpp"
#include "coilgun/optimization/problem.hpp"
#include "coilgun/optimization/termination.hpp"

#include <memory>
#include <type_traits>
#include <utility>

namespace coilgun::optimization {

class GeneticOptimizer {
public:
    GeneticOptimizer(VariableSchema schema, std::shared_ptr<BatchEvaluator> evaluator,
                     OptimizationConfig config = OptimizationConfig::defaults(),
                     TerminationConfig termination = {},
                     FeasibilityComparator comparator = FeasibilityComparator{});
    GeneticOptimizer(VariableSchema schema, BatchEvaluator& evaluator,
                     OptimizationConfig config = OptimizationConfig::defaults(),
                     TerminationConfig termination = {},
                     FeasibilityComparator comparator = FeasibilityComparator{});
    GeneticOptimizer(VariableSchema schema, const OptimizationProblem& problem,
                     OptimizationConfig config = OptimizationConfig::defaults(),
                     TerminationConfig termination = {},
                     FeasibilityComparator comparator = FeasibilityComparator{});
    template <typename Problem,
              std::enable_if_t<std::is_base_of_v<OptimizationProblem, Problem> &&
                                   std::is_base_of_v<BatchEvaluator, Problem>,
                               int> = 0>
    GeneticOptimizer(VariableSchema schema, Problem& problem,
                     OptimizationConfig config = OptimizationConfig::defaults(),
                     TerminationConfig termination = {},
                     FeasibilityComparator comparator = FeasibilityComparator{})
        : GeneticOptimizer(std::move(schema), static_cast<BatchEvaluator&>(problem), config,
                           std::move(termination), std::move(comparator)) {}

    [[nodiscard]] OptimizationResult optimize();
    [[nodiscard]] OptimizationResult run() { return optimize(); }

private:
    VariableSchema schema_;
    std::shared_ptr<BatchEvaluator> owned_evaluator_;
    BatchEvaluator* evaluator_ = nullptr;
    OptimizationConfig config_;
    TerminationConfig termination_;
    FeasibilityComparator comparator_;
};

using SingleObjectiveOptimizer = GeneticOptimizer;

OptimizationResult optimize_single_objective(
    const VariableSchema& schema, BatchEvaluator& evaluator,
    const OptimizationConfig& config = OptimizationConfig::defaults(),
    const TerminationConfig& termination = {},
    const FeasibilityComparator& comparator = FeasibilityComparator{});

} // namespace coilgun::optimization
