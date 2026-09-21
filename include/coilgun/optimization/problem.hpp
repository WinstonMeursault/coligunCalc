#pragma once

#include "coilgun/optimization/constraint.hpp"
#include "coilgun/optimization/objective.hpp"
#include "coilgun/optimization/types.hpp"
#include "coilgun/optimization/variables.hpp"

#include <functional>
#include <optional>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace coilgun::optimization {

/** A callable candidate-repair policy used by a ProblemSpec. */
class RepairPolicy {
public:
    using Function = std::function<CandidateVariables(const VariableSchema&, const CandidateVariables&)>;

    RepairPolicy() = default;
    RepairPolicy(Function function) : function_(std::move(function)) {}
    RepairPolicy(std::function<CandidateVariables(const CandidateVariables&)> function)
        : function_([function = std::move(function)](const VariableSchema&, const CandidateVariables& candidate) {
              return function(candidate);
          }) {}

    template <typename Callable,
              std::enable_if_t<!std::is_same_v<std::decay_t<Callable>, RepairPolicy> &&
                                   (std::is_invocable_r_v<CandidateVariables, Callable,
                                                          const VariableSchema&, const CandidateVariables&> ||
                                    std::is_invocable_r_v<CandidateVariables, Callable,
                                                          const CandidateVariables&>), int> = 0>
    RepairPolicy(Callable&& callable) {
        if constexpr (std::is_invocable_r_v<CandidateVariables, Callable,
                                            const VariableSchema&, const CandidateVariables&>) {
            function_ = std::forward<Callable>(callable);
        } else {
            function_ = [callable = std::forward<Callable>(callable)](
                            const VariableSchema&, const CandidateVariables& candidate) mutable {
                return callable(candidate);
            };
        }
    }

    [[nodiscard]] explicit operator bool() const noexcept { return static_cast<bool>(function_); }
    [[nodiscard]] CandidateVariables apply(const VariableSchema& schema,
                                           const CandidateVariables& candidate) const;
    [[nodiscard]] CandidateVariables operator()(const VariableSchema& schema,
                                                const CandidateVariables& candidate) const {
        return apply(schema, candidate);
    }

private:
    Function function_;
};

/** Immutable run contract for variables, objectives, constraints, and repair. */
class ProblemSpec {
public:
    ProblemSpec(VariableSchema schema,
                std::vector<ObjectiveDefinition> objectives,
                std::vector<ConstraintDefinition> constraints = {});
    template <typename Callable,
              std::enable_if_t<!std::is_same_v<std::decay_t<Callable>, ProblemSpec>, int> = 0>
    ProblemSpec(VariableSchema schema,
                std::vector<ObjectiveDefinition> objectives,
                Callable&& repair_policy)
        : ProblemSpec(std::move(schema), std::move(objectives), {},
                      RepairPolicy(std::forward<Callable>(repair_policy))) {}
    ProblemSpec(VariableSchema schema,
                std::vector<ObjectiveDefinition> objectives,
                std::vector<ConstraintDefinition> constraints,
                RepairPolicy repair_policy);

    [[nodiscard]] const VariableSchema& schema() const noexcept { return schema_; }
    [[nodiscard]] const std::vector<ObjectiveDefinition>& objectives() const noexcept { return objectives_; }
    [[nodiscard]] const std::vector<ConstraintDefinition>& constraints() const noexcept { return constraints_; }
    [[nodiscard]] const RepairPolicy& repair_policy() const noexcept { return repair_policy_; }
    [[nodiscard]] CandidateVariables repair(const CandidateVariables& candidate) const;

private:
    VariableSchema schema_;
    std::vector<ObjectiveDefinition> objectives_;
    std::vector<ConstraintDefinition> constraints_;
    RepairPolicy repair_policy_;
};

class OptimizationProblem {
public:
    OptimizationProblem() = default;
    explicit OptimizationProblem(ProblemSpec spec) : spec_(std::move(spec)) {}
    virtual ~OptimizationProblem() = default;
    virtual EvaluationResult evaluate(const CandidateVariables& variables) const = 0;

    /** Null for legacy schema-first problems; non-null for spec-aware problems. */
    [[nodiscard]] virtual const ProblemSpec* spec() const noexcept {
        return spec_ ? &*spec_ : nullptr;
    }

protected:
    void set_spec(ProblemSpec spec) { spec_ = std::move(spec); }

private:
    std::optional<ProblemSpec> spec_;
};

} // namespace coilgun::optimization
