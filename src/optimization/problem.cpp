#include "coilgun/optimization/problem.hpp"

#include <stdexcept>

namespace coilgun::optimization {

CandidateVariables RepairPolicy::apply(const VariableSchema& schema,
                                       const CandidateVariables& candidate) const {
    if (!function_) throw std::invalid_argument("repair policy must be callable");
    return function_(schema, candidate);
}

ProblemSpec::ProblemSpec(VariableSchema schema,
                         std::vector<ObjectiveDefinition> objectives,
                         std::vector<ConstraintDefinition> constraints)
    : ProblemSpec(std::move(schema), std::move(objectives), std::move(constraints),
                  RepairPolicy([](const VariableSchema& owner,
                                  const CandidateVariables& candidate) {
                      return owner.repair(candidate);
                  })) {}

ProblemSpec::ProblemSpec(VariableSchema schema,
                         std::vector<ObjectiveDefinition> objectives,
                         std::vector<ConstraintDefinition> constraints,
                         RepairPolicy repair_policy)
    : schema_(std::move(schema)), objectives_(std::move(objectives)),
      constraints_(std::move(constraints)), repair_policy_(std::move(repair_policy)) {
    if (objectives_.empty()) throw std::invalid_argument("ProblemSpec requires at least one objective");
    std::unordered_set<std::string> objective_ids;
    for (const auto& objective : objectives_) {
        objective.validate();
        if (!objective_ids.insert(objective.id).second)
            throw std::invalid_argument("duplicate objective ID: " + objective.id);
    }

    std::unordered_set<std::string> constraint_ids;
    for (const auto& constraint : constraints_) {
        constraint.validate();
        if (!constraint_ids.insert(constraint.id).second)
            throw std::invalid_argument("duplicate constraint ID: " + constraint.id);
    }

    if (!repair_policy_) throw std::invalid_argument("repair policy must be callable");
}

CandidateVariables ProblemSpec::repair(const CandidateVariables& candidate) const {
    // The policy may apply coupled/domain-specific changes, but the schema owns
    // the final dimensionality and variable-domain invariants.
    return schema_.repair(repair_policy_.apply(schema_, candidate));
}

} // namespace coilgun::optimization
