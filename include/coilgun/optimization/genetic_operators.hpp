#pragma once
#include "coilgun/optimization/population.hpp"

namespace coilgun::optimization {
CandidateVariables sbx_crossover(const CandidateVariables& parent_a, const CandidateVariables& parent_b,
                                 const VariableSchema& schema, RandomContext& rng,
                                 double crossover_rate = 0.9, double distribution_index = 2.0);
void polynomial_mutation(CandidateVariables& values, const VariableSchema& schema, RandomContext& rng,
                         double mutation_rate = 0.1, double distribution_index = 20.0);
}
