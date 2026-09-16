#pragma once

#include "coilgun/optimization/objective.hpp"
#include "coilgun/optimization/population.hpp"

#include <cstddef>
#include <vector>

namespace coilgun::optimization {

struct Nsga2Ranking {
    // Fronts and candidate-indexed arrays use the original input order.
    std::vector<std::vector<std::size_t>> fronts;
    std::vector<std::size_t> ranks;
    std::vector<double> crowding_distances;
};

// Rank candidates by constraint-domination and then Pareto-domination.
// Objective definitions are optional; when omitted, each ObjectiveValue's
// maximize flag is used and objective scales are one.
Nsga2Ranking nsga2_rank(const std::vector<Candidate>& candidates,
                        const std::vector<ObjectiveDefinition>& definitions = {},
                        const FeasibilityComparator& comparator = FeasibilityComparator{});

std::vector<std::vector<std::size_t>> non_dominated_sort(
    const std::vector<Candidate>& candidates,
    const std::vector<ObjectiveDefinition>& definitions = {},
    const FeasibilityComparator& comparator = FeasibilityComparator{});

std::vector<double> crowding_distances(
    const std::vector<Candidate>& candidates,
    const std::vector<std::size_t>& front,
    const std::vector<ObjectiveDefinition>& definitions = {});

// Select one mating parent using NSGA-II rank first and crowding distance as
// the tie-breaker. Sampling is performed uniformly from the population for
// each tournament contender, preserving the supplied random stream.
Candidate nsga2_tournament_select(
    const Population& population,
    const Nsga2Ranking& ranking,
    RandomContext& rng,
    std::size_t tournament_size = 2);

// Merge parents and offspring, rank the merged population, and retain at most
// target_size candidates. Ties are resolved in stable input order.
std::vector<Candidate> nsga2_select(
    const std::vector<Candidate>& parents,
    const std::vector<Candidate>& offspring,
    std::size_t target_size,
    const std::vector<ObjectiveDefinition>& definitions = {},
    const FeasibilityComparator& comparator = FeasibilityComparator{});

Population nsga2_select(
    const Population& parents,
    const Population& offspring,
    std::size_t target_size,
    const std::vector<ObjectiveDefinition>& definitions = {},
    const FeasibilityComparator& comparator = FeasibilityComparator{});

inline std::vector<Candidate> select_next_generation(
    const std::vector<Candidate>& parents,
    const std::vector<Candidate>& offspring,
    std::size_t target_size,
    const std::vector<ObjectiveDefinition>& definitions = {},
    const FeasibilityComparator& comparator = FeasibilityComparator{}) {
    return nsga2_select(parents, offspring, target_size, definitions, comparator);
}

} // namespace coilgun::optimization
