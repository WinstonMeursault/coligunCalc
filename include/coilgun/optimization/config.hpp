#pragma once

#include <cstddef>
#include <cstdint>

namespace coilgun::optimization {

enum class SelectionStrategy { Auto, SingleObjective, NSGA2 };

struct OptimizationConfig {
    static constexpr std::size_t population_size_default = 100;

    std::size_t population_size = population_size_default;
    std::size_t max_generations = 100;
    double crossover_rate = 0.9;
    double mutation_rate = 0.1;
    std::size_t elite_count = 1;
    std::uint64_t random_seed = 0;
    SelectionStrategy strategy = SelectionStrategy::Auto;

    static OptimizationConfig defaults();
    void validate() const;
};

} // namespace coilgun::optimization
