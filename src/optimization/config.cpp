#include "coilgun/optimization/config.hpp"

#include <cmath>
#include <stdexcept>

namespace coilgun::optimization {

OptimizationConfig OptimizationConfig::defaults() { return {}; }

void OptimizationConfig::validate() const {
    if (population_size == 0) throw std::invalid_argument("population_size must be positive");
    if (max_generations == 0) throw std::invalid_argument("max_generations must be positive");
    if (!std::isfinite(crossover_rate) || crossover_rate < 0.0 || crossover_rate > 1.0)
        throw std::invalid_argument("crossover_rate must be finite and in [0, 1]");
    if (!std::isfinite(mutation_rate) || mutation_rate < 0.0 || mutation_rate > 1.0)
        throw std::invalid_argument("mutation_rate must be finite and in [0, 1]");
    if (elite_count > population_size) throw std::invalid_argument("elite_count cannot exceed population_size");
}

} // namespace coilgun::optimization
