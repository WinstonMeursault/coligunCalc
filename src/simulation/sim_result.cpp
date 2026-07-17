/**
 * @file sim_result.cpp
 * @brief Simulation result data structure implementation.
 * @author Winston Meursault
 */

#include "coilgun/simulation/sim_result.hpp"

namespace coilgun::simulation {

SimResult SimResult::sampled(int every_n) const {
    SimResult r;
    r.summary = summary;
    r.history.reserve(history.size() / static_cast<std::size_t>(every_n) + 1);
    for (std::size_t i = 0; i < history.size(); i += static_cast<std::size_t>(every_n))
        r.history.push_back(history[i]);
    return r;
}

} // namespace coilgun::simulation
