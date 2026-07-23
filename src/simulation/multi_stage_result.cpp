/**
 * @file multi_stage_result.cpp
 * @brief Multi-stage simulation result data structure implementation.
 * @author Winston Meursault
 */

#include "coilgun/simulation/multi_stage_result.hpp"

#include <stdexcept>

namespace coilgun::simulation {

MultiStageResult MultiStageResult::sampled(int every_n) const {
    if (every_n <= 0) throw std::invalid_argument("every_n must be positive");
    MultiStageResult r;
    r.summary = summary;
    r.history.reserve(history.size() / static_cast<std::size_t>(every_n) + 1);
    for (std::size_t i = 0; i < history.size(); i += static_cast<std::size_t>(every_n))
        r.history.push_back(history[i]);
    return r;
}

} // namespace coilgun::simulation
