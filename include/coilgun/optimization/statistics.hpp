#pragma once
#include <cstdint>
namespace coilgun::optimization {
struct EvaluationStatistics {
    // evaluations counts uncached delegate evaluations; cache_hits counts results served by the cache.
    std::uint64_t evaluations = 0, successful_evaluations = 0, failed_evaluations = 0;
    std::uint64_t cache_hits = 0, fallbacks = 0, seed = 0;
    double elapsed_seconds = 0.0;
};
}
