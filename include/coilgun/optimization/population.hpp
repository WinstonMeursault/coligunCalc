#pragma once
#include "coilgun/optimization/types.hpp"
#include "coilgun/optimization/variables.hpp"
#include "coilgun/optimization/comparator.hpp"
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace coilgun::optimization {

class RandomContext {
public:
    explicit RandomContext(std::uint64_t seed = 0) : engine_(seed) {}
    double uniform(double lo = 0.0, double hi = 1.0) { return std::uniform_real_distribution<double>(lo, hi)(engine_); }
    std::size_t index(std::size_t n) { return n ? std::uniform_int_distribution<std::size_t>(0, n - 1)(engine_) : 0; }
    std::mt19937_64& engine() { return engine_; }
private: std::mt19937_64 engine_;
};

class Population {
public:
    using container_type = std::vector<Candidate>;
    using iterator = container_type::iterator;
    using const_iterator = container_type::const_iterator;
    static Population initialize(const VariableSchema& schema, std::size_t size, RandomContext rng);
    std::size_t size() const noexcept { return candidates_.size(); }
    bool empty() const noexcept { return candidates_.empty(); }
    Candidate& operator[](std::size_t i) { return candidates_[i]; }
    const Candidate& operator[](std::size_t i) const { return candidates_[i]; }
    void push_back(Candidate c) { candidates_.push_back(std::move(c)); }
    iterator begin() { return candidates_.begin(); }
    iterator end() { return candidates_.end(); }
    const_iterator begin() const { return candidates_.begin(); }
    const_iterator end() const { return candidates_.end(); }
    std::vector<Candidate> elites(std::size_t count, const FeasibilityComparator& comparator) const;
private:
    container_type candidates_;
};

Candidate tournament_select(const Population& population, const FeasibilityComparator& comparator,
                            RandomContext& rng, std::size_t tournament_size = 2);

} // namespace coilgun::optimization
