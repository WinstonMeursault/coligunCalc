#include "coilgun/optimization/genetic_operators.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace coilgun::optimization {

Population Population::initialize(const VariableSchema& schema, std::size_t size, RandomContext rng) {
    Population p;
    p.candidates_.reserve(size);
    for (std::size_t n = 0; n < size; ++n) {
        Candidate c; c.id = static_cast<CandidateId>(n);
        c.variables.values.reserve(schema.size());
        for (const auto& v : schema.variables()) {
            if (v.type == VariableType::Continuous) c.variables.values.push_back(rng.uniform(v.lower_bound, v.upper_bound));
            else c.variables.values.push_back(static_cast<double>(std::uniform_int_distribution<long long>(
                static_cast<long long>(v.lower_bound), static_cast<long long>(v.upper_bound))(rng.engine())));
        }
        c.variables = schema.repair(c.variables);
        p.push_back(std::move(c));
    }
    return p;
}

std::vector<Candidate> Population::elites(std::size_t count, const FeasibilityComparator& comparator) const {
    count = std::min(count, size());
    std::vector<Candidate> result(candidates_.begin(), candidates_.end());
    std::stable_sort(result.begin(), result.end(), [&](const Candidate& a, const Candidate& b) {
        return comparator.compare(a, b) < 0;
    });
    result.resize(count);
    return result;
}

Candidate tournament_select(const Population& population, const FeasibilityComparator& comparator,
                            RandomContext& rng, std::size_t tournament_size) {
    if (population.empty()) throw std::invalid_argument("cannot select from empty population");
    Candidate best = population[rng.index(population.size())];
    for (std::size_t i = 1; i < std::max<std::size_t>(1, tournament_size); ++i) {
        const Candidate& contender = population[rng.index(population.size())];
        if (comparator.better(contender, best)) best = contender;
    }
    return best;
}

CandidateVariables sbx_crossover(const CandidateVariables& a, const CandidateVariables& b,
                                 const VariableSchema& schema, RandomContext& rng,
                                 double rate, double di) {
    if (a.values.size() != schema.size() || b.values.size() != schema.size()) throw std::invalid_argument("parent dimensionality mismatch");
    if (rate < 0 || rate > 1 || di <= 0 || !std::isfinite(di)) throw std::invalid_argument("invalid crossover parameters");
    if (rng.uniform() > rate) return schema.repair(a);
    CandidateVariables child; child.values.resize(schema.size());
    for (std::size_t i = 0; i < schema.size(); ++i) {
        const auto& v = schema.at(i); const double x1 = a.values[i], x2 = b.values[i];
        if (v.type != VariableType::Continuous || std::abs(x1 - x2) < 1e-15) {
            child.values[i] = rng.uniform() < 0.5 ? x1 : x2; continue;
        }
        const double u = rng.uniform(); const double beta = u <= 0.5 ? std::pow(2*u, 1.0/(di+1)) : std::pow(1.0/(2*(1-u)), 1.0/(di+1));
        child.values[i] = 0.5 * ((x1 + x2) - beta * (x2 - x1));
    }
    return schema.repair(child);
}

void polynomial_mutation(CandidateVariables& values, const VariableSchema& schema, RandomContext& rng,
                         double rate, double di) {
    if (values.values.size() != schema.size()) throw std::invalid_argument("candidate dimensionality mismatch");
    if (rate < 0 || rate > 1 || di <= 0 || !std::isfinite(di)) throw std::invalid_argument("invalid mutation parameters");
    for (std::size_t i = 0; i < schema.size(); ++i) {
        const auto& v = schema.at(i); if (rng.uniform() > rate) continue;
        if (v.type == VariableType::Continuous) {
            const double y = std::clamp(values.values[i], v.lower_bound, v.upper_bound);
            const double d = v.upper_bound - v.lower_bound; if (d == 0) continue;
            const double u = rng.uniform(), delta = u < 0.5 ? std::pow(2*u, 1.0/(di+1))-1 : 1-std::pow(2*(1-u), 1.0/(di+1));
            values.values[i] = y + delta*d;
        } else {
            const long long lo = static_cast<long long>(v.lower_bound), hi = static_cast<long long>(v.upper_bound);
            values.values[i] = static_cast<double>(std::uniform_int_distribution<long long>(lo, hi)(rng.engine()));
        }
    }
    values = schema.repair(values);
}
}
