#include <doctest/doctest.h>
#include "coilgun/optimization/genetic_operators.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace coilgun::optimization;

TEST_CASE("population initialization is reproducible and repaired") {
    VariableSchema schema({VariableSpec::continuous("x", 0, 1), VariableSpec::integer("n", 1, 3),
                           VariableSpec::enumeration("mode", {"a", "b", "c"})});
    auto a = Population::initialize(schema, 12, RandomContext(42));
    auto b = Population::initialize(schema, 12, RandomContext(42));
    REQUIRE(a.size() == 12);
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        CHECK(a[i].variables.values == b[i].variables.values);
        CHECK(a[i].variables.values[0] >= 0);
        CHECK(a[i].variables.values[0] <= 1);
        CHECK(a[i].variables.values[1] == doctest::Approx(std::round(a[i].variables.values[1])));
    }
}

TEST_CASE("mixed variable crossover and mutation preserve schema types") {
    VariableSchema schema({VariableSpec::continuous("x", 0, 1), VariableSpec::integer("n", 1, 3),
                           VariableSpec::enumeration("mode", {"a", "b", "c"})});
    Candidate a, b; a.variables = CandidateVariables{{0.2, 1, 0}}; b.variables = CandidateVariables{{0.8, 3, 2}};
    RandomContext rng(7);
    auto child = sbx_crossover(a.variables, b.variables, schema, rng, 1.0, 2.0);
    polynomial_mutation(child, schema, rng, 1.0, 20.0);
    CHECK(child.values[0] >= 0); CHECK(child.values[0] <= 1);
    CHECK(child.values[1] >= 1); CHECK(child.values[1] <= 3);
    CHECK(child.values[1] == doctest::Approx(std::round(child.values[1])));
    CHECK(child.values[2] >= 0); CHECK(child.values[2] <= 2);
    CHECK(child.values[2] == doctest::Approx(std::round(child.values[2])));
}

TEST_CASE("zero probabilities leave genetic values unchanged") {
    VariableSchema schema({VariableSpec::continuous("x", 0, 1), VariableSpec::integer("n", 1, 3)});
    CandidateVariables x{{0.25, 2}}, y{{0.75, 3}};
    RandomContext rng(1);
    CHECK(sbx_crossover(x, y, schema, rng, 0.0, 2.0).values == x.values);
    auto before = x; polynomial_mutation(x, schema, rng, 0.0, 20.0);
    CHECK(x.values == before.values);
}

TEST_CASE("genetic operators validate finite rates and honor exact boundaries") {
    VariableSchema schema({VariableSpec::continuous("x", 0, 1)});
    CandidateVariables a{{0.25}}, b{{0.75}};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    RandomContext crossover_nan_rng(1);
    RandomContext mutation_nan_rng(2);
    CHECK_THROWS_AS(sbx_crossover(a, b, schema, crossover_nan_rng, nan, 2.0), std::invalid_argument);
    CHECK_THROWS_AS(polynomial_mutation(a, schema, mutation_nan_rng, nan, 20.0), std::invalid_argument);

    RandomContext zero_rng(3);
    RandomContext control_rng(3);
    CHECK(sbx_crossover(a, b, schema, zero_rng, 0.0, 2.0).values == a.values);
    CHECK(zero_rng.uniform() == control_rng.uniform());
    auto unchanged = a;
    polynomial_mutation(unchanged, schema, zero_rng, 0.0, 20.0);
    CHECK(unchanged.values == a.values);

    auto crossed = sbx_crossover(a, b, schema, control_rng, 1.0, 2.0);
    CHECK(crossed.values[0] != a.values[0]);
    auto mutated = a;
    polynomial_mutation(mutated, schema, control_rng, 1.0, 20.0);
    CHECK(mutated.values[0] != a.values[0]);
}

TEST_CASE("elite preservation keeps best candidates and population size") {
    Population p;
    for (int i = 0; i < 4; ++i) { Candidate c; c.id = static_cast<CandidateId>(i); c.objectives.push_back({"f", static_cast<double>(i), true}); c.evaluation_status = EvaluationStatus::Success; p.push_back(c); }
    FeasibilityComparator cmp;
    auto elites = p.elites(2, cmp);
    REQUIRE(elites.size() == 2);
    CHECK(elites[0].id == 3);
    CHECK(elites[1].id == 2);
}
