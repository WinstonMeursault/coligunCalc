#include <doctest/doctest.h>

#include "coilgun/optimization/variables.hpp"

#include <limits>
#include <stdexcept>

using namespace coilgun::optimization;

TEST_CASE("empty variable schemas encode and repair empty candidates") {
    const VariableSchema schema;
    const CandidateVariables empty;

    CHECK(schema.empty());
    CHECK(schema.size() == 0);
    CHECK(schema.encode(empty).values.empty());
    CHECK(schema.decode(empty).values.empty());
    CHECK(schema.repair(empty).values.empty());
}

TEST_CASE("mixed schemas repair values to deterministic bounded representations") {
    const VariableSchema schema({
        VariableSpec::continuous("length", 1.0, 2.0),
        VariableSpec::integer("turns", 2, 5),
        VariableSpec::enumeration("material", {"copper", "aluminium"}),
    });

    const auto repaired = schema.repair(CandidateVariables{{0.5, 3.5, 1.0}});
    CHECK(repaired.values == std::vector<double>{1.0, 4.0, 1.0});
    CHECK(schema.decode(schema.encode(repaired)).values == repaired.values);
}

TEST_CASE("continuous and integer repair handles infinities and NaN deterministically") {
    const VariableSchema schema({
        VariableSpec::continuous("x", -2.0, 3.0),
        VariableSpec::integer("n", -2, 3),
    });
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto repaired = schema.repair(CandidateVariables{{nan, std::numeric_limits<double>::infinity()}});

    CHECK(repaired.values[0] == doctest::Approx(-2.0));
    CHECK(repaired.values[1] == doctest::Approx(3.0));
}

TEST_CASE("enum repair rejects non-integral and out-of-range indices") {
    const VariableSchema schema({VariableSpec::enumeration("mode", {"a", "b"})});

    CHECK_THROWS_AS(schema.repair(CandidateVariables{{0.5}}), std::invalid_argument);
    CHECK_THROWS_AS(schema.repair(CandidateVariables{{2.0}}), std::out_of_range);
    CHECK_THROWS_AS(schema.repair(CandidateVariables{{-1.0}}), std::out_of_range);
}

TEST_CASE("schema construction rejects invalid bounds, enum definitions, and duplicate IDs") {
    CHECK_THROWS_AS(VariableSchema({VariableSpec::continuous("x", 2.0, 1.0)}), std::invalid_argument);
    CHECK_THROWS_AS(VariableSchema({VariableSpec::continuous("x", 0.0, 1.0),
                                    VariableSpec::integer("x", 0, 1)}), std::invalid_argument);
    CHECK_THROWS_AS(VariableSchema({VariableSpec::enumeration("mode", {})}), std::invalid_argument);
}

TEST_CASE("schema construction rejects non-finite or non-integral integer bounds") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(VariableSchema({VariableSpec{"fractional", VariableType::Integer, 0.5, 2.0, {}}}),
                    std::invalid_argument);
    CHECK_THROWS_AS(VariableSchema({VariableSpec{"infinite", VariableType::Integer, 0.0,
                                                  std::numeric_limits<double>::infinity(), {}}}),
                    std::invalid_argument);
    CHECK_THROWS_AS(VariableSchema({VariableSpec{"nan", VariableType::Integer, nan, 2.0, {}}}),
                    std::invalid_argument);
}

TEST_CASE("schema construction rejects unknown variable types") {
    const auto unknown = static_cast<VariableType>(99);
    CHECK_THROWS_AS(VariableSchema({VariableSpec{"unknown", unknown, 0.0, 1.0, {}}}),
                    std::invalid_argument);
}

TEST_CASE("continuous factory validates its bounds") {
    CHECK_THROWS_AS(VariableSpec::continuous("x", 2.0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(VariableSpec::continuous("x", -std::numeric_limits<double>::infinity(), 1.0),
                    std::invalid_argument);
}

TEST_CASE("repair rejects candidate dimensionality mismatch without mutating input or schema") {
    const VariableSchema schema({VariableSpec::continuous("x", 0.0, 1.0)});
    const CandidateVariables input{{2.0}};
    const auto specs_before = schema.variables();

    CHECK_THROWS_AS(schema.repair(CandidateVariables{}), std::invalid_argument);
    CHECK(schema.repair(input).values[0] == doctest::Approx(1.0));
    CHECK(input.values == std::vector<double>{2.0});
    CHECK(schema.variables() == specs_before);
}
