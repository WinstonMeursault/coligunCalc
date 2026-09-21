#include <doctest/doctest.h>
#include "coilgun/optimization/objective.hpp"
#include "coilgun/optimization/constraint.hpp"
#include "coilgun/optimization/comparator.hpp"
#include <limits>
#include <stdexcept>

using namespace coilgun::optimization;

TEST_CASE("objectives normalize scales and orient directions") {
    ObjectiveDefinition max{"speed", true, 10.0};
    ObjectiveDefinition min{"mass", false, 2.0};
    CHECK(max.normalize(25.0) == doctest::Approx(2.5));
    CHECK(max.oriented(25.0) == doctest::Approx(-2.5));
    CHECK(min.oriented(5.0) == doctest::Approx(2.5));
    ObjectiveDefinition zero{"x", true, 0};
    ObjectiveDefinition infinite{"x", true, std::numeric_limits<double>::infinity()};
    ObjectiveDefinition nan_scale{"x", true, std::numeric_limits<double>::quiet_NaN()};
    CHECK_THROWS_AS(zero.validate(), std::invalid_argument);
    CHECK_THROWS_AS(infinite.validate(), std::invalid_argument);
    CHECK_THROWS_AS(nan_scale.validate(), std::invalid_argument);
    CHECK_THROWS_AS(max.normalize(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    CHECK_THROWS_AS(max.normalize(-std::numeric_limits<double>::infinity()), std::invalid_argument);
}

TEST_CASE("constraints implement equal inequality and range violation") {
    CHECK(ConstraintDefinition{"e", ConstraintKind::Hard, ConstraintRelation::Equal, 3, 0, 2}.evaluate(4).violation == doctest::Approx(1));
    CHECK(ConstraintDefinition{"l", ConstraintKind::Hard, ConstraintRelation::LessEqual, 0, 3, 2}.evaluate(5).normalized_violation == doctest::Approx(1));
    CHECK(ConstraintDefinition{"g", ConstraintKind::Hard, ConstraintRelation::GreaterEqual, 3, 0, 2}.evaluate(2).normalized_violation == doctest::Approx(0.5));
    auto in = ConstraintDefinition{"r", ConstraintKind::Hard, ConstraintRelation::InRange, 1, 4, 1};
    CHECK(in.evaluate(3).satisfied);
    CHECK(in.evaluate(0).violation == doctest::Approx(1));
    ConstraintDefinition invalid{"r", ConstraintKind::Hard, ConstraintRelation::InRange, 4, 1, 1};
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);
}

TEST_CASE("constraints reject non-finite scales bounds and values") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS((ConstraintDefinition{"x", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                          0.0, 1.0, nan}).validate(), std::invalid_argument);
    CHECK_THROWS_AS((ConstraintDefinition{"x", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                          0.0, inf, 1.0}).validate(), std::invalid_argument);
    CHECK_THROWS_AS((ConstraintDefinition{"x", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                          -inf, 1.0, 1.0}).validate(), std::invalid_argument);
    const auto definition = ConstraintDefinition{"x", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                                  0.0, 1.0, 1.0};
    CHECK_THROWS_AS(definition.evaluate(nan), std::invalid_argument);
    CHECK_THROWS_AS(definition.evaluate(inf), std::invalid_argument);
}

TEST_CASE("hard feasibility aggregates independently from soft constraints") {
    auto h = ConstraintDefinition{"h", ConstraintKind::Hard, ConstraintRelation::LessEqual, 0, 1, 2}.evaluate(3);
    auto s = ConstraintDefinition{"s", ConstraintKind::Soft, ConstraintRelation::LessEqual, 0, 1, 1}.evaluate(4);
    CHECK(aggregate_normalized_violation({h, s}, ConstraintKind::Hard) == doctest::Approx(1));
    CHECK(aggregate_normalized_violation({h, s}, ConstraintKind::Soft) == doctest::Approx(3));
    CHECK_FALSE(is_feasible({h, s}));
    CHECK(is_feasible({s}));
}

static Candidate candidate(double objective, std::vector<ConstraintReport> constraints = {}) {
    Candidate c; c.objectives.push_back(ObjectiveValue{"score", objective, true}); c.constraints = std::move(constraints);
    c.evaluation_status = EvaluationStatus::Success;
    return c;
}

TEST_CASE("comparators prioritize feasibility and then objectives") {
    auto bad = candidate(100, {ConstraintDefinition{"h", ConstraintKind::Hard, ConstraintRelation::LessEqual, 0, 1, 1}.evaluate(2)});
    auto good = candidate(1);
    CHECK(FeasibilityComparator{}.better(good, bad));
    CHECK(FeasibilityComparator{}.better(bad, candidate(0, {bad.constraints.front()})));
    auto soft_bad = candidate(0, {ConstraintDefinition{"s", ConstraintKind::Soft, ConstraintRelation::LessEqual, 0, 1, 1}.evaluate(3)});
    CHECK(FeasibilityComparator{FeasibilityStrategy::Lexicographic}.better(good, soft_bad));
    CHECK(FeasibilityComparator{FeasibilityStrategy::Penalty, 2}.better(good, soft_bad));
}

TEST_CASE("lexicographic comparison honors per-constraint priorities") {
    const auto high_lhs = ConstraintDefinition{"high", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                                0.0, 0.0, 1.0, 1};
    const auto low_lhs = ConstraintDefinition{"low", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                               0.0, 0.0, 1.0, 2};
    const auto high_rhs = ConstraintDefinition{"high", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                                0.0, 0.0, 1.0, 1};
    const auto low_rhs = ConstraintDefinition{"low", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                               0.0, 0.0, 1.0, 2};
    const auto lhs = candidate(0.0, {high_lhs.evaluate(1.0), low_lhs.evaluate(100.0)});
    const auto rhs = candidate(0.0, {high_rhs.evaluate(2.0), low_rhs.evaluate(0.0)});
    CHECK(FeasibilityComparator{FeasibilityStrategy::Lexicographic}.better(rhs, lhs));

    const auto soft_high_lhs = ConstraintDefinition{"soft-high", ConstraintKind::Soft,
                                                     ConstraintRelation::LessEqual, 0.0, 0.0, 1.0, 1};
    const auto soft_low_lhs = ConstraintDefinition{"soft-low", ConstraintKind::Soft,
                                                    ConstraintRelation::LessEqual, 0.0, 0.0, 1.0, 2};
    const auto soft_high_rhs = ConstraintDefinition{"soft-high", ConstraintKind::Soft,
                                                     ConstraintRelation::LessEqual, 0.0, 0.0, 1.0, 1};
    const auto soft_low_rhs = ConstraintDefinition{"soft-low", ConstraintKind::Soft,
                                                    ConstraintRelation::LessEqual, 0.0, 0.0, 1.0, 2};
    const auto soft_lhs = candidate(0.0, {soft_high_lhs.evaluate(1.0), soft_low_lhs.evaluate(100.0)});
    const auto soft_rhs = candidate(0.0, {soft_high_rhs.evaluate(2.0), soft_low_rhs.evaluate(0.0)});
    CHECK(FeasibilityComparator{FeasibilityStrategy::Lexicographic}.better(soft_lhs, soft_rhs));

    const auto hard_bad = candidate(0.0, {ConstraintDefinition{"hard", ConstraintKind::Hard,
                                                                ConstraintRelation::LessEqual, 0.0, 0.0, 1.0, 2}
                                             .evaluate(1.0)});
    const auto soft_bad = candidate(0.0, {ConstraintDefinition{"soft", ConstraintKind::Soft,
                                                                ConstraintRelation::LessEqual, 0.0, 0.0, 1.0, 1}
                                             .evaluate(1.0)});
    CHECK(FeasibilityComparator{FeasibilityStrategy::Lexicographic}.better(soft_bad, hard_bad));
}

TEST_CASE("comparators reject invalid penalty weights") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS((FeasibilityComparator{FeasibilityStrategy::Penalty, -1.0}), std::invalid_argument);
    CHECK_THROWS_AS((FeasibilityComparator{FeasibilityStrategy::Penalty, nan}), std::invalid_argument);
    CHECK_THROWS_AS((FeasibilityComparator{FeasibilityStrategy::Penalty, inf}), std::invalid_argument);
}

TEST_CASE("comparators rank successful evaluations ahead of invalid and failed candidates") {
    auto successful = candidate(-100.0);
    successful.evaluation_status = EvaluationStatus::Success;
    auto invalid = candidate(100.0);
    invalid.evaluation_status = EvaluationStatus::Invalid;
    auto failed = candidate(100.0);
    failed.evaluation_status = EvaluationStatus::Failed;

    for (const auto strategy : {FeasibilityStrategy::FeasibilityFirst,
                                 FeasibilityStrategy::Penalty,
                                 FeasibilityStrategy::Lexicographic}) {
        const FeasibilityComparator comparator{strategy, 2.0};
        CHECK(comparator.better(successful, invalid));
        CHECK(comparator.better(successful, failed));
    }
}

TEST_CASE("finite hard violations outrank non-success status failures") {
    auto successful_bad = candidate(-100.0, {ConstraintDefinition{"hard", ConstraintKind::Hard,
                                                                   ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                                .evaluate(1.0)});
    successful_bad.evaluation_status = EvaluationStatus::Success;
    auto failed_empty = candidate(100.0);
    failed_empty.evaluation_status = EvaluationStatus::Failed;
    for (const auto strategy : {FeasibilityStrategy::FeasibilityFirst,
                                 FeasibilityStrategy::Penalty,
                                 FeasibilityStrategy::Lexicographic}) {
        CHECK(FeasibilityComparator{strategy, 2.0}.better(successful_bad, failed_empty));
    }
}

TEST_CASE("non-success statuses have deterministic violation ordering") {
    auto invalid = candidate(0.0);
    invalid.evaluation_status = EvaluationStatus::Invalid;
    auto failed = candidate(0.0);
    failed.evaluation_status = EvaluationStatus::Failed;
    auto unevaluated = candidate(0.0);
    unevaluated.evaluation_status = EvaluationStatus::Unevaluated;
    for (const auto strategy : {FeasibilityStrategy::FeasibilityFirst,
                                 FeasibilityStrategy::Penalty,
                                 FeasibilityStrategy::Lexicographic}) {
        const FeasibilityComparator comparator{strategy, 2.0};
        CHECK(comparator.better(invalid, failed));
        CHECK(comparator.better(failed, unevaluated));
    }
}

TEST_CASE("penalty comparator adds soft violation to a maximizing objective") {
    const auto objective_advantage = candidate(100.0, {ConstraintDefinition{"soft", ConstraintKind::Soft,
                                                                              ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                                           .evaluate(0.5)});
    const auto no_penalty = candidate(99.0);
    CHECK(FeasibilityComparator{FeasibilityStrategy::Penalty, 1.0}.better(objective_advantage, no_penalty));
}

TEST_CASE("penalty comparator adds soft violation to a minimizing objective") {
    auto objective_advantage = candidate(0.0, {ConstraintDefinition{"soft", ConstraintKind::Soft,
                                                                        ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                                     .evaluate(2.0)});
    objective_advantage.objectives.front().maximize = false;
    auto no_penalty = candidate(100.0);
    no_penalty.objectives.front().maximize = false;
    CHECK(FeasibilityComparator{FeasibilityStrategy::Penalty, 1.0}.better(objective_advantage, no_penalty));
}

TEST_CASE("penalty comparator honors explicit objective scaling") {
    const auto scaled_a = candidate(100.0, {ConstraintDefinition{"soft", ConstraintKind::Soft,
                                                                   ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                                .evaluate(1.0)});
    const auto scaled_b = candidate(98.0);
    const ObjectiveDefinition definition{"score", true, 100.0};
    const FeasibilityComparator comparator{FeasibilityStrategy::Penalty, 1.0};
    CHECK(comparator.better(scaled_b, scaled_a, definition));
    CHECK(comparator.better(scaled_a, scaled_b, ObjectiveDefinition{"score", true, 1.0}));
}

TEST_CASE("zero penalty weight disables soft penalties") {
    const auto soft_bad = candidate(2.0, {ConstraintDefinition{"soft", ConstraintKind::Soft,
                                                                 ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                              .evaluate(100.0)});
    const auto soft_good = candidate(1.0);
    CHECK(FeasibilityComparator{FeasibilityStrategy::Penalty, 0.0}.better(soft_bad, soft_good));
}

TEST_CASE("penalty keeps hard feasibility ahead of soft penalties and objectives") {
    const auto feasible = candidate(-1000.0, {ConstraintDefinition{"soft", ConstraintKind::Soft,
                                                                      ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                                   .evaluate(100.0)});
    const auto hard_infeasible = candidate(1000.0, {ConstraintDefinition{"hard", ConstraintKind::Hard,
                                                                           ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                                        .evaluate(1.0)});
    const FeasibilityComparator comparator{FeasibilityStrategy::Penalty, 1.0};
    CHECK(comparator.better(feasible, hard_infeasible));

    const auto less_hard_violation = candidate(-1000.0, {ConstraintDefinition{"hard", ConstraintKind::Hard,
                                                                                ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                                             .evaluate(1.0),
                                                           ConstraintDefinition{"soft", ConstraintKind::Soft,
                                                                                ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                                             .evaluate(100.0)});
    const auto more_hard_violation = candidate(1000.0, {ConstraintDefinition{"hard", ConstraintKind::Hard,
                                                                                ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                                             .evaluate(2.0)});
    CHECK(comparator.better(less_hard_violation, more_hard_violation));
}

TEST_CASE("penalty constraint-only comparison ignores soft violations") {
    auto lhs = candidate(0.0, {ConstraintDefinition{"soft", ConstraintKind::Soft,
                                                       ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                    .evaluate(100.0)});
    auto rhs = candidate(0.0);
    lhs.objectives.clear();
    rhs.objectives.clear();
    const FeasibilityComparator comparator{FeasibilityStrategy::Penalty, 10.0};
    CHECK(comparator.compare(lhs, rhs) == 0);
    CHECK(comparator.compare(rhs, lhs) == 0);
}

TEST_CASE("lexicographic status ordering precedes non-success constraint priorities") {
    auto invalid = candidate(0.0, {ConstraintDefinition{"hard", ConstraintKind::Hard,
                                                         ConstraintRelation::LessEqual, 0.0, 0.0, 1.0}
                                      .evaluate(2.0)});
    invalid.evaluation_status = EvaluationStatus::Invalid;
    auto failed = candidate(0.0);
    failed.evaluation_status = EvaluationStatus::Failed;

    const FeasibilityComparator comparator{FeasibilityStrategy::Lexicographic};
    CHECK(comparator.better(invalid, failed));
    CHECK_FALSE(comparator.better(failed, invalid));
}

TEST_CASE("lexicographic compares constraints for candidates with the same non-success status") {
    const auto hard = ConstraintDefinition{"hard", ConstraintKind::Hard,
                                           ConstraintRelation::LessEqual, 0.0, 0.0, 1.0, 1};
    const auto soft = ConstraintDefinition{"soft", ConstraintKind::Soft,
                                           ConstraintRelation::LessEqual, 0.0, 0.0, 1.0, 1};

    auto hard_better = candidate(100.0, {hard.evaluate(1.0), soft.evaluate(2.0)});
    auto hard_worse = candidate(0.0, {hard.evaluate(2.0), soft.evaluate(1.0)});
    hard_better.evaluation_status = EvaluationStatus::Failed;
    hard_worse.evaluation_status = EvaluationStatus::Failed;

    const FeasibilityComparator comparator{FeasibilityStrategy::Lexicographic};
    CHECK(comparator.better(hard_better, hard_worse));

    auto soft_better = candidate(100.0, {hard.evaluate(1.0), soft.evaluate(1.0)});
    auto soft_worse = candidate(0.0, {hard.evaluate(1.0), soft.evaluate(2.0)});
    soft_better.evaluation_status = EvaluationStatus::Failed;
    soft_worse.evaluation_status = EvaluationStatus::Failed;
    CHECK(comparator.better(soft_better, soft_worse));
}

TEST_CASE("same-status non-success candidates use raw oriented objective tie breaks") {
    for (const auto status : {EvaluationStatus::Invalid, EvaluationStatus::Failed,
                              EvaluationStatus::Unevaluated}) {
        for (const auto strategy : {FeasibilityStrategy::FeasibilityFirst,
                                    FeasibilityStrategy::Penalty,
                                    FeasibilityStrategy::Lexicographic}) {
            auto objectively_better = candidate(10.0);
            objectively_better.evaluation_status = status;
            if (strategy == FeasibilityStrategy::Penalty) {
                objectively_better.constraints.push_back(
                    ConstraintDefinition{"soft", ConstraintKind::Soft, ConstraintRelation::LessEqual,
                                         0.0, 0.0, 1.0}
                        .evaluate(100.0));
            }

            auto objectively_worse = candidate(2.0);
            objectively_worse.evaluation_status = status;

            const FeasibilityComparator comparator{strategy, 100.0};
            CHECK(comparator.better(objectively_better, objectively_worse));
            CHECK_FALSE(comparator.better(objectively_worse, objectively_better));
        }
    }
}
