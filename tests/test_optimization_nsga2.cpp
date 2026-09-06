#include <doctest/doctest.h>

#include "coilgun/optimization/nsga2.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace coilgun::optimization;

namespace {
Candidate candidate(CandidateId id, double first, double second, bool first_max = true,
                   bool second_max = true) {
    Candidate result;
    result.id = id;
    result.evaluation_status = EvaluationStatus::Success;
    result.objectives = {{"first", first, first_max}, {"second", second, second_max}};
    return result;
}

Candidate constrained(CandidateId id, double first, double second, double violation) {
    auto result = candidate(id, first, second);
    result.constraints.push_back({"hard", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                  violation, 0.0, 0.0, violation, violation, violation == 0.0, 0});
    return result;
}
}

TEST_CASE("NSGA-II sorts a known two-objective front and normalizes directions") {
    std::vector<Candidate> candidates{
        candidate(0, 1.0, 1.0), candidate(1, 2.0, 0.0), candidate(2, 0.0, 2.0),
        candidate(3, 0.5, 0.5),
    };

    const auto ranking = nsga2_rank(candidates);
    REQUIRE(ranking.fronts.size() == 2);
    CHECK(ranking.fronts[0] == std::vector<std::size_t>{0, 1, 2});
    CHECK(ranking.fronts[1] == std::vector<std::size_t>{3});
    CHECK(ranking.ranks[0] == 0);
    CHECK(ranking.ranks[3] == 1);

    const auto min_second = nsga2_rank(
        {candidate(0, 1.0, 1.0, true, false), candidate(1, 2.0, 0.0, true, false),
         candidate(2, 0.0, -1.0, true, false)});
    CHECK(min_second.fronts[0] == std::vector<std::size_t>{1, 2});
}

TEST_CASE("NSGA-II assigns infinite crowding at objective boundaries") {
    std::vector<Candidate> candidates{
        candidate(0, 0.0, 3.0), candidate(1, 1.0, 2.0), candidate(2, 2.0, 1.0),
        candidate(3, 3.0, 0.0),
    };
    const auto ranking = nsga2_rank(candidates);
    REQUIRE(ranking.fronts.size() == 1);
    CHECK(std::isinf(ranking.crowding_distances[ranking.fronts[0][0]]));
    CHECK(std::isinf(ranking.crowding_distances[ranking.fronts[0][3]]));
    CHECK(ranking.crowding_distances[ranking.fronts[0][1]] == doctest::Approx(4.0 / 3.0));
    CHECK(ranking.crowding_distances[ranking.fronts[0][2]] == doctest::Approx(4.0 / 3.0));
}

TEST_CASE("NSGA-II keeps duplicate objective values finite and stable") {
    std::vector<Candidate> candidates{
        candidate(10, 1.0, 1.0), candidate(11, 1.0, 1.0), candidate(12, 1.0, 1.0),
        candidate(13, 2.0, 0.0),
    };
    const auto ranking = nsga2_rank(candidates);
    REQUIRE(ranking.fronts.size() == 1);
    for (const auto distance : ranking.crowding_distances) CHECK(!std::isnan(distance));
    const auto selected = nsga2_select(candidates, {}, 2);
    REQUIRE(selected.size() == 2);
    CHECK(selected[0].id == 10);
    CHECK(selected[1].id == 11);
}

TEST_CASE("NSGA-II gives feasible candidates priority over infeasible candidates") {
    std::vector<Candidate> candidates{
        constrained(0, 100.0, 100.0, 0.0), constrained(1, 0.0, 0.0, 1.0),
        constrained(2, 1.0, 1.0, 0.2),
    };
    const auto ranking = nsga2_rank(candidates);
    REQUIRE(ranking.fronts.size() == 3);
    CHECK(ranking.fronts[0] == std::vector<std::size_t>{0});
    CHECK(ranking.fronts[1] == std::vector<std::size_t>{2});
    CHECK(ranking.fronts[2] == std::vector<std::size_t>{1});
}

TEST_CASE("NSGA-II ties infeasible candidates with equal total hard violation") {
    // Candidate 1 would dominate candidate 0 by objectives alone. Constraint
    // domination must treat equal-violation infeasible candidates as tied.
    const auto ranking = nsga2_rank({constrained(0, 1.0, 1.0, 1.0),
                                     constrained(1, 2.0, 2.0, 1.0)});
    REQUIRE(ranking.fronts.size() == 1);
    CHECK(ranking.fronts[0] == std::vector<std::size_t>{0, 1});
    CHECK(ranking.ranks[0] == 0);
    CHECK(ranking.ranks[1] == 0);
}

TEST_CASE("NSGA-II accepts failed candidates with empty objectives") {
    Candidate failed;
    failed.id = 0;
    failed.evaluation_status = EvaluationStatus::Failed;
    Candidate invalid;
    invalid.id = 1;
    invalid.evaluation_status = EvaluationStatus::Invalid;
    Candidate another_failed;
    another_failed.id = 2;
    another_failed.evaluation_status = EvaluationStatus::Failed;

    const auto failed_only = nsga2_rank({failed, invalid});
    REQUIRE(failed_only.fronts.size() == 1);
    CHECK(failed_only.fronts[0] == std::vector<std::size_t>{0, 1});

    const auto mixed = nsga2_rank({failed, candidate(3, 1.0, 2.0), invalid, another_failed});
    REQUIRE(mixed.fronts.size() == 2);
    CHECK(mixed.fronts[0] == std::vector<std::size_t>{1});
    CHECK(mixed.fronts[1] == std::vector<std::size_t>{0, 2, 3});
    for (const auto index : mixed.fronts[1])
        CHECK(std::isinf(mixed.crowding_distances[index]));
}

TEST_CASE("NSGA-II merges parents and offspring then truncates by rank and crowding") {
    std::vector<Candidate> parents{candidate(0, 0.0, 3.0), candidate(1, 1.0, 2.0)};
    std::vector<Candidate> offspring{candidate(2, 2.0, 1.0), candidate(3, 3.0, 0.0)};
    const auto selected = nsga2_select(parents, offspring, 3);
    REQUIRE(selected.size() == 3);
    CHECK(selected[0].id == 0);
    CHECK(selected[1].id == 3);
    CHECK(selected[2].id == 1);
}

TEST_CASE("NSGA-II validates fixed objective count") {
    CHECK_THROWS_AS(nsga2_rank({candidate(0, 1.0, 2.0)} ,
                               {ObjectiveDefinition{"only", true, 1.0}}), std::invalid_argument);
    CHECK_THROWS_AS(nsga2_rank({candidate(0, 1.0, 2.0), candidate(1, 1.0, 2.0, true, true)},
                               {ObjectiveDefinition{"a", true, 1.0}, ObjectiveDefinition{"b", true, 1.0},
                                ObjectiveDefinition{"c", true, 1.0}}),
                    std::invalid_argument);
}
