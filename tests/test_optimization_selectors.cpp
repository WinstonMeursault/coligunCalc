#include <doctest/doctest.h>

#include "coilgun/optimization/result.hpp"
#include "coilgun/optimization/selectors.hpp"

#include <functional>
#include <limits>
#include <stdexcept>

using namespace coilgun::optimization;

namespace {
Candidate candidate(CandidateId id, std::initializer_list<ObjectiveValue> objectives,
                    double violation = 0.0) {
    Candidate value;
    value.id = id;
    value.evaluation_status = EvaluationStatus::Success;
    value.objectives.assign(objectives.begin(), objectives.end());
    if (violation > 0.0) {
        value.constraints.push_back({"limit", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                     violation, 0.0, 1.0, violation, violation, false, 0});
    }
    return value;
}

class LastCandidateSelector final : public RepresentativeSelector {
public:
    std::optional<Candidate> select(const OptimizationResult& result) const override {
        if (result.pareto_front.empty()) return std::nullopt;
        return result.pareto_front.back();
    }
};

OptimizationResult sample_result() {
    OptimizationResult result;
    result.pareto_front = {
        candidate(1, {{"velocity", 10.0, true}, {"cost", 8.0, false}}),
        candidate(2, {{"velocity", 14.0, true}, {"cost", 12.0, false}}),
        candidate(3, {{"velocity", 12.0, true}, {"cost", 5.0, false}}),
    };
    return result;
}
}

TEST_CASE("representative selection is explicit and does not mutate the Pareto front") {
    const auto result = sample_result();
    const auto before_ids = std::vector<CandidateId>{1, 2, 3};
    const auto selected = result.select_representative(MaxObjective{"velocity"});
    REQUIRE(selected.has_value());
    CHECK(selected->id == 2);
    CHECK(result.pareto_front.size() == before_ids.size());
    for (std::size_t i = 0; i < before_ids.size(); ++i) CHECK(result.pareto_front[i].id == before_ids[i]);
    CHECK(result.best_by_objective.empty());
}

TEST_CASE("objective selectors honor direction and preserve first candidate on ties") {
    const auto result = sample_result();
    REQUIRE(result.select_representative(MaxObjective{"cost"})->id == 3);

    OptimizationResult ties;
    ties.pareto_front = {
        candidate(7, {{"score", 2.0, true}}),
        candidate(8, {{"score", 2.0, true}}),
    };
    CHECK(ties.select_representative(MaxObjective{"score"})->id == 7);
}

TEST_CASE("constraint selector minimizes normalized violation") {
    const auto result = sample_result();
    CHECK(result.select_representative(MinConstraintViolationMargin{})->id == 1);

    OptimizationResult infeasible;
    infeasible.pareto_front = {
        candidate(1, {{"score", 5.0, true}}, 0.8),
        candidate(2, {{"score", 4.0, true}}, 0.2),
    };
    CHECK(infeasible.select_representative(MinConstraintViolationMargin{})->id == 2);

    OptimizationResult soft_constraints;
    soft_constraints.pareto_front = {
        candidate(3, {{"score", 5.0, true}}),
        candidate(4, {{"score", 4.0, true}}),
    };
    soft_constraints.pareto_front[0].constraints.push_back(
        {"temperature", ConstraintKind::Soft, ConstraintRelation::LessEqual,
         1.0, 0.0, 0.0, 1.0, 0.4, false, 0});
    soft_constraints.pareto_front[1].constraints.push_back(
        {"temperature", ConstraintKind::Soft, ConstraintRelation::LessEqual,
         1.0, 0.0, 0.0, 1.0, 0.1, false, 0});
    CHECK(soft_constraints.select_representative(MinConstraintViolationMargin{})->id == 4);
}

TEST_CASE("ideal point distance and weighted score normalize objective ranges") {
    const auto result = sample_result();
    CHECK(result.select_representative(IdealPointDistance{})->id == 3);

    const auto weighted = result.select_representative(WeightedScore{{0.75, 0.25}});
    REQUIRE(weighted.has_value());
    CHECK(weighted->id == 2);
}

TEST_CASE("lexicographic selector uses ordered objective directions") {
    const auto result = sample_result();
    CHECK(result.select_representative(LexicographicObjectives{{"velocity", "cost"}})->id == 2);
    CHECK(result.select_representative(LexicographicObjectives{{"cost", "velocity"}})->id == 3);
}

TEST_CASE("custom selector callback can choose a representative") {
    const auto result = sample_result();
    CallbackSelector selector([](const OptimizationResult& value) {
        return value.pareto_front.at(1);
    });
    CHECK(result.select_representative(selector)->id == 2);
    CHECK(result.select_representative(LastCandidateSelector{})->id == 3);
}

TEST_CASE("single-objective fronts use the same result selection API") {
    OptimizationResult result;
    result.pareto_front = {candidate(9, {{"score", 42.0, true}})};
    CHECK(result.select_representative(MaxObjective{"score"})->id == 9);
    CHECK(result.select_representative(IdealPointDistance{})->id == 9);
    CHECK(result.select_representative(WeightedScore{{1.0}})->id == 9);
    CHECK(result.select_representative(LexicographicObjectives{{"score"}})->id == 9);
}

TEST_CASE("empty Pareto fronts return an empty representative") {
    OptimizationResult result;
    CHECK_FALSE(result.select_representative(MaxObjective{"score"}).has_value());
    CHECK_FALSE(result.select_representative(MinConstraintViolationMargin{}).has_value());
    CHECK_FALSE(result.select_representative(IdealPointDistance{0.0}).has_value());
    CHECK_FALSE(result.select_representative(WeightedScore{{}}).has_value());
    CHECK_FALSE(result.select_representative(LexicographicObjectives{{}}).has_value());
    CHECK_FALSE(result.select_representative(CallbackSelector{}).has_value());
}

TEST_CASE("invalid selector configuration is reported") {
    const auto result = sample_result();
    CHECK_THROWS_AS(result.select_representative(MaxObjective{"missing"}), std::invalid_argument);
    CHECK_THROWS_AS(result.select_representative(WeightedScore{{}}), std::invalid_argument);
}

TEST_CASE("normalized selectors reject reordered Pareto objective IDs") {
    OptimizationResult result;
    result.pareto_front = {
        candidate(1, {{"velocity", 10.0, true}, {"cost", 8.0, false}}),
        candidate(2, {{"cost", 5.0, false}, {"velocity", 14.0, true}}),
    };

    CHECK_THROWS_AS(result.select_representative(IdealPointDistance{}), std::invalid_argument);
    CHECK_THROWS_AS(result.select_representative(WeightedScore{{0.5, 0.5}}), std::invalid_argument);
}

TEST_CASE("normalized selectors reject Pareto objective direction changes") {
    OptimizationResult result;
    result.pareto_front = {
        candidate(1, {{"velocity", 10.0, true}, {"cost", 8.0, false}}),
        candidate(2, {{"velocity", 14.0, false}, {"cost", 5.0, false}}),
    };

    CHECK_THROWS_AS(result.select_representative(IdealPointDistance{}), std::invalid_argument);
    CHECK_THROWS_AS(result.select_representative(WeightedScore{{0.5, 0.5}}), std::invalid_argument);
}

TEST_CASE("normalized selectors handle opposite finite extrema without overflow") {
    OptimizationResult result;
    result.pareto_front = {
        candidate(1, {{"score", -std::numeric_limits<double>::max(), true}}),
        candidate(2, {{"score", std::numeric_limits<double>::max(), true}}),
    };

    CHECK(result.select_representative(IdealPointDistance{})->id == 2);
    CHECK(result.select_representative(WeightedScore{{1.0}})->id == 2);
}
