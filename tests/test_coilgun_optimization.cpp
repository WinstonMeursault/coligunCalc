#include <doctest/doctest.h>

#include "coilgun/optimization/coilgun_problem.hpp"
#include "coilgun/physics/constants.hpp"

#include <cmath>
#include <memory>

using namespace coilgun::optimization;
using coilgun::components::Armature;
using coilgun::components::DrivingCoil;
using coilgun::physics::ALUMINUM;
using coilgun::physics::COPPER;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::TriggerConfig;
using coilgun::simulation::TriggerMode;

namespace {
CoilgunOptimizationProblem::Config make_config() {
    CoilgunOptimizationProblem::Config config;
    config.coils.emplace_back(0.005, 0.010, 0.010, 12,
                              COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    config.coils.emplace_back(0.005, 0.010, 0.010, 12,
                              COPPER.resistivity_ref, 1e-6, 0.7, 0.03);
    config.armature = Armature(0.002, 0.008, 0.010,
                               ALUMINUM.resistivity_ref, ALUMINUM.density,
                               0.0, 0.005, 1, 1, 0.003);
    config.excitations = {{500.0, 500e-6, true}, {500.0, 500e-6, true}};
    config.triggers = {{TriggerMode::Position, 0.010}};
    config.dt = 1e-6;
    config.termination.max_steps = 8;
    config.termination.enable_velocity_check = false;
    config.objective_id = "muzzle_velocity";
    return config;
}
}

TEST_CASE("coilgun adapter decodes bindings and extracts simulation metrics") {
    VariableSchema schema({
        VariableSpec::integer("turns", 10, 20),
        VariableSpec::continuous("trigger", 0.005, 0.020),
        VariableSpec::continuous("voltage", 100.0, 1000.0),
    });
    auto config = make_config();
    config.bindings = {
        {"turns", CoilgunParameter::CoilTurns, 0},
        {"trigger", CoilgunParameter::TriggerValue, 0},
        {"voltage", CoilgunParameter::ExcitationVoltage, 1},
    };
    config.constraints.push_back({"efficiency", CoilgunMetric::Efficiency,
                                  ConstraintDefinition{"efficiency", ConstraintKind::Soft,
                                      ConstraintRelation::GreaterEqual, 0.0, 0.0, 1.0, 0}});

    CoilgunOptimizationProblem problem(schema, std::move(config));
    const auto result = problem.evaluate(CandidateVariables{{16.0, 0.010, 600.0}});

    REQUIRE(result.status == EvaluationStatus::Success);
    REQUIRE(result.objectives.size() == 1);
    CHECK(result.objectives.front().id == "muzzle_velocity");
    CHECK(std::isfinite(result.objectives.front().value));
    CHECK(result.constraints.size() == 1);
    CHECK(result.metadata.at("peak_current") != "");
    CHECK(result.metadata.at("efficiency") != "");
    CHECK(result.metadata.at("energy_loss") != "");
}

TEST_CASE("coilgun adapter isolates invalid candidates and falls back from GPU batch") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    CoilgunOptimizationProblem problem(schema, std::move(config));

    problem.set_gpu_batch_evaluator([](const std::vector<CandidateVariables>&,
                                       const EvaluationContext&) -> std::vector<EvaluationResult> {
        throw std::runtime_error("GPU unavailable");
    });
    const auto results = problem.evaluate_batch({CandidateVariables{{500.0}},
                                                  CandidateVariables{{700.0}}},
                                                EvaluationContext{});
    REQUIRE(results.size() == 2);
    CHECK(results[0].status == EvaluationStatus::Success);
    CHECK(results[1].status == EvaluationStatus::Success);
    CHECK(problem.last_batch_used_fallback());
}

TEST_CASE("coilgun adapter rejects malformed variables and non-finite GPU output") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    CoilgunOptimizationProblem problem(schema, std::move(config));

    const auto malformed = problem.evaluate(CandidateVariables{{NAN}});
    CHECK(malformed.status == EvaluationStatus::Invalid);

    problem.set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                       const EvaluationContext&) {
        std::vector<EvaluationResult> output(candidates.size(), EvaluationResult::success());
        output.front().objectives.push_back({"muzzle_velocity", INFINITY, true});
        return output;
    });
    const auto non_finite = problem.evaluate_batch({CandidateVariables{{500.0}},
                                                    CandidateVariables{{600.0}}},
                                                   EvaluationContext{});
    REQUIRE(non_finite.size() == 2);
    CHECK(non_finite.front().status == EvaluationStatus::Failed);
    CHECK(non_finite.front().diagnostics.front().code == "non_finite_result");
}
