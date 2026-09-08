#include <doctest/doctest.h>

#include "coilgun/optimization/coilgun_problem.hpp"
#include "coilgun/optimization/genetic_optimizer.hpp"
#include "coilgun/physics/constants.hpp"

#include <cmath>
#include <iomanip>
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
    config.coil_specs = {
        {0.005, 0.010, 0.010, 12, COPPER.resistivity_ref, 1e-6, 0.7, 0.0, false},
        {0.005, 0.010, 0.010, 12, COPPER.resistivity_ref, 1e-6, 0.7, 0.03, false},
    };
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
    REQUIRE(problem.statistics_snapshot());
    CHECK(problem.statistics_snapshot()->fallbacks == 1);
}

TEST_CASE("coilgun adapter preserves valid GPU invalid candidates alongside successes") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    CoilgunOptimizationProblem problem(schema, std::move(config));

    problem.set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                       const EvaluationContext&) {
        std::vector<EvaluationResult> output;
        output.reserve(candidates.size());
        output.push_back(EvaluationResult::invalid("candidate_invalid", "synthetic invalid candidate"));
        for (std::size_t i = 1; i < candidates.size(); ++i) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"muzzle_velocity", static_cast<double>(i), true});
            output.push_back(std::move(result));
        }
        return output;
    });

    const auto results = problem.evaluate_batch({CandidateVariables{{500.0}},
                                                  CandidateVariables{{600.0}},
                                                  CandidateVariables{{700.0}}}, {});
    REQUIRE(results.size() == 3);
    CHECK(results[0].status == EvaluationStatus::Invalid);
    CHECK(results[0].diagnostics.front().code == "candidate_invalid");
    CHECK(results[1].status == EvaluationStatus::Success);
    CHECK(results[2].status == EvaluationStatus::Success);
    CHECK_FALSE(problem.last_batch_used_fallback());
    REQUIRE(problem.statistics_snapshot());
    CHECK(problem.statistics_snapshot()->fallbacks == 0);
}

TEST_CASE("optimizer reports per-run actual GPU callback fallbacks and seed") {
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

    OptimizationConfig optimization;
    optimization.population_size = 2;
    optimization.max_generations = 1;
    optimization.random_seed = 314159;
    TerminationConfig termination;
    termination.max_generations = 1;

    const auto first = GeneticOptimizer(schema, problem, optimization, termination).optimize();
    const auto second = GeneticOptimizer(schema, problem, optimization, termination).optimize();

    CHECK(first.statistics.seed == optimization.random_seed);
    CHECK(second.statistics.seed == optimization.random_seed);
    CHECK(first.statistics.gpu_fallbacks == 1);
    CHECK(second.statistics.gpu_fallbacks == 1);
    REQUIRE(problem.statistics_snapshot());
    CHECK(problem.statistics_snapshot()->fallbacks == 2);
}

TEST_CASE("coilgun adapter dispatches injected batch evaluator through genetic optimizer") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    CoilgunOptimizationProblem problem(schema, std::move(config));

    std::size_t callback_calls = 0;
    problem.set_gpu_batch_evaluator([&callback_calls](const std::vector<CandidateVariables>& candidates,
                                                      const EvaluationContext&) {
        ++callback_calls;
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for (const auto& candidate : candidates) {
            (void)candidate;
            auto result = EvaluationResult::success();
            result.objectives.push_back({"muzzle_velocity", 1.0, true});
            results.push_back(std::move(result));
        }
        return results;
    });

    OptimizationConfig optimization;
    optimization.population_size = 2;
    optimization.max_generations = 1;
    TerminationConfig termination;
    termination.max_generations = 1;
    const auto result = GeneticOptimizer(schema, problem, optimization, termination).optimize();

    CHECK(result.termination.reason == TerminationReason::MaxGenerations);
    CHECK(callback_calls > 0);
}

TEST_CASE("coilgun adapter peak voltage includes initial excitation voltage") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    config.termination.max_steps = 0;
    CoilgunOptimizationProblem problem(schema, std::move(config));

    const auto result = problem.evaluate(CandidateVariables{{750.0}});

    REQUIRE(result.status == EvaluationStatus::Success);
    CHECK(std::stod(result.metadata.at("peak_voltage")) == doctest::Approx(750.0));
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
    CHECK(non_finite.front().status == EvaluationStatus::Success);
    CHECK(problem.last_batch_used_fallback());
}

TEST_CASE("coilgun adapter falls back for malformed GPU batch shape and constraints") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    CoilgunOptimizationProblem problem(schema, std::move(config));

    problem.set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                       const EvaluationContext&) {
        std::vector<EvaluationResult> output(candidates.size() - 1, EvaluationResult::success());
        return output;
    });
    const auto wrong_size = problem.evaluate_batch({CandidateVariables{{500.0}}, CandidateVariables{{600.0}}}, {});
    REQUIRE(wrong_size.size() == 2);
    CHECK(problem.last_batch_used_fallback());

    problem.set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                       const EvaluationContext&) {
        std::vector<EvaluationResult> output(candidates.size(), EvaluationResult::success());
        for (auto& result : output) {
            result.objectives.push_back({"muzzle_velocity", 1.0, true});
            result.constraints.push_back({"finite", ConstraintKind::Hard, ConstraintRelation::LessEqual,
                                          INFINITY, 0.0, 1.0, 0.0, 0.0, false, 0});
        }
        return output;
    });
    const auto nonfinite_constraint = problem.evaluate_batch({CandidateVariables{{500.0}}}, {});
    REQUIRE(nonfinite_constraint.size() == 1);
    CHECK(problem.last_batch_used_fallback());

    problem.set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                       const EvaluationContext&) {
        std::vector<EvaluationResult> output(candidates.size(), EvaluationResult::success());
        for (auto& result : output) result.objectives.push_back({"wrong_objective", 1.0, true});
        return output;
    });
    const auto wrong_objective = problem.evaluate_batch({CandidateVariables{{500.0}}}, {});
    REQUIRE(wrong_objective.size() == 1);
    CHECK(problem.last_batch_used_fallback());
    REQUIRE(problem.statistics_snapshot());
    CHECK(problem.statistics_snapshot()->fallbacks == 3);
}

TEST_CASE("fixed coils preserve their original inductance semantics") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto fixed = make_config();
    fixed.coils.erase(fixed.coils.begin() + 1, fixed.coils.end());
    fixed.coils[0] = DrivingCoil(0.005, 0.010, 0.010, 12,
                                 COPPER.resistivity_ref, 1e-6, 0.7, 0.0, true);
    fixed.excitations.resize(1);
    fixed.excitations[0] = {10000.0, 1e-3, true};
    fixed.triggers.clear();
    fixed.dt = 1e-5;
    fixed.termination.max_steps = 1000;
    fixed.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};

    auto explicit_specs = fixed;
    explicit_specs.coil_specs.push_back({0.005, 0.010, 0.010, 12,
                                         COPPER.resistivity_ref, 1e-6, 0.7, 0.0, true});

    CoilgunOptimizationProblem fixed_problem(schema, fixed);
    CoilgunOptimizationProblem specs_problem(schema, std::move(explicit_specs));
    const auto fixed_result = fixed_problem.evaluate(CandidateVariables{{600.0}});
    const auto specs_result = specs_problem.evaluate(CandidateVariables{{600.0}});

    REQUIRE(fixed_result.status == EvaluationStatus::Success);
    REQUIRE(specs_result.status == EvaluationStatus::Success);
    INFO("fixed=" << std::setprecision(17) << fixed_result.objectives.front().value << ", specs=" <<
         specs_result.objectives.front().value);
    CHECK(fixed_result.objectives.front().value == specs_result.objectives.front().value);
}

TEST_CASE("armature mass bindings are rejected during configuration validation") {
    VariableSchema schema({VariableSpec::continuous("mass", 0.001, 0.010)});
    auto config = make_config();
    config.bindings = {{"mass", CoilgunParameter::ArmatureMass, 0}};

    CHECK_THROWS_WITH_AS(CoilgunOptimizationProblem(schema, std::move(config)),
                         "ArmatureMass bindings are not supported",
                         std::invalid_argument);
}

TEST_CASE("coil geometry bindings require complete coil specs") {
    VariableSchema schema({VariableSpec::integer("turns", 10, 20)});
    auto config = make_config();
    config.bindings = {{"turns", CoilgunParameter::CoilTurns, 0}};

    CHECK_THROWS_WITH_AS(CoilgunOptimizationProblem(schema, std::move(config)),
                         "coil geometry bindings require coil_specs",
                         std::invalid_argument);
}
