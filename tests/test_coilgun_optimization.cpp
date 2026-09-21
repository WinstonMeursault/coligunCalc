#include <doctest/doctest.h>

#include "coilgun/optimization/coilgun_problem.hpp"
#include "coilgun/optimization/genetic_optimizer.hpp"
#include "coilgun/physics/constants.hpp"

#include <cmath>
#include <barrier>
#include <iomanip>
#include <memory>
#include <mutex>
#include <thread>

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

TEST_CASE("coilgun adapter rejects maximum-temperature constraints when thermal simulation is disabled") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.constraints.push_back({"maximum_temperature", CoilgunMetric::MaximumTemperature,
                                  ConstraintDefinition{"maximum_temperature", ConstraintKind::Hard,
                                      ConstraintRelation::LessEqual, 350.0, 0.0, 1000.0, 0}});

    CHECK_THROWS_WITH_AS(CoilgunOptimizationProblem(schema, std::move(config)),
                         "maximum-temperature constraints require thermal simulation to be enabled",
                         std::invalid_argument);
}

TEST_CASE("coilgun adapter checks maximum-temperature constraints after other constraints") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.constraints.push_back({"efficiency", CoilgunMetric::Efficiency,
                                  ConstraintDefinition{"efficiency", ConstraintKind::Soft,
                                      ConstraintRelation::GreaterEqual, 0.0, 0.0, 1.0, 0}});
    config.constraints.push_back({"maximum_temperature", CoilgunMetric::MaximumTemperature,
                                  ConstraintDefinition{"maximum_temperature", ConstraintKind::Hard,
                                      ConstraintRelation::LessEqual, 350.0, 0.0, 1000.0, 0}});

    CHECK_THROWS_WITH_AS(CoilgunOptimizationProblem(schema, std::move(config)),
                         "maximum-temperature constraints require thermal simulation to be enabled",
                         std::invalid_argument);
}

TEST_CASE("coilgun adapter accepts maximum-temperature constraints with thermal simulation enabled") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.enable_thermal = true;
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    config.constraints.push_back({"maximum_temperature", CoilgunMetric::MaximumTemperature,
                                  ConstraintDefinition{"maximum_temperature", ConstraintKind::Hard,
                                      ConstraintRelation::LessEqual, 350.0, 0.0, 1000.0, 0}});

    CoilgunOptimizationProblem problem(schema, std::move(config));
    const auto result = problem.evaluate(CandidateVariables{{500.0}});

    REQUIRE(result.status == EvaluationStatus::Success);
    REQUIRE(result.constraints.size() == 1);
    CHECK(result.constraints.front().id == "maximum_temperature");
    CHECK(std::isfinite(result.constraints.front().value));
}

TEST_CASE("coilgun adapter accepts non-temperature constraints without thermal simulation") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    config.constraints.push_back({"efficiency", CoilgunMetric::Efficiency,
                                  ConstraintDefinition{"efficiency", ConstraintKind::Soft,
                                      ConstraintRelation::GreaterEqual, 0.0, 0.0, 1.0, 0}});

    CoilgunOptimizationProblem problem(schema, std::move(config));
    const auto result = problem.evaluate(CandidateVariables{{500.0}});

    REQUIRE(result.status == EvaluationStatus::Success);
    REQUIRE(result.constraints.size() == 1);
    CHECK(result.constraints.front().id == "efficiency");
}

TEST_CASE("coil turns binding changes the physical objective") {
    VariableSchema schema({VariableSpec::integer("turns", 8, 24)});
    auto config = make_config();
    config.coil_specs = {
        {0.005, 0.010, 0.010, 12, COPPER.resistivity_ref, 1e-6, 0.7, 0.0, false},
        {0.005, 0.010, 0.010, 12, COPPER.resistivity_ref, 1e-6, 0.7, 0.03, false},
    };
    config.bindings = {{"turns", CoilgunParameter::CoilTurns, 0}};
    CoilgunOptimizationProblem problem(schema, std::move(config));

    const auto low_turns = problem.evaluate(CandidateVariables{{8.0}});
    const auto high_turns = problem.evaluate(CandidateVariables{{24.0}});

    REQUIRE(low_turns.status == EvaluationStatus::Success);
    REQUIRE(high_turns.status == EvaluationStatus::Success);
    const double low_velocity = low_turns.objectives.front().value;
    const double high_velocity = high_turns.objectives.front().value;
    REQUIRE(std::isfinite(low_velocity));
    REQUIRE(std::isfinite(high_velocity));
    CHECK(std::abs(low_velocity - high_velocity) > 0.0);
}

TEST_CASE("trigger value binding changes the physical objective") {
    VariableSchema schema({VariableSpec::continuous("trigger", 0.0, 100e-6)});
    auto config = make_config();
    config.triggers[0] = {TriggerMode::TimeDelay, 0.0};
    config.bindings = {{"trigger", CoilgunParameter::TriggerValue, 0}};
    CoilgunOptimizationProblem problem(schema, std::move(config));

    const auto immediate = problem.evaluate(CandidateVariables{{0.0}});
    const auto delayed = problem.evaluate(CandidateVariables{{100e-6}});

    REQUIRE(immediate.status == EvaluationStatus::Success);
    REQUIRE(delayed.status == EvaluationStatus::Success);
    const double immediate_velocity = immediate.objectives.front().value;
    const double delayed_velocity = delayed.objectives.front().value;
    REQUIRE(std::isfinite(immediate_velocity));
    REQUIRE(std::isfinite(delayed_velocity));
    CHECK(std::abs(immediate_velocity - delayed_velocity) > 0.0);
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

TEST_CASE("problem-only optimizer preserves the coilgun batch callback path") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    CoilgunOptimizationProblem problem(schema, std::move(config));
    std::size_t callback_calls = 0;
    problem.set_gpu_batch_evaluator([&callback_calls](
        const std::vector<CandidateVariables>& candidates,
        const EvaluationContext&) {
        ++callback_calls;
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for (const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"muzzle_velocity", candidate.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    });

    OptimizationConfig optimizer_config;
    optimizer_config.population_size = 4;
    optimizer_config.max_generations = 1;
    optimizer_config.crossover_rate = 0.0;
    optimizer_config.mutation_rate = 0.0;
    const auto result = GeneticOptimizer(problem, optimizer_config).optimize();

    CHECK(result.termination.reason == TerminationReason::MaxGenerations);
    CHECK(callback_calls == 1);
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

TEST_CASE("coilgun GPU batch fallback is attributed to the active run collector") {
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
    auto collector = std::make_shared<EvaluationStatisticsCollector>();
    const auto results = problem.evaluate_batch({CandidateVariables{{500.0}}},
                                                EvaluationContext{9, false, collector});

    REQUIRE(results.size() == 1);
    CHECK(results.front().status == EvaluationStatus::Success);
    CHECK(collector->snapshot().fallbacks == 1);
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

TEST_CASE("coilgun callback validation enforces the declared objective and constraint schema") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    config.constraints = {
        {"efficiency", CoilgunMetric::Efficiency,
         ConstraintDefinition{"efficiency", ConstraintKind::Soft,
                              ConstraintRelation::GreaterEqual, 0.0, 0.0, 1.0, 2}},
        {"energy_loss", CoilgunMetric::EnergyLoss,
         ConstraintDefinition{"energy_loss", ConstraintKind::Hard,
                              ConstraintRelation::LessEqual, 0.0, 1.0, 2.0, 3}},
    };
    CoilgunOptimizationProblem problem(schema, std::move(config));
    const CandidateVariables candidate{{500.0}};
    const auto valid = problem.evaluate(candidate);
    REQUIRE(valid.status == EvaluationStatus::Success);

    const auto run_malformed = [&](auto mutate) {
        problem.set_gpu_batch_evaluator([valid, mutate](const std::vector<CandidateVariables>& candidates,
                                                        const EvaluationContext&) {
            auto result = valid;
            mutate(result);
            return std::vector<EvaluationResult>(candidates.size(), std::move(result));
        });
        const auto results = problem.evaluate_batch({candidate}, {});
        REQUIRE(results.size() == 1);
        CHECK(results.front().status == EvaluationStatus::Success);
        CHECK(problem.last_batch_used_fallback());
    };

    SUBCASE("missing constraint") {
        run_malformed([](EvaluationResult& result) { result.constraints.pop_back(); });
    }
    SUBCASE("extra constraint") {
        run_malformed([](EvaluationResult& result) { result.constraints.push_back(result.constraints.back()); });
    }
    SUBCASE("constraint order") {
        run_malformed([](EvaluationResult& result) { std::swap(result.constraints[0], result.constraints[1]); });
    }
    SUBCASE("constraint identity and declaration") {
        run_malformed([](EvaluationResult& result) { result.constraints[0].id = "wrong"; });
        run_malformed([](EvaluationResult& result) { result.constraints[0].kind = ConstraintKind::Hard; });
        run_malformed([](EvaluationResult& result) {
            result.constraints[0].relation = ConstraintRelation::Equal;
        });
        run_malformed([](EvaluationResult& result) { result.constraints[0].lower_bound += 1.0; });
        run_malformed([](EvaluationResult& result) { result.constraints[0].priority += 1; });
    }
    SUBCASE("derived fields and objective direction") {
        run_malformed([](EvaluationResult& result) { result.constraints[0].normalized_violation += 1.0; });
        run_malformed([](EvaluationResult& result) { result.constraints[0].normalized_violation = INFINITY; });
        run_malformed([](EvaluationResult& result) { result.objectives[0].maximize = false; });
    }
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

TEST_CASE("coilgun problem cache identity isolates physical configurations and repeats hit") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto first_config = make_config();
    first_config.coils.erase(first_config.coils.begin() + 1, first_config.coils.end());
    first_config.excitations.resize(1);
    first_config.triggers.clear();
    first_config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    auto second_config = first_config;
    second_config.excitations[0].initial_voltage = 650.0;

    auto first = std::make_shared<CoilgunOptimizationProblem>(schema, std::move(first_config));
    auto second = std::make_shared<CoilgunOptimizationProblem>(schema, std::move(second_config));
    CHECK(first->cache_identity().namespace_id != second->cache_identity().namespace_id);

    auto cache = std::make_shared<InMemoryEvaluationCache>();
    CachedBatchEvaluator first_cached(first, cache);
    CachedBatchEvaluator second_cached(second, cache);
    const CandidateVariables candidate{{500.0}};
    REQUIRE(first_cached.evaluate_batch({candidate}, {}).front().status == EvaluationStatus::Success);
    REQUIRE(second_cached.evaluate_batch({candidate}, {}).front().status == EvaluationStatus::Success);
    CHECK(first_cached.statistics().cache_hits == 0);
    CHECK(second_cached.statistics().cache_hits == 0);
    REQUIRE(first_cached.evaluate_batch({candidate}, {}).front().status == EvaluationStatus::Success);
    CHECK(first_cached.statistics().cache_hits == 1);
}

TEST_CASE("callback-backed coilgun problem identities change per instance and replacement") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto make_problem = [&] {
        auto config = make_config();
        config.coils.erase(config.coils.begin() + 1, config.coils.end());
        config.excitations.resize(1);
        config.triggers.clear();
        config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
        return std::make_shared<CoilgunOptimizationProblem>(schema, std::move(config));
    };
    auto first = make_problem();
    auto second = make_problem();
    first->set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                     const EvaluationContext&) {
        return std::vector<EvaluationResult>(candidates.size(), EvaluationResult::success());
    });
    second->set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                      const EvaluationContext&) {
        return std::vector<EvaluationResult>(candidates.size(), EvaluationResult::success());
    });
    CHECK(first->cache_identity().namespace_id != second->cache_identity().namespace_id);
    const auto callback_identity = first->cache_identity().namespace_id;
    first->set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                      const EvaluationContext&) {
        return std::vector<EvaluationResult>(candidates.size(), EvaluationResult::success());
    });
    CHECK(first->cache_identity().namespace_id != callback_identity);
    first->clear_gpu_batch_evaluator();
    CHECK(first->cache_identity().namespace_id != callback_identity);
}

TEST_CASE("coilgun callback replacement keeps a shared cache generation paired") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    auto problem = std::make_shared<CoilgunOptimizationProblem>(schema, std::move(config));

    std::barrier callback_entered{2};
    std::barrier callback_release{2};
    problem->set_gpu_batch_evaluator([&](const std::vector<CandidateVariables>& candidates,
                                         const EvaluationContext&) {
        callback_entered.arrive_and_wait();
        callback_release.arrive_and_wait();
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for ([[maybe_unused]] const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"muzzle_velocity", 1.0, true});
            results.push_back(std::move(result));
        }
        return results;
    });
    const auto old_identity = problem->cache_identity();
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    CachedBatchEvaluator cached{problem, cache};
    const CandidateVariables candidate{{500.0}};
    const EvaluationContext context{31, false};
    std::vector<EvaluationResult> old_results;
    std::thread worker([&] { old_results = cached.evaluate_batch({candidate}, context); });
    callback_entered.arrive_and_wait();

    problem->set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                        const EvaluationContext&) {
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for ([[maybe_unused]] const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"muzzle_velocity", 2.0, true});
            results.push_back(std::move(result));
        }
        return results;
    });
    callback_release.arrive_and_wait();
    worker.join();

    REQUIRE(old_results.size() == 1);
    CHECK(old_results.front().objectives.front().value == 1.0);
    const auto old_cached = cache->get(make_cache_key(old_identity, candidate, context));
    REQUIRE(old_cached);
    CHECK(old_cached->objectives.front().value == 1.0);

    const auto new_results = cached.evaluate_batch({candidate}, context);
    REQUIRE(new_results.size() == 1);
    CHECK(new_results.front().objectives.front().value == 2.0);
    CHECK(cached.statistics().cache_hits == 0);
    const auto repeated_results = cached.evaluate_batch({candidate}, context);
    CHECK(repeated_results.front().objectives.front().value == 2.0);
    CHECK(cached.statistics().cache_hits == 1);
}

TEST_CASE("unmanaged coilgun problems reject escaping snapshots") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    CoilgunOptimizationProblem problem(schema, std::move(config));
    problem.set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                       const EvaluationContext&) {
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for ([[maybe_unused]] const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"muzzle_velocity", 3.0, true});
            results.push_back(std::move(result));
        }
        return results;
    });

    CHECK_THROWS_WITH_AS(problem.evaluation_snapshot(),
                         "BatchEvaluator::evaluation_snapshot requires shared ownership",
                         std::logic_error);
}

TEST_CASE("shared coilgun callback snapshots survive external owner reset") {
    VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 1000.0)});
    auto config = make_config();
    config.coils.erase(config.coils.begin() + 1, config.coils.end());
    config.excitations.resize(1);
    config.triggers.clear();
    config.bindings = {{"voltage", CoilgunParameter::ExcitationVoltage, 0}};
    auto problem = std::make_shared<CoilgunOptimizationProblem>(schema, std::move(config));
    problem->set_gpu_batch_evaluator([](const std::vector<CandidateVariables>& candidates,
                                        const EvaluationContext&) {
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for ([[maybe_unused]] const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"muzzle_velocity", 3.0, true});
            results.push_back(std::move(result));
        }
        return results;
    });

    auto snapshot = problem->evaluation_snapshot();
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    CachedBatchEvaluator cached(problem, cache);
    problem.reset();

    const auto results = snapshot.evaluate({CandidateVariables{{500.0}}}, {});
    REQUIRE(results.size() == 1);
    CHECK(results.front().status == EvaluationStatus::Success);
    CHECK(results.front().objectives.front().value == 3.0);
    const auto cached_results = cached.evaluate_batch({CandidateVariables{{500.0}}}, {});
    REQUIRE(cached_results.size() == 1);
    CHECK(cached_results.front().objectives.front().value == 3.0);
}
