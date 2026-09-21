/**
 * @file test_cuda_batch_evaluator.cpp
 * @brief Production CUDA optimization evaluator tests.
 */

#include <doctest/doctest.h>

#include "coilgun/coilgun_cuda.hpp"
#include "coilgun/optimization/cuda_batch_evaluator.hpp"
#include "gpu_numerical_tolerances.hpp"

#include <cmath>
#include <barrier>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <vector>

namespace coilgun::optimization {
struct CudaBatchEvaluatorTestAccess {
    using Execution = std::function<CudaExecutionResponse(
        const std::vector<CandidateVariables>&, const EvaluationContext&)>;

    static std::shared_ptr<CudaBatchEvaluator> make(
        const CoilgunOptimizationProblem& problem, simulation::cuda::GpuBackend backend,
        CudaFallbackOptions options, Execution execution) {
        return std::shared_ptr<CudaBatchEvaluator>(
            new CudaBatchEvaluator(problem, std::move(backend), options, std::move(execution)));
    }
};
}

namespace {

using namespace coilgun;
using optimization::CandidateVariables;
using optimization::CoilgunOptimizationProblem;
using optimization::CoilgunParameter;
using optimization::CoilgunVariableBinding;
using optimization::VariableSchema;
using simulation::TriggerConfig;
using simulation::TriggerMode;
using simulation::cuda::BackendMode;
using simulation::cuda::GpuBackend;
using optimization::CudaBatchEvaluator;
using optimization::StatisticsBatchEvaluator;
using optimization::CudaExecutionResponse;
using optimization::CudaExecutionRow;
using optimization::CudaFallbackOptions;
using optimization::CudaFallbackPolicy;

static_assert(!std::is_constructible_v<
              CudaBatchEvaluator, const CoilgunOptimizationProblem&, GpuBackend,
              CudaFallbackOptions, optimization::CudaBatchEvaluatorTestAccess::Execution>,
              "arbitrary CUDA result injection must not be a public constructor");

components::Armature armature() {
    return components::Armature(0.005, 0.025, 0.080, physics::ALUMINUM.resistivity_ref,
                                physics::ALUMINUM.density, 0.0, 0.120, 2, 2, 0.000);
}

std::vector<components::DrivingCoil> coils() {
    return {
        components::DrivingCoil(0.010, 0.030, 0.050, 150, physics::COPPER.resistivity_ref,
                                1e-6, 0.7, 0.015),
        components::DrivingCoil(0.010, 0.030, 0.050, 150, physics::COPPER.resistivity_ref,
                                1e-6, 0.7, 0.085),
    };
}

CoilgunOptimizationProblem problem_with_bindings(
    CoilgunParameter parameter, std::size_t index = 0,
    std::vector<optimization::CoilgunMetricConstraint> constraints = {},
    simulation::OptimizationLevel optimization_level = simulation::OptimizationLevel::Full,
    double initial_voltage = 280.0) {
    VariableSchema schema({optimization::VariableSpec::continuous("value", 0.0, 1000.0)});
    CoilgunOptimizationProblem::Config config;
    config.coils = coils();
    config.armature = armature();
    config.excitations = {{initial_voltage, 0.0008, true}, {215.0, 0.0008, true}};
    config.triggers = {{TriggerMode::TimeDelay, 8e-6}};
    config.bindings = {{"value", parameter, index}};
    config.constraints = std::move(constraints);
    config.optimization_level = optimization_level;
    config.dt = 1e-6;
    config.termination.max_steps = 8;
    config.termination.enable_velocity_check = false;
    config.termination.enable_bound_check = false;
    return CoilgunOptimizationProblem(std::move(schema), std::move(config));
}

CoilgunOptimizationProblem problem_with_all_supported_bindings() {
    VariableSchema schema({
        optimization::VariableSpec::continuous("voltage", 0.0, 1000.0),
        optimization::VariableSpec::continuous("capacitance", 0.0, 0.01),
        optimization::VariableSpec::continuous("trigger", 0.0, 1e-3),
    });
    CoilgunOptimizationProblem::Config config;
    config.coils = coils();
    config.armature = armature();
    config.excitations = {{280.0, 0.0008, true}, {215.0, 0.0008, true}};
    config.triggers = {{TriggerMode::TimeDelay, 8e-6}};
    config.bindings = {
        {"voltage", CoilgunParameter::ExcitationVoltage, 0},
        {"capacitance", CoilgunParameter::ExcitationCapacitance, 0},
        {"trigger", CoilgunParameter::TriggerValue, 0},
    };
    config.dt = 1e-6;
    config.termination.max_steps = 8;
    config.termination.enable_velocity_check = false;
    config.termination.enable_bound_check = false;
    return CoilgunOptimizationProblem(std::move(schema), std::move(config));
}

GpuBackend direct_backend() {
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    backend.max_batch_sims = 16;
    return backend;
}

GpuBackend fallback_backend() {
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    backend.max_batch_sims = 16;
    return backend;
}

CudaExecutionResponse gpu_response(std::vector<optimization::EvaluationResult> results) {
    CudaExecutionResponse response;
    response.report.requested_backend = BackendMode::Direct;
    response.report.backend = BackendMode::Direct;
    response.report.gpu_executed = true;
    response.report.transfer_time_ms = 2.0;
    response.report.gpu_time_ms = 5.0;
    for (std::size_t index = 0; index < results.size(); ++index)
        response.rows.push_back(CudaExecutionRow{index, std::move(results[index])});
    return response;
}

std::function<CudaExecutionResponse(const std::vector<CandidateVariables>&,
                                    const optimization::EvaluationContext&)>
hook_for(CudaExecutionResponse response) {
    return [response = std::move(response)](const std::vector<CandidateVariables>&,
                                             const optimization::EvaluationContext&) mutable {
        return std::move(response);
    };
}

} // namespace

TEST_CASE("CUDA batch evaluator maps supported rows in stable order and proves GPU execution" * doctest::timeout(120)) {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaBatchEvaluator evaluator(problem, direct_backend());

    const std::vector<CandidateVariables> candidates{
        CandidateVariables{{280.0}}, CandidateVariables{{360.0}}, CandidateVariables{{440.0}},
    };
    const auto results = evaluator.evaluate_batch(candidates, {});

    REQUIRE(results.size() == candidates.size());
    for (const auto& result : results) {
        REQUIRE(result.status == optimization::EvaluationStatus::Success);
        REQUIRE(result.objectives.size() == 1);
        CHECK(result.objectives.front().id == problem.config().objective_id);
        CHECK(std::isfinite(result.objectives.front().value));
    }
    CHECK(std::stod(results[0].metadata.at("peak_voltage")) >= 280.0);
    CHECK(std::stod(results[1].metadata.at("peak_voltage")) >= 360.0);
    CHECK(std::stod(results[2].metadata.at("peak_voltage")) >= 440.0);

    const auto snapshot = evaluator.execution_snapshot();
    CHECK(snapshot.report.requested_backend == BackendMode::Direct);
    CHECK(snapshot.report.backend != BackendMode::Fallback);
    CHECK(snapshot.report.gpu_executed);
    CHECK(snapshot.report.gpu_time_ms > 0.0);
    CHECK(snapshot.host_time_ms >= snapshot.report.gpu_time_ms);
    REQUIRE(evaluator.statistics_snapshot());
    const auto stats = *evaluator.statistics_snapshot();
    CHECK(stats.gpu_requested_evaluations == 3);
    CHECK(stats.gpu_executed_evaluations == 3);
    CHECK(stats.gpu_successful_evaluations == 3);
    CHECK(stats.gpu_failed_evaluations == 0);
    CHECK(stats.cpu_fallback_evaluations == 0);
    CHECK(stats.gpu_batches == 1);
    CHECK(stats.gpu_failed_batches == 0);
    CHECK(stats.gpu_fallbacks == 0);
    CHECK(stats.gpu_transfer_seconds >= 0.0);
    CHECK(stats.gpu_kernel_seconds >= 0.0);
    CHECK(stats.gpu_elapsed_seconds >= 0.0);
}

TEST_CASE("CUDA batch evaluator expands locally invalid rows without reordering valid rows" * doctest::timeout(120)) {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaBatchEvaluator evaluator(problem, direct_backend());

    const std::vector<CandidateVariables> candidates{
        CandidateVariables{{280.0}},
        CandidateVariables{{std::numeric_limits<double>::quiet_NaN()}},
        CandidateVariables{{440.0}},
    };
    const auto results = evaluator.evaluate_batch(candidates, {});

    REQUIRE(results.size() == candidates.size());
    CHECK(results[0].status == optimization::EvaluationStatus::Success);
    CHECK(results[1].status == optimization::EvaluationStatus::Invalid);
    REQUIRE(results[1].diagnostics.size() == 1);
    CHECK(results[1].diagnostics.front().code == "non_finite_variable");
    CHECK(results[2].status == optimization::EvaluationStatus::Success);
    CHECK(std::stod(results[0].metadata.at("peak_voltage")) >= 280.0);
    CHECK(std::stod(results[2].metadata.at("peak_voltage")) >= 440.0);
}

TEST_CASE("CUDA batch evaluator records no GPU timing for invalid-only batches") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaBatchEvaluator evaluator(problem, direct_backend());
    const auto results = evaluator.evaluate_batch({
        CandidateVariables{{std::numeric_limits<double>::quiet_NaN()}},
        CandidateVariables{{}},
    });

    REQUIRE(results.size() == 2);
    CHECK(results[0].status == optimization::EvaluationStatus::Invalid);
    CHECK(results[1].status == optimization::EvaluationStatus::Invalid);
    REQUIRE(evaluator.statistics_snapshot());
    const auto stats = *evaluator.statistics_snapshot();
    CHECK(stats.gpu_requested_evaluations == 0);
    CHECK(stats.gpu_batches == 0);
    CHECK(stats.gpu_failed_batches == 0);
    CHECK(stats.gpu_executed_evaluations == 0);
    CHECK(stats.gpu_transfer_seconds == 0.0);
    CHECK(stats.gpu_kernel_seconds == 0.0);
    CHECK(stats.gpu_elapsed_seconds == 0.0);
    CHECK(evaluator.execution_snapshot().host_time_ms == 0.0);
}

TEST_CASE("CUDA batch evaluator rejects zero and negative voltage rows locally" * doctest::timeout(120)) {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaBatchEvaluator evaluator(problem, direct_backend());
    const auto results = evaluator.evaluate_batch({
        CandidateVariables{{280.0}}, CandidateVariables{{0.0}},
        CandidateVariables{{-1.0}}, CandidateVariables{{440.0}},
    });

    REQUIRE(results.size() == 4);
    CHECK(results[0].status == optimization::EvaluationStatus::Success);
    CHECK(results[1].status == optimization::EvaluationStatus::Invalid);
    CHECK(results[2].status == optimization::EvaluationStatus::Invalid);
    CHECK(results[3].status == optimization::EvaluationStatus::Success);
    CHECK(std::stod(results[0].metadata.at("peak_voltage")) >= 280.0);
    CHECK(std::stod(results[3].metadata.at("peak_voltage")) >= 440.0);
}

TEST_CASE("CUDA batch evaluator maps capacitance and trigger bindings" * doctest::timeout(120)) {
    auto problem = problem_with_all_supported_bindings();
    CudaBatchEvaluator evaluator(problem, direct_backend());
    const auto results = evaluator.evaluate_batch({
        CandidateVariables{{280.0, 0.0008, 8e-6}},
        CandidateVariables{{280.0, 0.0004, 8e-6}},
        CandidateVariables{{280.0, 0.0008, 0.0}},
    });
    REQUIRE(results.size() == 3);
    for (const auto& result : results)
        CHECK(result.status == optimization::EvaluationStatus::Success);
    CHECK(results[0].metadata.at("terminal_velocity") != results[1].metadata.at("terminal_velocity"));
    CHECK(results[0].metadata.at("terminal_velocity") != results[2].metadata.at("terminal_velocity"));
}

TEST_CASE("CUDA batch evaluator rejects unsupported bindings and configurations") {
    CHECK_THROWS_WITH_AS(
        CudaBatchEvaluator(problem_with_bindings(CoilgunParameter::ArmaturePosition, 0), direct_backend()),
        "CUDA batch evaluator does not support candidate binding: ArmaturePosition", std::invalid_argument);

    auto peak_current_problem = problem_with_bindings(
        CoilgunParameter::ExcitationVoltage, 0,
        {{"peak_current", optimization::CoilgunMetric::PeakCurrent,
          optimization::ConstraintDefinition{"peak_current", optimization::ConstraintKind::Hard,
                                              optimization::ConstraintRelation::LessEqual, 1.0, 1.0, 1.0, 0}}});
    CHECK_THROWS_WITH_AS(CudaBatchEvaluator(peak_current_problem, direct_backend()),
                         "CUDA batch evaluator does not support PeakCurrent metrics", std::invalid_argument);
    CHECK_THROWS_WITH_AS(
        CudaBatchEvaluator(problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0, {},
                                                 simulation::OptimizationLevel::Full, 0.0), direct_backend()),
        "CUDA batch evaluator excitation initial_voltage must be finite and positive; capacitance must be finite and positive",
        std::invalid_argument);
    CHECK_THROWS_WITH_AS(
        CudaBatchEvaluator(problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0, {},
                                                 simulation::OptimizationLevel::Full, -1.0), direct_backend()),
        "CUDA batch evaluator excitation initial_voltage must be finite and positive; capacitance must be finite and positive",
        std::invalid_argument);
}

TEST_CASE("CUDA batch evaluator keeps invalid capacitance rows local" * doctest::timeout(120)) {
    auto problem = problem_with_all_supported_bindings();
    CudaBatchEvaluator evaluator(problem, direct_backend());
    const auto results = evaluator.evaluate_batch({
        CandidateVariables{{280.0, 0.0008, 8e-6}},
        CandidateVariables{{280.0, 0.0, 8e-6}},
        CandidateVariables{{360.0, 0.0008, 8e-6}},
    });

    REQUIRE(results.size() == 3);
    CHECK(results[0].status == optimization::EvaluationStatus::Success);
    CHECK(results[1].status == optimization::EvaluationStatus::Invalid);
    CHECK(results[1].diagnostics.front().code == "invalid_configuration");
    CHECK(results[2].status == optimization::EvaluationStatus::Success);
}

TEST_CASE("CUDA batch evaluator exposes backend fallback instead of CPU results" * doctest::timeout(120)) {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaBatchEvaluator evaluator(problem, fallback_backend());
    const auto results = evaluator.evaluate_batch({CandidateVariables{{280.0}}});

    REQUIRE(results.size() == 1);
    CHECK(results.front().status == optimization::EvaluationStatus::Failed);
    REQUIRE(results.front().diagnostics.size() == 1);
    CHECK(results.front().diagnostics.front().code == "gpu_backend_fallback");
    const auto snapshot = evaluator.execution_snapshot();
    CHECK(snapshot.report.backend == BackendMode::Fallback);
    CHECK_FALSE(snapshot.report.gpu_executed);
    REQUIRE(evaluator.statistics_snapshot());
    CHECK(evaluator.statistics_snapshot()->fallbacks == 1);
}

TEST_CASE("CUDA batch evaluator exposes fallback policy without exposing execution injection") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaFallbackOptions options;
    options.policy = CudaFallbackPolicy::WholeBatchCpu;
    CudaBatchEvaluator evaluator(problem, direct_backend(), options);
    CHECK(evaluator.fallback_options().policy == CudaFallbackPolicy::WholeBatchCpu);
}

TEST_CASE("CUDA batch evaluator composes with statistics wrapper without double counting" * doctest::timeout(120)) {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    auto evaluator = std::make_shared<CudaBatchEvaluator>(problem, direct_backend());
    StatisticsBatchEvaluator wrapped(evaluator);
    auto collector = std::make_shared<optimization::EvaluationStatisticsCollector>();
    const auto results = wrapped.evaluate_batch({CandidateVariables{{280.0}}},
                                                 optimization::EvaluationContext{19, false, collector});

    REQUIRE(results.size() == 1);
    CHECK(results.front().status == optimization::EvaluationStatus::Success);
    CHECK(collector->snapshot().evaluations == 1);
    CHECK(collector->snapshot().successful_evaluations == 1);
    CHECK(collector->snapshot().failed_evaluations == 0);
    CHECK(collector->snapshot().fallbacks == 0);
    REQUIRE(evaluator->statistics_snapshot());
    CHECK(evaluator->statistics_snapshot()->evaluations == 0);
    CHECK(evaluator->statistics_snapshot()->successful_evaluations == 0);
    CHECK(evaluator->statistics_snapshot()->failed_evaluations == 0);
}

TEST_CASE("shared CUDA evaluator snapshots retain the evaluator owner") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    auto result = optimization::EvaluationResult::success();
    result.objectives.push_back({"muzzle_velocity", 8.0, true});
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), {}, hook_for(gpu_response({std::move(result)})));
    auto snapshot = evaluator->evaluation_snapshot();
    evaluator.reset();

    const auto results = snapshot.evaluate({CandidateVariables{{280.0}}}, {});
    REQUIRE(results.size() == 1);
    CHECK(results.front().status == optimization::EvaluationStatus::Success);
    CHECK(results.front().objectives.front().value == 8.0);
}

TEST_CASE("CUDA batch evaluator isolates synthetic row failure with per-candidate fallback") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    auto response = gpu_response({optimization::EvaluationResult::success(),
                                  optimization::EvaluationResult::failed("gpu_row", "device row failed"),
                                  optimization::EvaluationResult::success()});
    CudaFallbackOptions options;
    options.policy = CudaFallbackPolicy::PerCandidateCpu;
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), options, hook_for(std::move(response)));
    const auto results = evaluator->evaluate_batch({CandidateVariables{{280.0}},
                                                     CandidateVariables{{360.0}},
                                                     CandidateVariables{{440.0}}});
    REQUIRE(results.size() == 3);
    CHECK(results[0].status == optimization::EvaluationStatus::Success);
    CHECK(results[1].status == optimization::EvaluationStatus::Success);
    CHECK(results[1].metadata.at("gpu_failure_code") == "gpu_row");
    CHECK(results[1].metadata.at("cpu_fallback") == "true");
    CHECK(results[2].status == optimization::EvaluationStatus::Success);
    REQUIRE(evaluator->statistics_snapshot());
    const auto stats = *evaluator->statistics_snapshot();
    CHECK(stats.gpu_requested_evaluations == 3);
    CHECK(stats.gpu_executed_evaluations == 3);
    CHECK(stats.gpu_successful_evaluations == 2);
    CHECK(stats.gpu_failed_evaluations == 1);
    CHECK(stats.cpu_fallback_evaluations == 1);
    CHECK(stats.gpu_batches == 1);
    CHECK(stats.gpu_failed_batches == 0);
    CHECK(stats.gpu_fallbacks == 0);
    CHECK(stats.gpu_transfer_seconds == doctest::Approx(0.002));
    CHECK(stats.gpu_kernel_seconds == doctest::Approx(0.005));
    CHECK(stats.gpu_elapsed_seconds >= 0.0);
}

TEST_CASE("CUDA batch evaluator strict row failure remains failed") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    auto response = gpu_response({optimization::EvaluationResult::success(),
                                  optimization::EvaluationResult::failed("gpu_row", "device row failed")});
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), {}, hook_for(std::move(response)));
    const auto results = evaluator->evaluate_batch({CandidateVariables{{280.0}}, CandidateVariables{{360.0}}});
    REQUIRE(results.size() == 2);
    CHECK(results[0].status == optimization::EvaluationStatus::Success);
    CHECK(results[1].status == optimization::EvaluationStatus::Failed);
    REQUIRE(evaluator->statistics_snapshot());
    CHECK(evaluator->statistics_snapshot()->gpu_failed_evaluations == 1);
    CHECK(evaluator->statistics_snapshot()->cpu_fallback_evaluations == 0);
}

TEST_CASE("CUDA batch evaluator strict policy converts every non-success GPU row to Failed") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    auto response = gpu_response({optimization::EvaluationResult::invalid("gpu_invalid", "device rejected row")});
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), {}, hook_for(std::move(response)));
    const auto results = evaluator->evaluate_batch({CandidateVariables{{280.0}}});
    REQUIRE(results.size() == 1);
    CHECK(results.front().status == optimization::EvaluationStatus::Failed);
    REQUIRE(results.front().diagnostics.size() == 1);
    CHECK(results.front().diagnostics.front().code == "gpu_invalid");
}

TEST_CASE("CUDA batch evaluator whole-batch fallback counts one event and preserves invalid rows") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaExecutionResponse response;
    response.report.requested_backend = BackendMode::Direct;
    response.report.backend = BackendMode::Fallback;
    response.report.fallback_reason = "synthetic unavailable device";
    CudaFallbackOptions options;
    options.policy = CudaFallbackPolicy::WholeBatchCpu;
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, fallback_backend(), options, hook_for(std::move(response)));
    const auto results = evaluator->evaluate_batch({CandidateVariables{{280.0}},
                                                     CandidateVariables{{std::numeric_limits<double>::quiet_NaN()}},
                                                     CandidateVariables{{440.0}}});
    REQUIRE(results.size() == 3);
    CHECK(results[0].status == optimization::EvaluationStatus::Success);
    CHECK(results[1].status == optimization::EvaluationStatus::Invalid);
    CHECK(results[2].status == optimization::EvaluationStatus::Success);
    CHECK(results[0].metadata.at("cpu_fallback") == "true");
    CHECK(results[2].metadata.at("cpu_fallback") == "true");
    REQUIRE(evaluator->statistics_snapshot());
    const auto stats = *evaluator->statistics_snapshot();
    CHECK(stats.gpu_requested_evaluations == 2);
    CHECK(stats.gpu_executed_evaluations == 0);
    CHECK(stats.gpu_failed_batches == 1);
    CHECK(stats.gpu_fallbacks == 1);
    CHECK(stats.cpu_fallback_evaluations == 2);
}

TEST_CASE("CUDA batch evaluator reports protocol mismatch without CPU repair") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaExecutionResponse response = gpu_response({optimization::EvaluationResult::success()});
    CudaFallbackOptions options;
    options.policy = CudaFallbackPolicy::WholeBatchCpu;
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), options, hook_for(std::move(response)));
    const auto results = evaluator->evaluate_batch({CandidateVariables{{280.0}}, CandidateVariables{{360.0}}});
    REQUIRE(results.size() == 2);
    for (const auto& result : results) {
        CHECK(result.status == optimization::EvaluationStatus::Failed);
        REQUIRE(result.diagnostics.size() == 1);
        CHECK(result.diagnostics.front().code == "gpu_protocol_error");
    }
    REQUIRE(evaluator->statistics_snapshot());
    CHECK(evaluator->statistics_snapshot()->gpu_failed_batches == 1);
    CHECK(evaluator->statistics_snapshot()->cpu_fallback_evaluations == 0);
}

TEST_CASE("CUDA batch evaluator rejects out-of-order rows without CPU repair") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaExecutionResponse response = gpu_response({optimization::EvaluationResult::success(),
                                                   optimization::EvaluationResult::success()});
    std::swap(response.rows[0].index, response.rows[1].index);
    CudaFallbackOptions options;
    options.policy = CudaFallbackPolicy::WholeBatchCpu;
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), options, hook_for(std::move(response)));
    const auto results = evaluator->evaluate_batch({CandidateVariables{{280.0}}, CandidateVariables{{360.0}}});
    REQUIRE(results.size() == 2);
    for (const auto& result : results) {
        CHECK(result.status == optimization::EvaluationStatus::Failed);
        REQUIRE(result.diagnostics.size() == 1);
        CHECK(result.diagnostics.front().code == "gpu_protocol_error");
    }
    REQUIRE(evaluator->statistics_snapshot());
    CHECK(evaluator->statistics_snapshot()->gpu_failed_batches == 1);
    CHECK(evaluator->statistics_snapshot()->cpu_fallback_evaluations == 0);
}

TEST_CASE("CUDA batch evaluator treats thrown batches as one strict fallback event") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), {},
        [](const std::vector<CandidateVariables>&, const optimization::EvaluationContext&) -> CudaExecutionResponse {
            throw std::runtime_error("synthetic CUDA launch failure");
        });

    const auto results = evaluator->evaluate_batch({CandidateVariables{{280.0}}, CandidateVariables{{360.0}}});
    REQUIRE(results.size() == 2);
    for (const auto& result : results) {
        CHECK(result.status == optimization::EvaluationStatus::Failed);
        REQUIRE(result.diagnostics.size() == 1);
        CHECK(result.diagnostics.front().code == "gpu_execution_failed");
    }
    REQUIRE(evaluator->statistics_snapshot());
    const auto stats = *evaluator->statistics_snapshot();
    CHECK(stats.gpu_batches == 1);
    CHECK(stats.gpu_failed_batches == 1);
    CHECK(stats.gpu_fallbacks == 1);
    CHECK(stats.cpu_fallback_evaluations == 0);
    CHECK(stats.fallbacks == 1);
}

TEST_CASE("CUDA batch evaluator whole-batch throw fallback evaluates eligible rows exactly once") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    CudaFallbackOptions options;
    options.policy = CudaFallbackPolicy::WholeBatchCpu;
    std::size_t hook_calls = 0;
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), options,
        [&hook_calls](const std::vector<CandidateVariables>&, const optimization::EvaluationContext&) -> CudaExecutionResponse {
            ++hook_calls;
            throw std::runtime_error("synthetic CUDA launch failure");
        });

    const auto results = evaluator->evaluate_batch({CandidateVariables{{280.0}},
                                                     CandidateVariables{{std::numeric_limits<double>::quiet_NaN()}},
                                                     CandidateVariables{{440.0}}});
    REQUIRE(results.size() == 3);
    CHECK(hook_calls == 1);
    CHECK(results[0].status == optimization::EvaluationStatus::Success);
    CHECK(results[1].status == optimization::EvaluationStatus::Invalid);
    CHECK(results[2].status == optimization::EvaluationStatus::Success);
    CHECK(results[0].metadata.at("cpu_fallback") == "true");
    CHECK(results[2].metadata.at("cpu_fallback") == "true");
    REQUIRE(evaluator->statistics_snapshot());
    const auto stats = *evaluator->statistics_snapshot();
    CHECK(stats.gpu_requested_evaluations == 2);
    CHECK(stats.gpu_executed_evaluations == 0);
    CHECK(stats.gpu_failed_batches == 1);
    CHECK(stats.gpu_fallbacks == 1);
    CHECK(stats.cpu_fallback_evaluations == 2);
}

TEST_CASE("CUDA metrics compose through statistics and cache wrappers without duplication") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    auto execution = [](const std::vector<CandidateVariables>& candidates,
                        const optimization::EvaluationContext&) {
        CudaExecutionResponse response;
        response.report.requested_backend = BackendMode::Direct;
        response.report.backend = BackendMode::Direct;
        response.report.gpu_executed = true;
        response.report.transfer_time_ms = 4.0;
        response.report.gpu_time_ms = 6.0;
        for (std::size_t i = 0; i < candidates.size(); ++i)
            response.rows.push_back({i, optimization::EvaluationResult::success()});
        return response;
    };
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), CudaFallbackOptions{}, execution);
    auto tracked = std::make_shared<StatisticsBatchEvaluator>(evaluator);
    auto cached = std::make_shared<optimization::CachedBatchEvaluator>(
        tracked, std::make_shared<optimization::InMemoryEvaluationCache>());
    auto collector = std::make_shared<optimization::EvaluationStatisticsCollector>();
    const optimization::EvaluationContext context{31, false, collector};

    const auto first = cached->evaluate_batch({CandidateVariables{{280.0}}}, context);
    const auto second = cached->evaluate_batch({CandidateVariables{{280.0}}}, context);
    REQUIRE(first.size() == 1);
    REQUIRE(second.size() == 1);
    CHECK(first.front().status == optimization::EvaluationStatus::Success);
    CHECK(second.front().status == optimization::EvaluationStatus::Success);

    const auto run = collector->snapshot();
    CHECK(run.cache_hits == 1);
    CHECK(run.gpu_requested_evaluations == 1);
    CHECK(run.gpu_executed_evaluations == 1);
    CHECK(run.gpu_successful_evaluations == 1);
    CHECK(run.gpu_failed_evaluations == 0);
    CHECK(run.cpu_fallback_evaluations == 0);
    CHECK(run.gpu_batches == 1);
    CHECK(run.gpu_failed_batches == 0);
    CHECK(run.gpu_fallbacks == 0);
    CHECK(run.gpu_transfer_seconds == doctest::Approx(0.004));
    CHECK(run.gpu_kernel_seconds == doctest::Approx(0.006));
    REQUIRE(cached->statistics_snapshot());
    const auto aggregate = *cached->statistics_snapshot();
    CHECK(aggregate.gpu_requested_evaluations == 1);
    CHECK(aggregate.gpu_batches == 1);
    CHECK(aggregate.gpu_transfer_seconds == doctest::Approx(0.004));
}

TEST_CASE("concurrent CUDA evaluator calls isolate run-local GPU metrics") {
    auto problem = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    std::barrier calls(2);
    auto execution = [&calls](const std::vector<CandidateVariables>& candidates,
                              const optimization::EvaluationContext&) {
        calls.arrive_and_wait();
        CudaExecutionResponse response;
        response.report.requested_backend = BackendMode::Direct;
        response.report.backend = BackendMode::Direct;
        response.report.gpu_executed = true;
        response.report.transfer_time_ms = 1.0;
        response.report.gpu_time_ms = 2.0;
        for (std::size_t i = 0; i < candidates.size(); ++i)
            response.rows.push_back({i, optimization::EvaluationResult::success()});
        return response;
    };
    auto evaluator = optimization::CudaBatchEvaluatorTestAccess::make(
        problem, direct_backend(), CudaFallbackOptions{}, execution);
    auto first_collector = std::make_shared<optimization::EvaluationStatisticsCollector>();
    auto second_collector = std::make_shared<optimization::EvaluationStatisticsCollector>();
    std::thread first([&] {
        evaluator->evaluate_batch({CandidateVariables{{280.0}}},
                                   {41, false, first_collector});
    });
    std::thread second([&] {
        evaluator->evaluate_batch({CandidateVariables{{360.0}}},
                                   {42, false, second_collector});
    });
    first.join();
    second.join();

    for (const auto& collector : {first_collector, second_collector}) {
        const auto stats = collector->snapshot();
        CHECK(stats.gpu_requested_evaluations == 1);
        CHECK(stats.gpu_executed_evaluations == 1);
        CHECK(stats.gpu_successful_evaluations == 1);
        CHECK(stats.gpu_batches == 1);
        CHECK(stats.gpu_transfer_seconds == doctest::Approx(0.001));
        CHECK(stats.gpu_kernel_seconds == doctest::Approx(0.002));
        CHECK(std::isfinite(stats.gpu_elapsed_seconds));
        CHECK(stats.gpu_elapsed_seconds >= 0.0);
    }
    REQUIRE(evaluator->statistics_snapshot());
    CHECK(evaluator->statistics_snapshot()->gpu_batches == 2);
    CHECK(evaluator->statistics_snapshot()->gpu_requested_evaluations == 2);
}

TEST_CASE("CUDA batch evaluator matches CPU Reference terminal velocity and constraints" * doctest::timeout(120)) {
    const optimization::CoilgunMetricConstraint peak_voltage{
        "peak_voltage_limit", optimization::CoilgunMetric::PeakVoltage,
        optimization::ConstraintDefinition{"peak_voltage_limit", optimization::ConstraintKind::Hard,
                                            optimization::ConstraintRelation::LessEqual, 500.0, 500.0, 500.0, 0}};
    auto reference_problem = problem_with_bindings(
        CoilgunParameter::ExcitationVoltage, 0, {peak_voltage}, simulation::OptimizationLevel::Reference);
    auto full_problem = problem_with_bindings(
        CoilgunParameter::ExcitationVoltage, 0, {peak_voltage}, simulation::OptimizationLevel::Full);
    CudaBatchEvaluator evaluator(full_problem, direct_backend());
    const CandidateVariables candidate{{360.0}};
    const auto cpu = reference_problem.evaluate(candidate);
    const auto gpu = evaluator.evaluate_batch({candidate});

    REQUIRE(cpu.status == optimization::EvaluationStatus::Success);
    REQUIRE(gpu.size() == 1);
    REQUIRE(gpu.front().status == optimization::EvaluationStatus::Success);
    REQUIRE(cpu.objectives.size() == 1);
    REQUIRE(gpu.front().objectives.size() == 1);
    CHECK(gpu_test::numerically_equal(gpu.front().objectives.front().value,
                                      cpu.objectives.front().value,
                                      gpu_test::tolerance_for(simulation::cuda::GpuOptLevel::Full)));
    REQUIRE(cpu.constraints.size() == 1);
    REQUIRE(gpu.front().constraints.size() == 1);
    CHECK(gpu.front().constraints.front().id == cpu.constraints.front().id);
    CHECK(gpu_test::numerically_equal(gpu.front().constraints.front().value,
                                      cpu.constraints.front().value,
                                      gpu_test::tolerance_for(simulation::cuda::GpuOptLevel::Full)));
    CHECK(gpu.front().metadata.count("cuda_gpu_executed") == 1);
}

TEST_CASE("CUDA evaluator owns its problem for snapshots and CPU fallback") {
    auto source = problem_with_bindings(CoilgunParameter::ExcitationVoltage, 0);
    auto problem = std::make_unique<CoilgunOptimizationProblem>(
        source.schema(), source.config());
    auto evaluator = std::make_shared<CudaBatchEvaluator>(
        *problem, fallback_backend(), CudaFallbackOptions{CudaFallbackPolicy::WholeBatchCpu});
    auto snapshot = evaluator->evaluation_snapshot();
    problem.reset();
    evaluator.reset();

    const auto results = snapshot.evaluate({CandidateVariables{{280.0}}}, {});
    REQUIRE(results.size() == 1);
    CHECK(results.front().status == optimization::EvaluationStatus::Success);
    CHECK(results.front().metadata.at("cpu_fallback") == "true");
}

TEST_CASE("CUDA production cache identity isolates problem backend and fallback configuration" * doctest::timeout(120)) {
    auto first_problem_object = problem_with_bindings(
        CoilgunParameter::ExcitationVoltage, 0, {}, simulation::OptimizationLevel::Full, 280.0);
    auto second_problem_object = problem_with_bindings(
        CoilgunParameter::ExcitationVoltage, 0, {}, simulation::OptimizationLevel::Full, 281.0);
    auto first_problem = std::make_shared<CoilgunOptimizationProblem>(
        first_problem_object.schema(), first_problem_object.config());
    auto second_problem = std::make_shared<CoilgunOptimizationProblem>(
        second_problem_object.schema(), second_problem_object.config());
    auto first_evaluator = std::make_shared<CudaBatchEvaluator>(*first_problem, direct_backend());
    auto second_evaluator = std::make_shared<CudaBatchEvaluator>(*second_problem, direct_backend());
    CHECK(first_evaluator->cache_identity().namespace_id != second_evaluator->cache_identity().namespace_id);

    auto cache = std::make_shared<optimization::InMemoryEvaluationCache>();
    optimization::CachedBatchEvaluator first(first_evaluator, cache);
    optimization::CachedBatchEvaluator second(second_evaluator, cache);
    const CandidateVariables candidate{{360.0}};
    REQUIRE(first.evaluate_batch({candidate}, {}).front().status == optimization::EvaluationStatus::Success);
    REQUIRE(second.evaluate_batch({candidate}, {}).front().status == optimization::EvaluationStatus::Success);
    CHECK(first.statistics().cache_hits == 0);
    CHECK(second.statistics().cache_hits == 0);
    REQUIRE(first.evaluate_batch({candidate}, {}).front().status == optimization::EvaluationStatus::Success);
    CHECK(first.statistics().cache_hits == 1);

    auto alternate_backend = direct_backend();
    alternate_backend.threads_per_block = 256;
    CudaBatchEvaluator alternate(*first_problem, alternate_backend);
    CHECK(alternate.cache_identity().namespace_id != first_evaluator->cache_identity().namespace_id);
    CudaBatchEvaluator alternate_policy(*first_problem, direct_backend(),
                                        CudaFallbackOptions{CudaFallbackPolicy::PerCandidateCpu});
    CHECK(alternate_policy.cache_identity().namespace_id != first_evaluator->cache_identity().namespace_id);
}
