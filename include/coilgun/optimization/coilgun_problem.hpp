#pragma once

#include "coilgun/components/armature.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/optimization/constraint.hpp"
#include "coilgun/optimization/evaluator.hpp"
#include "coilgun/optimization/problem.hpp"
#include "coilgun/optimization/variables.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/multi_stage_sim.hpp"

#include <functional>
#include <cstdint>
#include <atomic>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

namespace coilgun::optimization {

class CudaBatchEvaluator;

enum class CoilgunParameter {
    CoilInnerRadius, CoilOuterRadius, CoilLength, CoilTurns, CoilPosition,
    ExcitationVoltage, ExcitationCapacitance, TriggerValue,
    ArmaturePosition, ArmatureVelocity, ArmatureMass
};

enum class CoilgunMetric {
    TerminalVelocity, MuzzleVelocity = TerminalVelocity,
    MaximumTemperature, PeakCurrent, PeakVoltage, Efficiency, EnergyLoss
};

struct CoilgunVariableBinding {
    std::string variable_id;
    CoilgunParameter parameter = CoilgunParameter::CoilPosition;
    std::size_t index = 0;
};

struct CoilgunExcitationConfig {
    double initial_voltage = 0.0;
    double capacitance = 0.0;
    bool crowbar = true;
};

/** Complete geometry/material data used when a coil geometry is variable. */
struct CoilgunCoilSpec {
    double inner_radius = 0.0;
    double outer_radius = 0.0;
    double length = 0.0;
    int turns = 0;
    double resistivity = 0.0;
    double wire_area = 0.0;
    double fill_factor = 0.0;
    double position = 0.0;
    bool force_exact_self_inductance = false;
    components::DrivingCoil make_coil() const;
};

struct CoilgunMetricConstraint {
    std::string id;
    CoilgunMetric metric = CoilgunMetric::TerminalVelocity;
    ConstraintDefinition definition;
};

class CoilgunOptimizationProblem final : public OptimizationProblem, public Evaluator, public BatchEvaluator {
public:
    struct Config {
        // `coils` is the convenient fixed-geometry path. `coil_specs`, when
        // supplied, retains all material parameters for geometry remapping.
        std::vector<components::DrivingCoil> coils;
        std::vector<CoilgunCoilSpec> coil_specs;
        std::optional<components::Armature> armature;
        std::vector<CoilgunExcitationConfig> excitations;
        std::vector<simulation::TriggerConfig> triggers;
        std::vector<CoilgunVariableBinding> bindings;
        std::vector<CoilgunMetricConstraint> constraints;
        std::string objective_id = "muzzle_velocity";
        double dt = 0.0;
        bool enable_thermal = false;
        simulation::OptimizationLevel optimization_level = simulation::OptimizationLevel::Full;
        simulation::TerminationPolicy termination = simulation::TerminationPolicy::defaults();
    };

    using GpuBatchEvaluator = std::function<std::vector<EvaluationResult>(
        const std::vector<CandidateVariables>&, const EvaluationContext&)>;

    CoilgunOptimizationProblem(VariableSchema schema, Config config);

    EvaluationResult evaluate(const CandidateVariables& variables) const override;
    EvaluationResult evaluate(const CandidateVariables& variables,
                              const EvaluationContext&) const override;
    std::vector<EvaluationResult> evaluate_batch(
        const std::vector<CandidateVariables>& variables,
        const EvaluationContext& context = {}) override;
    std::vector<EvaluationResult> evaluate_batch(
        const std::vector<CandidateVariables>& variables,
        const EvaluationContext& context = {}) const;
    std::optional<std::vector<EvaluationResult>> evaluate_batch_const(
        const std::vector<CandidateVariables>& variables,
        const EvaluationContext& context = {}) const override;

    void set_gpu_batch_evaluator(GpuBatchEvaluator evaluator);
    void clear_gpu_batch_evaluator();
    [[nodiscard]] bool has_gpu_batch_evaluator() const noexcept;
    [[nodiscard]] EvaluationCacheIdentity cache_identity() const override;
    [[nodiscard]] bool last_batch_used_fallback() const noexcept;
    [[nodiscard]] const VariableSchema& schema() const noexcept { return OptimizationProblem::spec()->schema(); }
    [[nodiscard]] const ProblemSpec& problem_spec() const noexcept { return *OptimizationProblem::spec(); }
    [[nodiscard]] const ProblemSpec* spec() const noexcept override { return OptimizationProblem::spec(); }
    [[nodiscard]] const Config& config() const noexcept { return config_; }
    [[nodiscard]] std::optional<EvaluationStatistics> statistics_snapshot() const override;

private:
    friend class CudaBatchEvaluator;
    EvaluationResult evaluate_cpu(const CandidateVariables&) const;
    EvaluationResult result_from_simulation(
        const simulation::MultiStageResult&, const components::Armature&,
        const std::vector<CoilgunExcitationConfig>&) const;
    std::vector<components::DrivingCoil> make_coils(const std::vector<double>&) const;
    double metric_value(CoilgunMetric, const simulation::MultiStageResult&, const components::Armature&,
                        const std::vector<CoilgunExcitationConfig>&) const;
    EvaluationResult invalid_result(const std::string&, const std::string&) const;
    [[nodiscard]] std::string canonical_fingerprint() const;
    [[nodiscard]] std::string canonical_fingerprint(bool has_gpu_evaluator,
                                                     std::uint64_t callback_identity) const;
    std::vector<EvaluationResult> evaluate_batch_with_callback(
        const std::vector<CandidateVariables>&, const EvaluationContext&,
        GpuBatchEvaluator) const;
    BatchEvaluationSnapshot make_evaluation_snapshot() override;

    Config config_;
    GpuBatchEvaluator gpu_evaluator_;
    std::uint64_t gpu_callback_identity_ = 0;
    mutable std::mutex gpu_callback_mutex_;
    mutable std::mutex statistics_mutex_;
    mutable bool last_batch_used_fallback_ = false;
    mutable EvaluationStatistics statistics_;
    static std::atomic<std::uint64_t> next_gpu_callback_identity_;
};

using CoilgunOptimizationConfig = CoilgunOptimizationProblem::Config;

} // namespace coilgun::optimization
