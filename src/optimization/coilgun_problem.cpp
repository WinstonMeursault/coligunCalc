#include "coilgun/optimization/coilgun_problem.hpp"

#include "coilgun/physics/constants.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

namespace coilgun::optimization {
namespace {

bool finite_vector(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); });
}

bool valid_gpu_success(const EvaluationResult& result,
                       const CoilgunOptimizationProblem::Config& config) {
    if (result.objectives.size() != 1) return false;
    const auto& objective = result.objectives.front();
    if (objective.id != config.objective_id || !objective.maximize ||
        !std::isfinite(objective.value))
        return false;
    if (result.constraints.size() != config.constraints.size()) return false;
    for (std::size_t index = 0; index < config.constraints.size(); ++index) {
        const auto& expected = config.constraints[index];
        const auto& actual = result.constraints[index];
        if (actual.id != expected.id || actual.kind != expected.definition.kind ||
            actual.relation != expected.definition.relation ||
            actual.lower_bound != expected.definition.lower_bound ||
            actual.upper_bound != expected.definition.upper_bound ||
            actual.priority != expected.definition.priority ||
            !std::isfinite(actual.value) || !std::isfinite(actual.lower_bound) ||
            !std::isfinite(actual.upper_bound) || !std::isfinite(actual.violation) ||
            !std::isfinite(actual.normalized_violation))
            return false;
        try {
            const auto derived = expected.definition.evaluate(actual.value);
            if (actual.violation != derived.violation ||
                actual.normalized_violation != derived.normalized_violation ||
                actual.satisfied != derived.satisfied)
                return false;
        } catch (...) {
            return false;
        }
    }
    return true;
}

std::string number(double value) {
    std::ostringstream stream;
    stream << std::setprecision(17) << value;
    return stream.str();
}

class CanonicalWriter {
public:
    void string(std::string_view value) {
        output_ += std::to_string(value.size());
        output_.push_back(':');
        output_.append(value);
        output_.push_back(';');
    }

    void boolean(bool value) { string(value ? "1" : "0"); }
    void integer(std::int64_t value) { string(std::to_string(value)); }
    void size(std::size_t value) { string(std::to_string(value)); }
    void real(double value) {
        std::uint64_t bits = 0;
        static_assert(sizeof(bits) == sizeof(value));
        std::memcpy(&bits, &value, sizeof(bits));
        string(std::to_string(bits));
    }

    [[nodiscard]] std::string finish() && { return std::move(output_); }

private:
    std::string output_;
};

template <typename Enum>
void enum_value(CanonicalWriter& writer, Enum value) {
    writer.integer(static_cast<std::underlying_type_t<Enum>>(value));
}

void write_coil(CanonicalWriter& writer, const components::DrivingCoil& coil) {
    writer.real(coil.inner_radius());
    writer.real(coil.outer_radius());
    writer.real(coil.length());
    writer.integer(coil.turns());
    writer.real(coil.resistivity());
    writer.real(coil.wire_area());
    writer.real(coil.fill_factor());
    writer.real(coil.position());
    writer.boolean(coil.force_exact_self_inductance());
    writer.real(coil.resistance());
    writer.real(coil.self_inductance());
}

void write_coil_spec(CanonicalWriter& writer, const CoilgunCoilSpec& coil) {
    writer.real(coil.inner_radius);
    writer.real(coil.outer_radius);
    writer.real(coil.length);
    writer.integer(coil.turns);
    writer.real(coil.resistivity);
    writer.real(coil.wire_area);
    writer.real(coil.fill_factor);
    writer.real(coil.position);
    writer.boolean(coil.force_exact_self_inductance);
}

} // namespace

std::atomic<std::uint64_t> CoilgunOptimizationProblem::next_gpu_callback_identity_{0};

components::DrivingCoil CoilgunCoilSpec::make_coil() const {
    return components::DrivingCoil(inner_radius, outer_radius, length, turns, resistivity,
                                    wire_area, fill_factor, position, force_exact_self_inductance);
}

CoilgunOptimizationProblem::CoilgunOptimizationProblem(VariableSchema schema, Config config)
    : config_(std::move(config)) {
    if (!config_.armature) throw std::invalid_argument("coilgun optimization requires an armature");
    if (config_.coils.empty() && config_.coil_specs.empty())
        throw std::invalid_argument("coilgun optimization requires at least one coil");
    if (!config_.coil_specs.empty() && config_.coil_specs.size() != config_.coils.size() &&
        !config_.coils.empty())
        throw std::invalid_argument("coils and coil_specs must have equal sizes");
    const auto stage_count = config_.coil_specs.empty() ? config_.coils.size() : config_.coil_specs.size();
    if (config_.excitations.size() != stage_count)
        throw std::invalid_argument("excitation count must equal coil count");
    if (config_.triggers.size() + 1 != stage_count)
        throw std::invalid_argument("trigger count must equal coil count minus one");
    if (!std::isfinite(config_.dt) || config_.dt <= 0.0)
        throw std::invalid_argument("dt must be finite and positive");
    if (config_.objective_id.empty()) throw std::invalid_argument("objective ID must not be empty");
    for (const auto& binding : config_.bindings) {
        bool found = false;
        for (const auto& variable : schema.variables()) found |= variable.id == binding.variable_id;
        if (!found) throw std::invalid_argument("binding references unknown variable: " + binding.variable_id);
        switch (binding.parameter) {
        case CoilgunParameter::CoilInnerRadius:
        case CoilgunParameter::CoilOuterRadius:
        case CoilgunParameter::CoilLength:
        case CoilgunParameter::CoilTurns:
        case CoilgunParameter::CoilPosition:
            if (binding.index >= stage_count) throw std::out_of_range("coil binding index out of range");
            if (config_.coil_specs.empty())
                throw std::invalid_argument("coil geometry bindings require coil_specs");
            break;
        case CoilgunParameter::ExcitationVoltage:
        case CoilgunParameter::ExcitationCapacitance:
            if (binding.index >= config_.excitations.size()) throw std::out_of_range("excitation binding index out of range");
            break;
        case CoilgunParameter::TriggerValue:
            if (binding.index >= config_.triggers.size()) throw std::out_of_range("trigger binding index out of range");
            break;
        case CoilgunParameter::ArmatureMass:
            throw std::invalid_argument("ArmatureMass bindings are not supported");
        default: break;
        }
    }
    for (auto& constraint : config_.constraints) {
        if (constraint.id.empty()) constraint.id = constraint.definition.id;
        if (constraint.definition.id.empty()) constraint.definition.id = constraint.id;
        if (constraint.id.empty()) throw std::invalid_argument("constraint ID must not be empty");
        // The adapter's public metric ID is the sole identifier exposed in
        // reports and therefore the one published by its ProblemSpec.
        constraint.definition.id = constraint.id;
        if (!config_.enable_thermal && constraint.metric == CoilgunMetric::MaximumTemperature)
            throw std::invalid_argument("maximum-temperature constraints require thermal simulation to be enabled");
        constraint.definition.validate();
    }
    std::vector<ConstraintDefinition> constraint_definitions;
    constraint_definitions.reserve(config_.constraints.size());
    for (const auto& constraint : config_.constraints)
        constraint_definitions.push_back(constraint.definition);
    set_spec(ProblemSpec(std::move(schema),
                         {ObjectiveDefinition{config_.objective_id, true, 1.0}},
                         std::move(constraint_definitions)));
}

EvaluationResult CoilgunOptimizationProblem::invalid_result(const std::string& code,
                                                             const std::string& message) const {
    return EvaluationResult::invalid(code, message);
}

EvaluationResult CoilgunOptimizationProblem::result_from_simulation(
    const simulation::MultiStageResult& result, const components::Armature& armature,
    const std::vector<CoilgunExcitationConfig>& excitations) const {
    const double velocity = metric_value(CoilgunMetric::TerminalVelocity, result, armature, excitations);
    if (!std::isfinite(velocity))
        return EvaluationResult::failed("non_finite_result", "simulation returned a non-finite terminal velocity");
    auto output = EvaluationResult::success();
    output.objectives.push_back({config_.objective_id, velocity, true});
    for (const CoilgunMetricConstraint& constraint : config_.constraints) {
        const double value = metric_value(constraint.metric, result, armature, excitations);
        if (!std::isfinite(value))
            return EvaluationResult::failed("non_finite_result", "simulation returned a non-finite metric");
        auto report = constraint.definition.evaluate(value);
        report.id = constraint.id;
        output.constraints.push_back(std::move(report));
        output.metadata[constraint.id] = number(value);
    }
    const std::pair<const char*, CoilgunMetric> metadata[] = {
        {"terminal_velocity", CoilgunMetric::TerminalVelocity}, {"peak_current", CoilgunMetric::PeakCurrent},
        {"peak_voltage", CoilgunMetric::PeakVoltage}, {"maximum_temperature", CoilgunMetric::MaximumTemperature},
        {"efficiency", CoilgunMetric::Efficiency}, {"energy_loss", CoilgunMetric::EnergyLoss}};
    for (const auto& [name, metric] : metadata) {
        const double value = metric_value(metric, result, armature, excitations);
        if (!std::isfinite(value))
            return EvaluationResult::failed("non_finite_result", "simulation returned a non-finite metric");
        output.metadata[name] = number(value);
    }
    return output;
}

std::vector<components::DrivingCoil> CoilgunOptimizationProblem::make_coils(
    const std::vector<double>& values) const {
    const bool has_geometry_binding = std::any_of(
        config_.bindings.begin(), config_.bindings.end(), [](const auto& binding) {
            switch (binding.parameter) {
            case CoilgunParameter::CoilInnerRadius:
            case CoilgunParameter::CoilOuterRadius:
            case CoilgunParameter::CoilLength:
            case CoilgunParameter::CoilTurns:
            case CoilgunParameter::CoilPosition:
                return true;
            default:
                return false;
            }
        });
    if (!has_geometry_binding && !config_.coils.empty()) return config_.coils;
    if (config_.coil_specs.empty()) return config_.coils;

    std::vector<CoilgunCoilSpec> specs;
    specs = config_.coil_specs;
    auto excitations = config_.excitations;
    auto triggers = config_.triggers;
    auto set_value = [&](const CoilgunVariableBinding& binding, double value) {
        switch (binding.parameter) {
        case CoilgunParameter::CoilInnerRadius: specs[binding.index].inner_radius = value; break;
        case CoilgunParameter::CoilOuterRadius: specs[binding.index].outer_radius = value; break;
        case CoilgunParameter::CoilLength: specs[binding.index].length = value; break;
        case CoilgunParameter::CoilTurns: specs[binding.index].turns = static_cast<int>(std::llround(value)); break;
        case CoilgunParameter::CoilPosition: specs[binding.index].position = value; break;
        case CoilgunParameter::ExcitationVoltage: excitations[binding.index].initial_voltage = value; break;
        case CoilgunParameter::ExcitationCapacitance: excitations[binding.index].capacitance = value; break;
        case CoilgunParameter::TriggerValue: triggers[binding.index].value = value; break;
        default: break;
        }
    };
    for (const auto& binding : config_.bindings) {
        const auto it = std::find_if(schema().variables().begin(), schema().variables().end(),
                                     [&](const auto& variable) { return variable.id == binding.variable_id; });
        set_value(binding, values[static_cast<std::size_t>(it - schema().variables().begin())]);
    }
    std::vector<components::DrivingCoil> coils;
    coils.reserve(specs.size());
    for (const auto& spec : specs) coils.push_back(spec.make_coil());
    return coils;
}

double CoilgunOptimizationProblem::metric_value(
    CoilgunMetric metric, const simulation::MultiStageResult& result,
    const components::Armature& armature,
    const std::vector<CoilgunExcitationConfig>& excitations) const {
    switch (metric) {
    case CoilgunMetric::TerminalVelocity: return result.summary.muzzle_velocity;
    case CoilgunMetric::MaximumTemperature: {
        double maximum = physics::T_REFERENCE;
        for (const auto& step : result.history)
            for (double temperature : step.state.filament_temperatures) maximum = std::max(maximum, temperature);
        return maximum;
    }
    case CoilgunMetric::PeakCurrent: return result.summary.peak_coil_current;
    case CoilgunMetric::PeakVoltage: {
        double maximum = 0.0;
        for (const auto& excitation : excitations)
            maximum = std::max(maximum, std::abs(excitation.initial_voltage));
        for (const auto& step : result.history)
            for (double voltage : step.cap_voltages) maximum = std::max(maximum, std::abs(voltage));
        return maximum;
    }
    case CoilgunMetric::Efficiency: return result.summary.efficiency;
    case CoilgunMetric::EnergyLoss: {
        double input = 0.0;
        for (const auto& excitation : excitations)
            input += 0.5 * excitation.capacitance * excitation.initial_voltage * excitation.initial_voltage;
        const double kinetic = 0.5 * armature.mass() * result.summary.muzzle_velocity * result.summary.muzzle_velocity;
        return input - kinetic;
    }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

EvaluationResult CoilgunOptimizationProblem::evaluate_cpu(const CandidateVariables& candidate) const {
    if (candidate.values.size() != schema().size())
        return invalid_result("variable_count", "candidate variable count does not match schema");
    if (!finite_vector(candidate.values))
        return invalid_result("non_finite_variable", "candidate variables must be finite");
    CandidateVariables decoded;
    try { decoded = schema().decode(candidate); }
    catch (const std::exception& error) { return invalid_result("variable_decode", error.what()); }

    auto armature = *config_.armature;
    auto values = decoded.values;
    for (const auto& binding : config_.bindings) {
        const auto it = std::find_if(schema().variables().begin(), schema().variables().end(),
                                     [&](const auto& variable) { return variable.id == binding.variable_id; });
        const double value = values[static_cast<std::size_t>(it - schema().variables().begin())];
        if (binding.parameter == CoilgunParameter::ArmaturePosition)
            armature.update_position(value - armature.position());
        else if (binding.parameter == CoilgunParameter::ArmatureVelocity)
            armature.set_velocity(value);
    }
    try {
        auto coils = make_coils(values);
        auto effective_excitations = config_.excitations;
        auto effective_triggers = config_.triggers;
        for (const auto& binding : config_.bindings) {
            const auto it = std::find_if(schema().variables().begin(), schema().variables().end(),
                                         [&](const auto& variable) { return variable.id == binding.variable_id; });
            const double value = values[static_cast<std::size_t>(it - schema().variables().begin())];
            if (binding.parameter == CoilgunParameter::ExcitationVoltage)
                effective_excitations[binding.index].initial_voltage = value;
            else if (binding.parameter == CoilgunParameter::ExcitationCapacitance)
                effective_excitations[binding.index].capacitance = value;
            else if (binding.parameter == CoilgunParameter::TriggerValue)
                effective_triggers[binding.index].value = value;
        }
        std::vector<std::unique_ptr<simulation::Excitation>> excitations;
        for (const auto& specification : effective_excitations) {
            if (!std::isfinite(specification.initial_voltage) || !std::isfinite(specification.capacitance) || specification.capacitance <= 0.0)
                return invalid_result("invalid_excitation", "excitation parameters must be finite and capacitance positive");
            if (specification.crowbar)
                excitations.push_back(std::make_unique<simulation::CrowbarExcitation>(specification.initial_voltage, specification.capacitance));
            else
                excitations.push_back(std::make_unique<simulation::CapacitorExcitation>(specification.initial_voltage, specification.capacitance));
        }
        const auto simulated_armature = armature;
        simulation::MultiStageSim<simulation::EulerStepper> sim(
            std::move(coils), std::move(armature), std::move(excitations), std::move(effective_triggers),
            config_.dt, config_.enable_thermal, config_.optimization_level);
        const auto& result = sim.run(config_.termination);
        return result_from_simulation(result, simulated_armature, effective_excitations);
    } catch (const std::exception& error) {
        return EvaluationResult::failed("simulation_exception", error.what());
    } catch (...) {
        return EvaluationResult::failed("simulation_exception", "unknown simulation exception");
    }
}

EvaluationResult CoilgunOptimizationProblem::evaluate(const CandidateVariables& variables) const {
    return evaluate_cpu(variables);
}

EvaluationResult CoilgunOptimizationProblem::evaluate(const CandidateVariables& variables,
                                                       const EvaluationContext&) const {
    return evaluate_cpu(variables);
}

void CoilgunOptimizationProblem::set_gpu_batch_evaluator(GpuBatchEvaluator evaluator) {
    std::lock_guard lock(gpu_callback_mutex_);
    gpu_evaluator_ = std::move(evaluator);
    gpu_callback_identity_ = ++next_gpu_callback_identity_;
}

void CoilgunOptimizationProblem::clear_gpu_batch_evaluator() {
    std::lock_guard lock(gpu_callback_mutex_);
    gpu_evaluator_ = {};
    gpu_callback_identity_ = ++next_gpu_callback_identity_;
}

bool CoilgunOptimizationProblem::has_gpu_batch_evaluator() const noexcept {
    std::lock_guard lock(gpu_callback_mutex_);
    return static_cast<bool>(gpu_evaluator_);
}

std::string CoilgunOptimizationProblem::canonical_fingerprint() const {
    std::lock_guard lock(gpu_callback_mutex_);
    return canonical_fingerprint(static_cast<bool>(gpu_evaluator_), gpu_callback_identity_);
}

std::string CoilgunOptimizationProblem::canonical_fingerprint(
    bool has_gpu_evaluator, std::uint64_t callback_identity) const {
    CanonicalWriter writer;
    writer.string("coilgun-optimization-problem-result-v2");

    writer.size(schema().variables().size());
    for (const auto& variable : schema().variables()) {
        writer.string(variable.id);
        enum_value(writer, variable.type);
        writer.real(variable.lower_bound);
        writer.real(variable.upper_bound);
        writer.size(variable.enum_values.size());
        for (const auto& value : variable.enum_values) writer.string(value);
    }

    writer.size(config_.bindings.size());
    for (const auto& binding : config_.bindings) {
        writer.string(binding.variable_id);
        enum_value(writer, binding.parameter);
        writer.size(binding.index);
    }

    writer.size(config_.coils.size());
    for (const auto& coil : config_.coils) write_coil(writer, coil);
    writer.size(config_.coil_specs.size());
    for (const auto& coil : config_.coil_specs) write_coil_spec(writer, coil);

    writer.boolean(config_.armature.has_value());
    if (config_.armature) {
        const auto& armature = *config_.armature;
        writer.real(armature.inner_radius());
        writer.real(armature.outer_radius());
        writer.real(armature.length());
        writer.real(armature.resistivity());
        writer.real(armature.material_density());
        writer.real(armature.velocity());
        writer.real(armature.mass());
        writer.integer(armature.axial_filaments());
        writer.integer(armature.radial_filaments());
        writer.real(armature.position());
        enum_value(writer, armature.material());
        writer.boolean(armature.force_exact_self_inductance());
        writer.size(armature.resistances().size());
        for (const double value : armature.resistances()) writer.real(value);
        writer.size(armature.inductances().size());
        for (const double value : armature.inductances()) writer.real(value);
        writer.size(armature.masses().size());
        for (const double value : armature.masses()) writer.real(value);
    }

    writer.size(config_.excitations.size());
    for (const auto& excitation : config_.excitations) {
        writer.real(excitation.initial_voltage);
        writer.real(excitation.capacitance);
        writer.boolean(excitation.crowbar);
    }
    writer.size(config_.triggers.size());
    for (const auto& trigger : config_.triggers) {
        enum_value(writer, trigger.mode);
        writer.real(trigger.value);
    }

    writer.string(config_.objective_id);
    writer.size(config_.constraints.size());
    for (const auto& constraint : config_.constraints) {
        writer.string(constraint.id);
        enum_value(writer, constraint.metric);
        writer.string(constraint.definition.id);
        enum_value(writer, constraint.definition.kind);
        enum_value(writer, constraint.definition.relation);
        writer.real(constraint.definition.lower_bound);
        writer.real(constraint.definition.upper_bound);
        writer.real(constraint.definition.scale);
        writer.integer(constraint.definition.priority);
    }
    writer.real(config_.dt);
    writer.boolean(config_.enable_thermal);
    enum_value(writer, config_.optimization_level);
    writer.integer(config_.termination.max_steps);
    writer.integer(config_.termination.velocity_decay_steps);
    writer.real(config_.termination.accel_threshold);
    writer.boolean(config_.termination.enable_velocity_check);
    writer.boolean(config_.termination.enable_bound_check);
    writer.real(config_.termination.barrel_end_position);

    writer.boolean(has_gpu_evaluator);
    writer.size(callback_identity);
    return std::move(writer).finish();
}

EvaluationCacheIdentity CoilgunOptimizationProblem::cache_identity() const {
    return {"coilgun.optimization.problem:" + canonical_fingerprint(), "2"};
}

BatchEvaluationSnapshot CoilgunOptimizationProblem::make_evaluation_snapshot() {
    GpuBatchEvaluator gpu_evaluator;
    std::uint64_t callback_identity = 0;
    {
        std::lock_guard lock(gpu_callback_mutex_);
        gpu_evaluator = gpu_evaluator_;
        callback_identity = gpu_callback_identity_;
    }
    const auto identity = EvaluationCacheIdentity{
        "coilgun.optimization.problem:" +
            canonical_fingerprint(static_cast<bool>(gpu_evaluator), callback_identity),
        "2"};
    return {identity, [this, gpu_evaluator = std::move(gpu_evaluator)](
                           const std::vector<CandidateVariables>& variables,
                           const EvaluationContext& context) {
        return evaluate_batch_with_callback(variables, context, gpu_evaluator);
    }};
}

std::vector<EvaluationResult> CoilgunOptimizationProblem::evaluate_batch(
    const std::vector<CandidateVariables>& variables, const EvaluationContext& context) {
    return static_cast<const CoilgunOptimizationProblem&>(*this).evaluate_batch(variables, context);
}

std::vector<EvaluationResult> CoilgunOptimizationProblem::evaluate_batch(
    const std::vector<CandidateVariables>& variables, const EvaluationContext& context) const {
    GpuBatchEvaluator gpu_evaluator;
    {
        std::lock_guard lock(gpu_callback_mutex_);
        gpu_evaluator = gpu_evaluator_;
    }
    return evaluate_batch_with_callback(variables, context, std::move(gpu_evaluator));
}

std::vector<EvaluationResult> CoilgunOptimizationProblem::evaluate_batch_with_callback(
    const std::vector<CandidateVariables>& variables, const EvaluationContext& context,
    GpuBatchEvaluator gpu_evaluator) const {
    {
        std::lock_guard lock(statistics_mutex_);
        last_batch_used_fallback_ = false;
    }
    if (variables.empty()) return {};
    if (gpu_evaluator) {
        try {
            auto results = gpu_evaluator(variables, context);
            if (results.size() == variables.size()) {
                bool malformed = false;
                for (const auto& result : results) {
                    if (result.status == EvaluationStatus::Success &&
                        !valid_gpu_success(result, config_))
                        malformed = true;
                    if (result.status == EvaluationStatus::Unevaluated)
                        malformed = true;
                }
                if (!malformed) return results;
            }
        } catch (...) {
        }
        {
            std::lock_guard lock(statistics_mutex_);
            last_batch_used_fallback_ = true;
            ++statistics_.fallbacks;
        }
        if (context.statistics) context.statistics->add_fallbacks();
    }
    std::vector<EvaluationResult> results;
    results.reserve(variables.size());
    for (const auto& candidate : variables) results.push_back(evaluate_cpu(candidate));
    return results;
}

std::optional<std::vector<EvaluationResult>> CoilgunOptimizationProblem::evaluate_batch_const(
    const std::vector<CandidateVariables>& variables, const EvaluationContext& context) const {
    return evaluate_batch(variables, context);
}

std::optional<EvaluationStatistics> CoilgunOptimizationProblem::statistics_snapshot() const {
    std::lock_guard lock(statistics_mutex_);
    return statistics_;
}

bool CoilgunOptimizationProblem::last_batch_used_fallback() const noexcept {
    std::lock_guard lock(statistics_mutex_);
    return last_batch_used_fallback_;
}

} // namespace coilgun::optimization
