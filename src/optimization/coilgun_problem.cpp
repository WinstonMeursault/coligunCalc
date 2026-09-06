#include "coilgun/optimization/coilgun_problem.hpp"

#include "coilgun/physics/constants.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace coilgun::optimization {
namespace {

bool finite_vector(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); });
}

std::string number(double value) {
    std::ostringstream stream;
    stream << std::setprecision(17) << value;
    return stream.str();
}

} // namespace

components::DrivingCoil CoilgunCoilSpec::make_coil() const {
    return components::DrivingCoil(inner_radius, outer_radius, length, turns, resistivity,
                                    wire_area, fill_factor, position, force_exact_self_inductance);
}

CoilgunOptimizationProblem::CoilgunOptimizationProblem(VariableSchema schema, Config config)
    : schema_(std::move(schema)), config_(std::move(config)) {
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
        for (const auto& variable : schema_.variables()) found |= variable.id == binding.variable_id;
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
        constraint.definition.validate();
    }
}

EvaluationResult CoilgunOptimizationProblem::invalid_result(const std::string& code,
                                                             const std::string& message) const {
    return EvaluationResult::invalid(code, message);
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
        const auto it = std::find_if(schema_.variables().begin(), schema_.variables().end(),
                                     [&](const auto& variable) { return variable.id == binding.variable_id; });
        set_value(binding, values[static_cast<std::size_t>(it - schema_.variables().begin())]);
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
    if (candidate.values.size() != schema_.size())
        return invalid_result("variable_count", "candidate variable count does not match schema");
    if (!finite_vector(candidate.values))
        return invalid_result("non_finite_variable", "candidate variables must be finite");
    CandidateVariables decoded;
    try { decoded = schema_.decode(candidate); }
    catch (const std::exception& error) { return invalid_result("variable_decode", error.what()); }

    auto armature = *config_.armature;
    auto values = decoded.values;
    for (const auto& binding : config_.bindings) {
        const auto it = std::find_if(schema_.variables().begin(), schema_.variables().end(),
                                     [&](const auto& variable) { return variable.id == binding.variable_id; });
        const double value = values[static_cast<std::size_t>(it - schema_.variables().begin())];
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
            const auto it = std::find_if(schema_.variables().begin(), schema_.variables().end(),
                                         [&](const auto& variable) { return variable.id == binding.variable_id; });
            const double value = values[static_cast<std::size_t>(it - schema_.variables().begin())];
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
        const double velocity = metric_value(CoilgunMetric::TerminalVelocity, result, simulated_armature, effective_excitations);
        if (!std::isfinite(velocity)) return EvaluationResult::failed("non_finite_result", "simulation returned a non-finite terminal velocity");
        auto output = EvaluationResult::success();
        output.objectives.push_back({config_.objective_id, velocity, true});
        const auto& original_armature = simulated_armature;
        for (const CoilgunMetricConstraint& constraint : config_.constraints) {
            const double value = metric_value(constraint.metric, result, original_armature, effective_excitations);
            if (!std::isfinite(value)) return EvaluationResult::failed("non_finite_result", "simulation returned a non-finite metric");
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
            const double value = metric_value(metric, result, original_armature, effective_excitations);
            if (!std::isfinite(value)) return EvaluationResult::failed("non_finite_result", "simulation returned a non-finite metric");
            output.metadata[name] = number(value);
        }
        return output;
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

std::vector<EvaluationResult> CoilgunOptimizationProblem::evaluate_batch(
    const std::vector<CandidateVariables>& variables, const EvaluationContext& context) const {
    last_batch_used_fallback_ = false;
    if (variables.empty()) return {};
    if (gpu_evaluator_) {
        try {
            auto results = gpu_evaluator_(variables, context);
            if (results.size() == variables.size()) {
                bool malformed = false;
                for (auto& result : results) {
                    if (result.status == EvaluationStatus::Success) {
                        bool finite = !result.objectives.empty();
                        for (const auto& objective : result.objectives) finite &= std::isfinite(objective.value);
                        for (const auto& constraint : result.constraints)
                            finite &= std::isfinite(constraint.value) && std::isfinite(constraint.violation);
                        if (!finite) result = EvaluationResult::failed(
                            "non_finite_result", "GPU evaluator returned a non-finite result");
                    }
                }
                for (const auto& result : results)
                    malformed |= result.status == EvaluationStatus::Unevaluated;
                if (!malformed) return results;
            }
        } catch (...) {
        }
        last_batch_used_fallback_ = true;
    }
    std::vector<EvaluationResult> results;
    results.reserve(variables.size());
    for (const auto& candidate : variables) results.push_back(evaluate_cpu(candidate));
    return results;
}

} // namespace coilgun::optimization
