#include "coilgun/optimization/variables.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace coilgun::optimization {
namespace {

double repair_bounded(double value, double lower, double upper) {
    if (std::isnan(value) || value == -std::numeric_limits<double>::infinity()) return lower;
    if (value == std::numeric_limits<double>::infinity()) return upper;
    return std::fmin(upper, std::fmax(lower, value));
}

} // namespace

VariableSpec VariableSpec::continuous(std::string id, double lower, double upper) {
    return VariableSpec{std::move(id), VariableType::Continuous, lower, upper, {}};
}

VariableSpec VariableSpec::integer(std::string id, long long lower, long long upper) {
    return VariableSpec{std::move(id), VariableType::Integer, static_cast<double>(lower),
                        static_cast<double>(upper), {}};
}

VariableSpec VariableSpec::enumeration(std::string id, std::vector<std::string> values) {
    if (values.empty()) throw std::invalid_argument("enum variable must have at least one value");
    return VariableSpec{std::move(id), VariableType::Enum, 0.0,
                        static_cast<double>(values.size() - 1), std::move(values)};
}

VariableSchema::VariableSchema(std::vector<VariableSpec> variables) : variables_(std::move(variables)) {
    std::unordered_set<std::string> ids;
    for (const auto& variable : variables_) {
        if (variable.id.empty()) throw std::invalid_argument("variable ID must not be empty");
        if (!ids.insert(variable.id).second) throw std::invalid_argument("duplicate variable ID: " + variable.id);
        if (!std::isfinite(variable.lower_bound) || !std::isfinite(variable.upper_bound) ||
            variable.lower_bound > variable.upper_bound) {
            throw std::invalid_argument("variable bounds must be finite and ordered");
        }
        if (variable.type == VariableType::Enum &&
            (variable.enum_values.empty() || variable.lower_bound != 0.0 ||
             variable.upper_bound != static_cast<double>(variable.enum_values.size() - 1))) {
            throw std::invalid_argument("invalid enum variable definition");
        }
        if (variable.type != VariableType::Enum && !variable.enum_values.empty()) {
            throw std::invalid_argument("non-enum variable cannot have enum values");
        }
    }
}

const VariableSpec& VariableSchema::at(std::size_t index) const { return variables_.at(index); }

CandidateVariables VariableSchema::repair(const CandidateVariables& candidate) const {
    if (candidate.values.size() != variables_.size())
        throw std::invalid_argument("candidate variable count does not match schema");
    std::vector<double> repaired;
    repaired.reserve(candidate.values.size());
    for (std::size_t i = 0; i < variables_.size(); ++i) {
        const auto& variable = variables_[i];
        const double value = candidate.values[i];
        if (variable.type == VariableType::Continuous) {
            repaired.push_back(repair_bounded(value, variable.lower_bound, variable.upper_bound));
        } else if (variable.type == VariableType::Integer) {
            repaired.push_back(std::round(repair_bounded(value, variable.lower_bound, variable.upper_bound)));
        } else {
            if (!std::isfinite(value) || std::trunc(value) != value)
                throw std::invalid_argument("enum variable requires a finite integer index");
            if (value < variable.lower_bound || value > variable.upper_bound)
                throw std::out_of_range("enum variable index is out of range");
            repaired.push_back(value);
        }
    }
    return CandidateVariables{std::move(repaired)};
}

CandidateVariables VariableSchema::encode(const CandidateVariables& candidate) const { return repair(candidate); }
CandidateVariables VariableSchema::decode(const CandidateVariables& encoded) const { return repair(encoded); }

} // namespace coilgun::optimization
