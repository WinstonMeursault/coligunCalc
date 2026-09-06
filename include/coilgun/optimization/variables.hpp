#pragma once

#include "coilgun/optimization/types.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace coilgun::optimization {

enum class VariableType { Continuous, Integer, Enum };

struct VariableSpec {
    std::string id;
    VariableType type = VariableType::Continuous;
    double lower_bound = 0.0;
    double upper_bound = 0.0;
    std::vector<std::string> enum_values;

    static VariableSpec continuous(std::string id, double lower, double upper);
    static VariableSpec integer(std::string id, long long lower, long long upper);
    static VariableSpec enumeration(std::string id, std::vector<std::string> values);

    friend bool operator==(const VariableSpec&, const VariableSpec&) = default;
};

using Variable = VariableSpec;

class VariableSchema {
public:
    VariableSchema() = default;
    explicit VariableSchema(std::vector<VariableSpec> variables);

    [[nodiscard]] std::size_t size() const noexcept { return variables_.size(); }
    [[nodiscard]] bool empty() const noexcept { return variables_.empty(); }
    [[nodiscard]] const std::vector<VariableSpec>& variables() const noexcept { return variables_; }
    [[nodiscard]] const VariableSpec& at(std::size_t index) const;
    [[nodiscard]] CandidateVariables encode(const CandidateVariables& candidate) const;
    [[nodiscard]] CandidateVariables decode(const CandidateVariables& encoded) const;
    [[nodiscard]] CandidateVariables repair(const CandidateVariables& candidate) const;

private:
    std::vector<VariableSpec> variables_;
};

} // namespace coilgun::optimization
