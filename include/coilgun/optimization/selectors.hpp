#pragma once

#include "coilgun/optimization/result.hpp"

#include <functional>
#include <string>
#include <vector>

namespace coilgun::optimization {

class RepresentativeSelector {
public:
    virtual ~RepresentativeSelector() = default;
    virtual std::optional<Candidate> select(const OptimizationResult& result) const = 0;
};

class MaxObjective final : public RepresentativeSelector {
public:
    explicit MaxObjective(std::string objective_id) : objective_id_(std::move(objective_id)) {}
    std::optional<Candidate> select(const OptimizationResult& result) const override;
private:
    std::string objective_id_;
};

class MinConstraintViolationMargin final : public RepresentativeSelector {
public:
    std::optional<Candidate> select(const OptimizationResult& result) const override;
};

class IdealPointDistance final : public RepresentativeSelector {
public:
    explicit IdealPointDistance(double distance_power = 2.0) : distance_power_(distance_power) {}
    std::optional<Candidate> select(const OptimizationResult& result) const override;
private:
    double distance_power_;
};

class WeightedScore final : public RepresentativeSelector {
public:
    explicit WeightedScore(std::vector<double> weights) : weights_(std::move(weights)) {}
    std::optional<Candidate> select(const OptimizationResult& result) const override;
private:
    std::vector<double> weights_;
};

class LexicographicObjectives final : public RepresentativeSelector {
public:
    explicit LexicographicObjectives(std::vector<std::string> objective_ids)
        : objective_ids_(std::move(objective_ids)) {}
    std::optional<Candidate> select(const OptimizationResult& result) const override;
private:
    std::vector<std::string> objective_ids_;
};

class CallbackSelector final : public RepresentativeSelector {
public:
    using Callback = std::function<std::optional<Candidate>(const OptimizationResult&)>;
    CallbackSelector() = default;
    explicit CallbackSelector(Callback callback) : callback_(std::move(callback)) {}
    std::optional<Candidate> select(const OptimizationResult& result) const override;
private:
    Callback callback_;
};

using CustomSelector = CallbackSelector;

}
