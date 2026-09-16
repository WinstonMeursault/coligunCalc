#pragma once

#include "coilgun/optimization/types.hpp"
#include <string>

namespace coilgun::optimization {

struct ObjectiveDefinition {
    std::string id;
    bool maximize = true;
    double scale = 1.0;

    void validate() const;
    double normalize(double value) const;
    double oriented(double value) const;
};

}
