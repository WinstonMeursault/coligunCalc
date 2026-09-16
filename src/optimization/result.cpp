#include "coilgun/optimization/result.hpp"
#include "coilgun/optimization/selectors.hpp"

namespace coilgun::optimization {

std::optional<Candidate> OptimizationResult::select_representative(
    const RepresentativeSelector& selector) const {
    return selector.select(*this);
}

}
