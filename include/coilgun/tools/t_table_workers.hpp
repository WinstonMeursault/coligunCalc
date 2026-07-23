#pragma once

#include <algorithm>

namespace coilgun::tools {

constexpr unsigned int t_table_worker_count(unsigned int hardware) noexcept {
    const unsigned int usable = hardware > 4 ? hardware - 4 : 1;
    return std::min(usable, 32u);
}

} // namespace coilgun::tools
