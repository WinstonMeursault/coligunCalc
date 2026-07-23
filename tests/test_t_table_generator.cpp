#include <doctest/doctest.h>

#include "coilgun/tools/t_table_workers.hpp"

TEST_CASE("T table worker count does not underflow") {
    CHECK(coilgun::tools::t_table_worker_count(0) == 1);
    CHECK(coilgun::tools::t_table_worker_count(1) == 1);
    CHECK(coilgun::tools::t_table_worker_count(2) == 1);
    CHECK(coilgun::tools::t_table_worker_count(3) == 1);
    CHECK(coilgun::tools::t_table_worker_count(64) == 32);
}
