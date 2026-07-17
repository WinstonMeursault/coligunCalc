#include <algorithm>
#include <cmath>
#include <string>

#include <doctest/doctest.h>

#include "coilgun/physics/cache.hpp"
#include "coilgun/physics/lookup_tables.hpp"

using coilgun::physics::inductance_shape_factor;
using coilgun::physics::LRUCache;

TEST_CASE("inductance_shape_factor known points") {
    // Reference values from scipy.integrate.quad (Bessel/Struve kernel).
    // Dense table at dq=dp=0.002, bilinear interpolation error ~1e-5.
    CHECK(inductance_shape_factor(1.0, 2.0)
          == doctest::Approx(0.32284315).epsilon(1e-6));
    CHECK(inductance_shape_factor(0.5, 2.5)
          == doctest::Approx(0.23433846).epsilon(1e-6));
    CHECK(inductance_shape_factor(0.5, 3.0)
          == doctest::Approx(0.45174081).epsilon(1e-6));
    CHECK(inductance_shape_factor(1.0, 2.5)
          == doctest::Approx(0.80910223).epsilon(1e-6));
}

TEST_CASE("inductance_shape_factor bilinear interpolation") {
    double T00 = inductance_shape_factor(0.5, 2.5);
    double T10 = inductance_shape_factor(1.0, 2.5);
    double T01 = inductance_shape_factor(0.5, 3.0);
    double T11 = inductance_shape_factor(1.0, 3.0);

    double T_interp = inductance_shape_factor(0.75, 2.75);

    CHECK(T_interp >= std::min({T00, T10, T01, T11}));
    CHECK(T_interp <= std::max({T00, T10, T01, T11}));
}

TEST_CASE("inductance_shape_factor boundary clamping") {
    double val_min = inductance_shape_factor(0.01, 1.04);
    double val_ref_min = inductance_shape_factor(0.05, 1.05);
    CHECK(val_min == doctest::Approx(val_ref_min).epsilon(1e-12));

    double val_max = inductance_shape_factor(5.0, 5.0);
    double val_ref_max = inductance_shape_factor(4.0, 4.0);
    CHECK(val_max == doctest::Approx(val_ref_max).epsilon(1e-12));
}

TEST_CASE("inductance_shape_factor monotonic in q") {
    double p = 2.0;
    CHECK(inductance_shape_factor(2.0, p) > inductance_shape_factor(1.0, p));
    CHECK(inductance_shape_factor(3.0, p) > inductance_shape_factor(2.0, p));
}

TEST_CASE("inductance_shape_factor monotonic in p") {
    double q = 1.0;
    CHECK(inductance_shape_factor(q, 3.0) > inductance_shape_factor(q, 2.0));
}

TEST_CASE("LRUCache put and get") {
    LRUCache<int, std::string, 4> cache;
    cache.put(1, "one");
    cache.put(2, "two");

    std::string val;
    CHECK(cache.get(1, val));
    CHECK(val == "one");
    CHECK(cache.get(2, val));
    CHECK(val == "two");
    CHECK_FALSE(cache.get(3, val));
    CHECK(cache.size() == 2);
}

TEST_CASE("LRUCache eviction") {
    LRUCache<int, int, 3> cache;
    cache.put(1, 10);
    cache.put(2, 20);
    cache.put(3, 30);
    cache.put(4, 40);

    int val;
    CHECK_FALSE(cache.get(1, val));
    CHECK(cache.get(2, val));
    CHECK(val == 20);
    CHECK(cache.get(3, val));
    CHECK(val == 30);
    CHECK(cache.get(4, val));
    CHECK(val == 40);
    CHECK(cache.size() == 3);
}

TEST_CASE("LRUCache lru order") {
    LRUCache<int, int, 3> cache;
    cache.put(1, 10);
    cache.put(2, 20);
    cache.put(3, 30);

    int val;
    cache.get(1, val); // makes 1 MRU, 2 LRU

    cache.put(4, 40);  // should evict 2

    CHECK(cache.size() == 3);
    CHECK(cache.get(1, val));
    CHECK(val == 10);
    CHECK_FALSE(cache.get(2, val));
    CHECK(cache.get(3, val));
    CHECK(val == 30);
    CHECK(cache.get(4, val));
    CHECK(val == 40);
}

TEST_CASE("LRUCache update existing key") {
    LRUCache<int, int, 3> cache;
    cache.put(1, 10);
    cache.put(2, 20);
    cache.put(1, 100);

    int val;
    CHECK(cache.get(1, val));
    CHECK(val == 100);
    CHECK(cache.size() == 2);
}

TEST_CASE("LRUCache clear") {
    LRUCache<int, int, 5> cache;
    cache.put(1, 10);
    cache.put(2, 20);
    cache.clear();
    CHECK(cache.size() == 0);

    int val;
    CHECK_FALSE(cache.get(1, val));
}
