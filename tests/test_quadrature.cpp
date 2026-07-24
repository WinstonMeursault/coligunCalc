#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <vector>

#include "coilgun/physics/quadrature.hpp"
#include "doctest/doctest.h"

using coilgun::physics::gauss_legendre;
using coilgun::physics::gauss_laguerre;
using coilgun::physics::gauss_legendre_cached;
using coilgun::physics::gauss_laguerre_cached;

namespace {

template <typename Rule>
void check_cached_identity(Rule rule, const std::vector<int>& orders) {
    for (int n : orders) {
        const auto& first = rule(n);
        const auto& second = rule(n);

        CHECK(&first == &second);
        CHECK(first.nodes.data() == second.nodes.data());
        CHECK(first.weights.data() == second.weights.data());
    }
}

template <typename Rule>
void check_concurrent_reads(Rule rule, int n) {
    constexpr std::size_t thread_count = 8;
    std::atomic<bool> start{false};
    std::array<const coilgun::physics::QuadratureNodes*, thread_count> addresses{};
    std::array<double, thread_count> checksums{};
    std::array<std::thread, thread_count> threads;

    for (std::size_t thread = 0; thread < thread_count; ++thread) {
        threads[thread] = std::thread([&, thread] {
            while (!start.load(std::memory_order_acquire)) {
                std::this_thread::yield();
            }
            const auto& q = rule(n);
            addresses[thread] = &q;
            for (std::size_t i = 0; i < q.nodes.size(); ++i) {
                checksums[thread] += q.nodes[i] * q.weights[i];
            }
        });
    }

    start.store(true, std::memory_order_release);
    for (auto& thread : threads) {
        thread.join();
    }

    for (std::size_t thread = 1; thread < thread_count; ++thread) {
        CHECK(addresses[thread] == addresses[0]);
        CHECK(checksums[thread] == checksums[0]);
    }
}

} // namespace

TEST_CASE("gauss_legendre(4) integrates x^6 exactly") {
    auto q = gauss_legendre(4);

    CHECK(q.nodes.size() == 4);
    CHECK(q.weights.size() == 4);

    double integral = 0.0;
    for (std::size_t i = 0; i < q.nodes.size(); ++i) {
        double x = q.nodes[i];
        integral += q.weights[i] * std::pow(x, 6);
    }

    // ∫_{-1}^{1} x^6 dx = 2/7 ≈ 0.2857142857142857
    CHECK(integral == doctest::Approx(2.0 / 7.0));
}

TEST_CASE("gauss_legendre(9) integrates x^10 and x^17 exactly") {
    auto q = gauss_legendre(9);
    double integral;

    // ∫_{-1}^{1} x^{10} dx = 2/11 ≈ 0.181818...
    integral = 0.0;
    for (std::size_t i = 0; i < q.nodes.size(); ++i) {
        integral += q.weights[i] * std::pow(q.nodes[i], 10);
    }
    CHECK(integral == doctest::Approx(2.0 / 11.0));

    // ∫_{-1}^{1} x^{16} dx = 2/17 ≈ 0.117647... — 9-point exact up to degree 2n-1=17
    integral = 0.0;
    for (std::size_t i = 0; i < q.nodes.size(); ++i) {
        integral += q.weights[i] * std::pow(q.nodes[i], 16);
    }
    CHECK(integral == doctest::Approx(2.0 / 17.0));
}

TEST_CASE("gauss_laguerre(15) integrates x^2 e^{-x} → Γ(3)=2") {
    auto q = gauss_laguerre(15);

    CHECK(q.nodes.size() == 15);
    CHECK(q.weights.size() == 15);

    double integral = 0.0;
    for (std::size_t i = 0; i < q.nodes.size(); ++i) {
        double x = q.nodes[i];
        integral += q.weights[i] * x * x; // w_i already absorbs e^{-x}
    }

    CHECK(integral == doctest::Approx(2.0));
}

TEST_CASE("gauss_legendre weight sum equals interval length") {
    for (int n : {4, 9, 16, 32}) {
        auto q = gauss_legendre(n);
        double sum = 0.0;
        for (auto w : q.weights) {
            sum += w;
        }
        CHECK(sum == doctest::Approx(2.0)); // interval [-1, 1] has length 2
    }
}

TEST_CASE("gauss_laguerre weight sum equals 1") {
    for (int n : {15, 30}) {
        auto q = gauss_laguerre(n);
        double sum = 0.0;
        for (auto w : q.weights) {
            sum += w;
        }
        CHECK(sum == doctest::Approx(1.0));
    }
}

TEST_CASE("gauss_legendre rejects unsupported n") {
    CHECK_THROWS_AS(gauss_legendre(3), std::invalid_argument);
    CHECK_THROWS_AS(gauss_legendre(0), std::invalid_argument);
    CHECK_THROWS_WITH(gauss_legendre(3),
                      "gauss_legendre: unsupported node count 3 (supported: 4, 9, 16, 32)");
}

TEST_CASE("gauss_laguerre rejects unsupported n") {
    CHECK_THROWS_AS(gauss_laguerre(10), std::invalid_argument);
    CHECK_THROWS_AS(gauss_laguerre(0), std::invalid_argument);
    CHECK_THROWS_WITH(gauss_laguerre(10),
                      "gauss_laguerre: unsupported node count 10 (supported: 15, 30)");
}

TEST_CASE("all supported quadrature rules preserve exact nodes and weights") {
    const std::array<double, 4> gl4_nodes = {
        -0.8611363115940526, -0.3399810435848563,
         0.3399810435848563,  0.8611363115940526};
    const std::array<double, 4> gl4_weights = {
        0.3478548451374538, 0.6521451548625461,
        0.6521451548625461, 0.3478548451374538};
    const auto& gl4 = gauss_legendre(4);
    CHECK(gl4.nodes.size() == gl4_nodes.size());
    CHECK(gl4.weights.size() == gl4_weights.size());
    CHECK(std::equal(gl4.nodes.begin(), gl4.nodes.end(), gl4_nodes.begin()));
    CHECK(std::equal(gl4.weights.begin(), gl4.weights.end(), gl4_weights.begin()));

    const std::array<double, 9> gl9_nodes = {
        -0.9681602395076261, -0.8360311073266358,
        -0.6133714327005904, -0.3242534234038089, 0.0,
         0.3242534234038089,  0.6133714327005904,
         0.8360311073266358,  0.9681602395076261};
    const std::array<double, 9> gl9_weights = {
        0.0812743883615744, 0.1806481606948574,
        0.2606106964029354, 0.3123470770400029,
        0.3302393550012598, 0.3123470770400029,
        0.2606106964029354, 0.1806481606948574,
        0.0812743883615744};
    const auto& gl9 = gauss_legendre(9);
    CHECK(gl9.nodes.size() == gl9_nodes.size());
    CHECK(gl9.weights.size() == gl9_weights.size());
    CHECK(std::equal(gl9.nodes.begin(), gl9.nodes.end(), gl9_nodes.begin()));
    CHECK(std::equal(gl9.weights.begin(), gl9.weights.end(), gl9_weights.begin()));

    const std::array<double, 16> gl16_nodes = {
        -0.9894009349916499, -0.9445750230732326,
        -0.8656312023878318, -0.7554044083550030,
        -0.6178762444026438, -0.4580167776572274,
        -0.2816035507792589, -0.0950125098376374,
         0.0950125098376374,  0.2816035507792589,
         0.4580167776572274,  0.6178762444026438,
         0.7554044083550030,  0.8656312023878318,
         0.9445750230732326,  0.9894009349916499};
    const std::array<double, 16> gl16_weights = {
        0.0271524594117541, 0.0622535239386479,
        0.0951585116824928, 0.1246289712555339,
        0.1495959888165767, 0.1691565193950025,
        0.1826034150449236, 0.1894506104550685,
        0.1894506104550685, 0.1826034150449236,
        0.1691565193950025, 0.1495959888165767,
        0.1246289712555339, 0.0951585116824928,
        0.0622535239386479, 0.0271524594117541};
    const auto& gl16 = gauss_legendre(16);
    CHECK(gl16.nodes.size() == gl16_nodes.size());
    CHECK(gl16.weights.size() == gl16_weights.size());
    CHECK(std::equal(gl16.nodes.begin(), gl16.nodes.end(), gl16_nodes.begin()));
    CHECK(std::equal(gl16.weights.begin(), gl16.weights.end(), gl16_weights.begin()));

    const std::array<double, 32> gl32_nodes = {
        -0.9972638618494816, -0.9856115115452684,
        -0.9647622555875064, -0.9349060759377397,
        -0.8963211557660521, -0.8493676137325701,
        -0.7944837959679424, -0.7321821187402897,
        -0.6630442669302152, -0.5877157572407623,
        -0.5068999089322294, -0.4213512761306353,
        -0.3318686022821276, -0.2392873622521371,
        -0.1444719615827965, -0.0483076656877383,
         0.0483076656877383,  0.1444719615827965,
         0.2392873622521371,  0.3318686022821276,
         0.4213512761306353,  0.5068999089322294,
         0.5877157572407623,  0.6630442669302152,
         0.7321821187402897,  0.7944837959679424,
         0.8493676137325701,  0.8963211557660521,
         0.9349060759377397,  0.9647622555875064,
         0.9856115115452684,  0.9972638618494816};
    const std::array<double, 32> gl32_weights = {
        0.0070186100094701, 0.0162743947309057,
        0.0253920653092621, 0.0342738629130214,
        0.0428358980222267, 0.0509980592623762,
        0.0586840934785355, 0.0658222227763618,
        0.0723457941088485, 0.0781938957870703,
        0.0833119242269468, 0.0876520930044038,
        0.0911738786957639, 0.0938443990808046,
        0.0956387200792749, 0.0965400885147278,
        0.0965400885147278, 0.0956387200792749,
        0.0938443990808046, 0.0911738786957639,
        0.0876520930044038, 0.0833119242269468,
        0.0781938957870703, 0.0723457941088485,
        0.0658222227763618, 0.0586840934785355,
        0.0509980592623762, 0.0428358980222267,
        0.0342738629130214, 0.0253920653092621,
        0.0162743947309057, 0.0070186100094701};
    const auto& gl32 = gauss_legendre(32);
    CHECK(gl32.nodes.size() == gl32_nodes.size());
    CHECK(gl32.weights.size() == gl32_weights.size());
    CHECK(std::equal(gl32.nodes.begin(), gl32.nodes.end(), gl32_nodes.begin()));
    CHECK(std::equal(gl32.weights.begin(), gl32.weights.end(), gl32_weights.begin()));

    const std::array<double, 15> gl15_nodes = {
        0.0933078120172819, 0.4926917403018839,
        1.2155954120709496, 2.2699495262037437,
        3.6676227217514370, 5.4253366274135540,
        7.5659162266130680, 10.1202285680191140,
        13.1302824821757240, 16.6544077083299600,
        20.7764788994487700, 25.6238942267287800,
        31.4075191697539400, 38.5306833064860100,
        48.0260855726857940};
    const std::array<double, 15> gl15_weights = {
        0.2182348859400943, 0.3422101779228796,
        0.2630275779416778, 0.1264258181059299,
        0.0402068649210004, 0.0085638778036118,
        0.0012124361472142, 0.0001116743923443,
        0.0000064599267620, 0.0000002226316907,
        0.0000000042274304, 0.0000000000392190,
        0.0000000000001457, 0.0000000000000001,
        0.0000000000000000};
    const auto& gl15 = gauss_laguerre(15);
    CHECK(gl15.nodes.size() == gl15_nodes.size());
    CHECK(gl15.weights.size() == gl15_weights.size());
    CHECK(std::equal(gl15.nodes.begin(), gl15.nodes.end(), gl15_nodes.begin()));
    CHECK(std::equal(gl15.weights.begin(), gl15.weights.end(), gl15_weights.begin()));

    const std::array<double, 30> gl30_nodes = {
        0.0474071805408053, 0.2499239167531594,
        0.6148334543927684, 1.1431958256661015,
        1.8364545546225723, 2.6965218745572160,
        3.7258145077795093, 4.9272937658498820,
        6.3045155909650740, 7.8616932933702600,
        9.6037759854792630, 11.5365465979561400,
        13.6667446930642350, 16.0022211889810680,
        18.5521348401431500, 21.3272043217831280,
        24.3400357645326930, 27.6055547967809600,
        31.1415867011112370, 34.9696520082490700,
        39.1160849490678900, 43.6136529084848300,
        48.5039861638042000, 53.8413854065075060,
        59.6991218592355000, 66.1806177944384900,
        73.4412385955598800, 81.7368105067276800,
        91.5564665225368400, 104.1575244310588900};
    const std::array<double, 30> gl30_weights = {
        0.1160440860204389, 0.2208511247506771,
        0.2413998275878537, 0.1946367684464171,
        0.1237284159668765, 0.0636787803689866,
        0.0268604752733797, 0.0093380708816039,
        0.0026806968913368, 0.0006351291219409,
        0.0001239074599069, 0.0000198287884390,
        0.0000025893509291, 0.0000002740942841,
        0.0000000233283117, 0.0000000015807456,
        0.0000000000842748, 0.0000000000034852,
        0.0000000000001099, 0.0000000000000026,
        0.0000000000000000, 0.0000000000000000,
        0.0000000000000000, 0.0000000000000000,
        0.0000000000000000, 0.0000000000000000,
        0.0000000000000000, 0.0000000000000000,
        0.0000000000000000, 0.0000000000000000};
    const auto& gl30 = gauss_laguerre(30);
    CHECK(gl30.nodes.size() == gl30_nodes.size());
    CHECK(gl30.weights.size() == gl30_weights.size());
    CHECK(std::equal(gl30.nodes.begin(), gl30.nodes.end(), gl30_nodes.begin()));
    CHECK(std::equal(gl30.weights.begin(), gl30.weights.end(), gl30_weights.begin()));
}

TEST_CASE("quadrature rules are immutable cached objects") {
    static_assert(std::is_const_v<std::remove_reference_t<decltype(gauss_legendre_cached(4))>>);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(gauss_laguerre_cached(15))>>);

    check_cached_identity(
        gauss_legendre_cached,
        {4, 9, 16, 32});
    check_cached_identity(
        gauss_laguerre_cached,
        {15, 30});

    const auto& original = gauss_legendre_cached(9);
    auto copy = gauss_legendre(9);
    copy.nodes[0] = 0.0;
    copy.weights[0] = 0.0;
    CHECK(original.nodes[0] != 0.0);
    CHECK(original.weights[0] != 0.0);
}

TEST_CASE("quadrature rules support concurrent read access") {
    check_concurrent_reads(
        gauss_legendre_cached, 32);
    check_concurrent_reads(
        gauss_laguerre_cached, 30);
}

TEST_CASE("quadrature cache repeated-call benchmark") {
    constexpr std::size_t repetitions = 100000;
    volatile std::size_t sink = 0;

    const auto cached_start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repetitions; ++i) {
        const auto& q = gauss_legendre_cached(9);
        sink += q.nodes.size() + q.weights.size();
    }
    const auto cached_stop = std::chrono::steady_clock::now();

    const auto copied_start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repetitions; ++i) {
        auto q = gauss_legendre(9);
        sink += q.nodes.size() + q.weights.size();
    }
    const auto copied_stop = std::chrono::steady_clock::now();

    const double cached_ms = std::chrono::duration<double, std::milli>(
        cached_stop - cached_start).count();
    const double copied_ms = std::chrono::duration<double, std::milli>(
        copied_stop - copied_start).count();
    std::cout << "quadrature benchmark repetitions=" << repetitions
              << " cached_ms=" << cached_ms
              << " copied_ms=" << copied_ms << '\n';
    CHECK(sink == 4 * repetitions * 9);
}
