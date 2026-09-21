#include <doctest/doctest.h>

#include "coilgun/optimization/evaluator.hpp"
#include "coilgun/optimization/genetic_optimizer.hpp"

#include <atomic>
#include <barrier>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <thread>
#include <vector>

using namespace coilgun::optimization;

namespace {
VariableSchema schema() {
    return VariableSchema({VariableSpec::continuous("x", 0.0, 1.0)});
}

VariableSchema fixed_schema() {
    return VariableSchema({VariableSpec::continuous("x", 0.0, 0.0)});
}

OptimizationConfig config(std::uint64_t seed, std::size_t population = 1) {
    OptimizationConfig value;
    value.population_size = population;
    value.max_generations = 1;
    value.elite_count = 0;
    value.crossover_rate = 0.0;
    value.mutation_rate = 0.0;
    value.random_seed = seed;
    return value;
}

class ValueEvaluator final : public BatchEvaluator {
public:
    std::atomic<std::size_t> calls{0};

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        calls.fetch_add(1, std::memory_order_relaxed);
        std::vector<EvaluationResult> results;
        results.reserve(values.size());
        for (const auto& value : values) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"score", value.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    }
};

class InterleavedEvaluator final : public BatchEvaluator {
public:
    explicit InterleavedEvaluator(std::barrier<>& barrier) : barrier_(barrier) {}

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext& context) override {
        barrier_.arrive_and_wait();
        if (!context.statistics) return {};
        const auto cache_hits = context.seed == 11 ? 2U : 5U;
        const auto fallbacks = context.seed == 11 ? 1U : 3U;
        context.statistics->add_cache_hits(cache_hits);
        context.statistics->add_fallbacks(fallbacks);
        std::vector<EvaluationResult> results;
        for (const auto& value : values) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"score", value.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    }

private:
    std::barrier<>& barrier_;
};

class AggregateBlockingEvaluator final : public BatchEvaluator {
public:
    explicit AggregateBlockingEvaluator(std::barrier<>& barrier) : barrier_(barrier) {}

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                 const EvaluationContext&) override {
        barrier_.arrive_and_wait();
        std::vector<EvaluationResult> results;
        for (const auto& value : values) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"score", value.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    }

private:
    std::barrier<>& barrier_;
};
}

TEST_CASE("sequential runs sharing a cached evaluator keep independent run snapshots") {
    auto delegate = std::make_shared<ValueEvaluator>();
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    auto cached = std::make_shared<CachedBatchEvaluator>(delegate, cache);

    const auto first = GeneticOptimizer(schema(), cached, config(17)).optimize();
    const auto first_snapshot = first.statistics;
    const auto second = GeneticOptimizer(schema(), cached, config(17)).optimize();

    CHECK(first.statistics.cache_hits == first_snapshot.cache_hits);
    CHECK(first_snapshot.cache_hits == 0);
    CHECK(second.statistics.cache_hits == 1);
    CHECK(first.statistics.evaluations == 1);
    CHECK(second.statistics.evaluations == 1);
    CHECK(delegate->calls == 1);
}

TEST_CASE("interleaved concurrent runs receive only their own collector metrics") {
    std::barrier barrier(2);
    auto evaluator = std::make_shared<InterleavedEvaluator>(barrier);
    OptimizationResult first;
    OptimizationResult second;

    std::thread first_thread([&] {
        first = GeneticOptimizer(schema(), evaluator, config(11)).optimize();
    });
    std::thread second_thread([&] {
        second = GeneticOptimizer(schema(), evaluator, config(22)).optimize();
    });
    first_thread.join();
    second_thread.join();

    CHECK(first.statistics.seed == 11);
    CHECK(first.statistics.cache_hits == 2);
    CHECK(first.statistics.gpu_fallbacks == 1);
    CHECK(second.statistics.seed == 22);
    CHECK(second.statistics.cache_hits == 5);
    CHECK(second.statistics.gpu_fallbacks == 3);
}

TEST_CASE("optimizer core statistics remain complete without evaluator aggregates") {
    auto evaluator = std::make_shared<ValueEvaluator>();
    TerminationConfig termination;
    termination.max_evaluations = 1;
    const auto result = GeneticOptimizer(schema(), evaluator, config(31, 3), termination).optimize();

    CHECK(result.statistics.seed == 31);
    CHECK(result.statistics.evaluations == 1);
    CHECK(result.statistics.successful_evaluations == 1);
    CHECK(result.statistics.failed_evaluations == 0);
    CHECK(result.statistics.skipped_due_to_budget == 2);
    CHECK(result.statistics.generations == 1);
    CHECK(result.statistics.elapsed_seconds >= 0.0);
}

TEST_CASE("collector keeps all duration fields finite and non-negative") {
    EvaluationStatisticsCollector collector;
    EvaluationStatistics delta;
    delta.elapsed_seconds = -1.0;
    delta.gpu_transfer_seconds = std::numeric_limits<double>::quiet_NaN();
    delta.gpu_kernel_seconds = -2.0;
    delta.gpu_elapsed_seconds = 0.25;
    collector.add(delta);
    collector.add_elapsed_seconds(-3.0);
    collector.add_gpu_transfer_seconds(std::numeric_limits<double>::infinity());
    collector.add_gpu_kernel_seconds(0.5);

    const auto snapshot = collector.snapshot();
    CHECK(snapshot.elapsed_seconds == 0.0);
    CHECK(snapshot.gpu_transfer_seconds == 0.0);
    CHECK(snapshot.gpu_kernel_seconds == doctest::Approx(0.5));
    CHECK(snapshot.gpu_elapsed_seconds == doctest::Approx(0.25));
    CHECK(std::isfinite(snapshot.elapsed_seconds));
    CHECK(std::isfinite(snapshot.gpu_transfer_seconds));
    CHECK(std::isfinite(snapshot.gpu_kernel_seconds));
    CHECK(std::isfinite(snapshot.gpu_elapsed_seconds));
}

TEST_CASE("shared in-memory cache supports concurrent hits") {
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    auto result = EvaluationResult::success();
    result.objectives.push_back({"score", 1.0, true});
    const auto key = make_cache_key(CandidateVariables{{1.0}}, EvaluationContext{});
    cache->put(key, result);

    std::barrier barrier(2);
    std::optional<EvaluationResult> first;
    std::optional<EvaluationResult> second;
    std::thread first_thread([&] {
        barrier.arrive_and_wait();
        first = cache->get(key);
    });
    std::thread second_thread([&] {
        barrier.arrive_and_wait();
        second = cache->get(key);
    });
    first_thread.join();
    second_thread.join();

    REQUIRE(first);
    REQUIRE(second);
    CHECK(first->objectives.front().value == 1.0);
    CHECK(second->objectives.front().value == 1.0);
}

TEST_CASE("concurrent optimizer runs share cached storage without sharing hit counters") {
    auto delegate = std::make_shared<ValueEvaluator>();
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    auto cached = std::make_shared<CachedBatchEvaluator>(delegate, cache);
    auto value = EvaluationResult::success();
    value.objectives.push_back({"score", 0.0, true});
    cache->put(make_cache_key(CandidateVariables{{0.0}}, EvaluationContext{44, false}), value);

    std::barrier start(3);
    OptimizationResult first;
    OptimizationResult second;
    std::thread first_thread([&] {
        start.arrive_and_wait();
        first = GeneticOptimizer(fixed_schema(), cached, config(44)).optimize();
    });
    std::thread second_thread([&] {
        start.arrive_and_wait();
        second = GeneticOptimizer(fixed_schema(), cached, config(44)).optimize();
    });
    start.arrive_and_wait();
    first_thread.join();
    second_thread.join();

    CHECK(first.statistics.cache_hits == 1);
    CHECK(second.statistics.cache_hits == 1);
    CHECK(delegate->calls == 0);
}

TEST_CASE("shared evaluator aggregate snapshots are safe during concurrent batches") {
    std::barrier barrier(3);
    auto delegate = std::make_shared<AggregateBlockingEvaluator>(barrier);
    auto tracked = std::make_shared<StatisticsBatchEvaluator>(delegate);
    std::thread first([&] { tracked->evaluate_batch({CandidateVariables{{0.1}}}, {}); });
    std::thread second([&] { tracked->evaluate_batch({CandidateVariables{{0.2}}}, {}); });
    barrier.arrive_and_wait();
    const auto snapshot_while_running = tracked->statistics_snapshot();
    first.join();
    second.join();

    REQUIRE(snapshot_while_running);
    CHECK(snapshot_while_running->evaluations <= 2);
    CHECK(tracked->statistics().evaluations == 2);
    CHECK(tracked->statistics_snapshot()->successful_evaluations == 2);
}
