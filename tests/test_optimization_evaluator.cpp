#include <doctest/doctest.h>

#include "coilgun/optimization/cache.hpp"
#include "coilgun/optimization/evaluator.hpp"
#include "coilgun/optimization/statistics.hpp"

#include <atomic>
#include <barrier>
#include <stdexcept>
#include <thread>

using namespace coilgun::optimization;

namespace {
CandidateVariables cv(double value) { return CandidateVariables(std::vector<double>{value}); }
class IncrementingEvaluator final : public Evaluator {
public:
    EvaluationResult evaluate(const CandidateVariables& variables,
                              const EvaluationContext& context) const override {
        if (variables.values.front() < 0.0) {
            throw std::runtime_error("bad candidate");
        }
        auto result = EvaluationResult::success();
        result.metadata["seed"] = std::to_string(context.seed);
        result.objectives.push_back({"value", variables.values.front(), true});
        return result;
    }
};

class ThrowingBatchEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        if (candidates.front().values.front() < 0.0) {
            throw std::runtime_error("bad batch candidate");
        }
        auto result = EvaluationResult::success();
        result.objectives.push_back({"value", candidates.front().values.front(), true});
        return {result};
    }
};

class CandidateThrowingBatchEvaluator final : public BatchEvaluator {
public:
    std::size_t calls = 0;

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        ++calls;
        for (const auto& candidate : candidates) {
            if (candidate.values.front() < 0.0) throw std::runtime_error("bad candidate");
        }
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for (const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"value", candidate.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    }
};

class UnknownThrowingBatchEvaluator final : public BatchEvaluator {
public:
    std::size_t calls = 0;

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        ++calls;
        for (const auto& candidate : candidates) {
            if (candidate.values.front() < 0.0) throw 42;
        }
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for (const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"value", candidate.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    }
};

class EmptyBatchEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&,
                                                 const EvaluationContext&) override {
        return {};
    }
};

class RecordingBatchEvaluator final : public BatchEvaluator {
public:
    std::size_t calls = 0;
    std::vector<double> batches;

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        ++calls;
        std::vector<EvaluationResult> results;
        for (const auto& candidate : candidates) {
            batches.push_back(candidate.values.front());
            auto result = EvaluationResult::success();
            result.objectives.push_back({"value", candidate.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    }
};

class UnevaluatedBatchEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        return std::vector<EvaluationResult>(candidates.size());
    }
};

class ShortBatchEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        if (candidates.empty()) return {};
        auto result = EvaluationResult::success();
        result.objectives.push_back({"value", candidates.front().values.front(), true});
        return {std::move(result)};
    }
};

class FallbackReportingEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        ++statistics_.fallbacks;
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for (const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"value", candidate.values.front(), true});
            results.push_back(std::move(result));
        }
        return results;
    }

    std::optional<EvaluationStatistics> statistics_snapshot() const override { return statistics_; }

private:
    EvaluationStatistics statistics_;
};

class IdentifiedBatchEvaluator final : public BatchEvaluator {
public:
    IdentifiedBatchEvaluator(EvaluationCacheIdentity identity, std::string objective_id,
                             double objective_value)
        : identity_(std::move(identity)), objective_id_(std::move(objective_id)),
          objective_value_(objective_value) {}

    [[nodiscard]] EvaluationCacheIdentity cache_identity() const override { return identity_; }

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        ++calls;
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for ([[maybe_unused]] const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({objective_id_, objective_value_, true});
            results.push_back(std::move(result));
        }
        return results;
    }

    std::size_t calls = 0;

private:
    EvaluationCacheIdentity identity_;
    std::string objective_id_;
    double objective_value_;
};

class GenerationSwitchBatchEvaluator final : public BatchEvaluator {
public:
    explicit GenerationSwitchBatchEvaluator(std::barrier<>& switch_barrier)
        : switch_barrier_(switch_barrier) {}

    [[nodiscard]] EvaluationCacheIdentity cache_identity() const override {
        const auto generation = generation_.load();
        switch_barrier_.arrive_and_wait();
        return {"generation-switch-" + std::to_string(generation), "v1"};
    }

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        const auto generation = generation_.load();
        switch_barrier_.arrive_and_wait();
        return results_for(candidates, generation);
    }

protected:
    BatchEvaluationSnapshot make_evaluation_snapshot() override {
        const auto generation = generation_.load();
        switch_barrier_.arrive_and_wait();
        return {{"generation-switch-" + std::to_string(generation), "v1"},
                [this, generation](const std::vector<CandidateVariables>& candidates,
                                   const EvaluationContext&) {
                    switch_barrier_.arrive_and_wait();
                    return results_for(candidates, generation);
                }};
    }

public:
    void replace_callback_generation() { generation_.store(2); }

private:
    static std::vector<EvaluationResult> results_for(const std::vector<CandidateVariables>& candidates,
                                                     int generation) {
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for ([[maybe_unused]] const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"generation", static_cast<double>(generation), true});
            results.push_back(std::move(result));
        }
        return results;
    }

    std::barrier<>& switch_barrier_;
    std::atomic<int> generation_{1};
};

class LifetimeBarrierEvaluator final : public BatchEvaluator {
public:
    LifetimeBarrierEvaluator(std::barrier<>& entered, std::barrier<>& release,
                             std::atomic<bool>& destroyed)
        : entered_(entered), release_(release), destroyed_(destroyed) {}
    ~LifetimeBarrierEvaluator() override { destroyed_.store(true); }

    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& candidates,
                                                 const EvaluationContext&) override {
        entered_.arrive_and_wait();
        release_.arrive_and_wait();
        std::vector<EvaluationResult> results;
        results.reserve(candidates.size());
        for ([[maybe_unused]] const auto& candidate : candidates) {
            auto result = EvaluationResult::success();
            result.objectives.push_back({"value", value_, true});
            results.push_back(std::move(result));
        }
        return results;
    }

private:
    std::barrier<>& entered_;
    std::barrier<>& release_;
    std::atomic<bool>& destroyed_;
    double value_ = 7.0;
};
}

TEST_CASE("serial adapter preserves order, isolates failures, and accepts empty batches") {
    SerialBatchEvaluator adapter{std::make_shared<IncrementingEvaluator>()};
    const EvaluationContext context{77, false};
    const auto empty = adapter.evaluate_batch(std::vector<CandidateVariables>{}, context);
    CHECK(empty.empty());

    const std::vector<CandidateVariables> candidates{cv(3.0), cv(-1.0), cv(5.0)};
    const auto results = adapter.evaluate_batch(candidates, context);
    REQUIRE(results.size() == 3);
    CHECK(results[0].objectives[0].value == 3.0);
    CHECK(results[1].status == EvaluationStatus::Failed);
    CHECK(results[1].diagnostics[0].code == "evaluation_exception");
    CHECK(results[2].objectives[0].value == 5.0);
}

TEST_CASE("cache key is stable and cached adapter reports hits") {
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    auto serial = std::make_shared<SerialBatchEvaluator>(std::make_shared<IncrementingEvaluator>());
    CachedBatchEvaluator adapter{serial, cache};
    EvaluationContext context{11, false};
    const std::vector<CandidateVariables> candidates{cv(2.0), cv(2.0)};
    const auto first = adapter.evaluate_batch(candidates, context);
    const auto second = adapter.evaluate_batch(std::vector<CandidateVariables>{cv(2.0)}, context);
    REQUIRE(first.size() == 2);
    CHECK(second[0].objectives[0].value == 2.0);
    CHECK(adapter.statistics().cache_hits == 1);
    CHECK(make_cache_key(CandidateVariables{{2.0}}, context) == make_cache_key(CandidateVariables{{2.0}}, context));
    CHECK(make_cache_key(CandidateVariables{{2.0}}, EvaluationContext{11, false}) !=
          make_cache_key(CandidateVariables{{2.0}}, EvaluationContext{12, false}));
    CHECK(make_cache_key(CandidateVariables{{2.0}}, EvaluationContext{11, false}) !=
          make_cache_key(CandidateVariables{{2.0}}, EvaluationContext{11, true}));
}

TEST_CASE("cache identity is deterministic and rejects an empty namespace") {
    const EvaluationCacheIdentity identity{"solver.cpu", "v1"};
    const EvaluationContext context{11, false};
    const auto variables = CandidateVariables{{2.0}};

    CHECK(make_cache_key(identity, variables, context) == make_cache_key(identity, variables, context));
    CHECK(make_cache_key(identity, variables, context) !=
          make_cache_key(EvaluationCacheIdentity{"solver.cpu", "v2"}, variables, context));
    CHECK(make_cache_key(identity, variables, context) !=
          make_cache_key(EvaluationCacheIdentity{"solver.gpu", "v1"}, variables, context));
    CHECK_THROWS_AS(EvaluationCacheIdentity("", "v1"), std::invalid_argument);
    CHECK_THROWS_AS(EvaluationCacheIdentity("solver.cpu", ""), std::invalid_argument);
}

TEST_CASE("cached evaluators sharing storage are isolated by namespace") {
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    auto first_delegate = std::make_shared<IdentifiedBatchEvaluator>(
        EvaluationCacheIdentity{"solver.first", "v1"}, "first", 1.0);
    auto second_delegate = std::make_shared<IdentifiedBatchEvaluator>(
        EvaluationCacheIdentity{"solver.second", "v1"}, "second", 2.0);
    CachedBatchEvaluator first{first_delegate, cache};
    CachedBatchEvaluator second{second_delegate, cache};
    const EvaluationContext context{22, false};

    const auto first_result = first.evaluate_batch({cv(3.0)}, context);
    const auto second_result = second.evaluate_batch({cv(3.0)}, context);

    REQUIRE(first_result.front().objectives.size() == 1);
    CHECK(first_result.front().objectives.front().id == "first");
    CHECK(first_result.front().objectives.front().value == 1.0);
    REQUIRE(second_result.front().objectives.size() == 1);
    CHECK(second_result.front().objectives.front().id == "second");
    CHECK(second_result.front().objectives.front().value == 2.0);
    CHECK(second.statistics().cache_hits == 0);

    const auto repeated_second_result = second.evaluate_batch({cv(3.0)}, context);
    CHECK(repeated_second_result.front().objectives.front().id == "second");
    CHECK(first_delegate->calls == 1);
    CHECK(second_delegate->calls == 1);
    CHECK(second.statistics().cache_hits == 1);
}

TEST_CASE("cached evaluators sharing a namespace are isolated by version") {
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    auto version_one_delegate = std::make_shared<IdentifiedBatchEvaluator>(
        EvaluationCacheIdentity{"solver.shared", "v1"}, "v1", 1.0);
    auto version_two_delegate = std::make_shared<IdentifiedBatchEvaluator>(
        EvaluationCacheIdentity{"solver.shared", "v2"}, "v2", 2.0);
    CachedBatchEvaluator version_one{version_one_delegate, cache};
    CachedBatchEvaluator version_two{version_two_delegate, cache};
    const EvaluationContext context{23, false};

    version_one.evaluate_batch({cv(4.0)}, context);
    const auto version_two_result = version_two.evaluate_batch({cv(4.0)}, context);

    REQUIRE(version_two_result.front().objectives.size() == 1);
    CHECK(version_two_result.front().objectives.front().id == "v2");
    CHECK(version_two_result.front().objectives.front().value == 2.0);
    CHECK(version_one_delegate->calls == 1);
    CHECK(version_two_delegate->calls == 1);
    CHECK(version_two.statistics().cache_hits == 0);
}

TEST_CASE("cached evaluation snapshots keep callback identity paired with its result") {
    std::barrier switch_barrier{2};
    auto delegate = std::make_shared<GenerationSwitchBatchEvaluator>(switch_barrier);
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    CachedBatchEvaluator cached{delegate, cache};
    const EvaluationContext context{24, false};

    std::vector<EvaluationResult> first_results;
    std::thread worker([&] { first_results = cached.evaluate_batch({cv(7.0)}, context); });
    switch_barrier.arrive_and_wait();
    delegate->replace_callback_generation();
    switch_barrier.arrive_and_wait();
    worker.join();

    REQUIRE(first_results.size() == 1);
    CHECK(first_results.front().objectives.front().value == 1.0);

    auto old_identity = std::make_shared<IdentifiedBatchEvaluator>(
        EvaluationCacheIdentity{"generation-switch-1", "v1"}, "generation", 1.0);
    CachedBatchEvaluator old_generation{old_identity, cache};
    const auto old_result = old_generation.evaluate_batch({cv(7.0)}, context);
    REQUIRE(old_result.size() == 1);
    CHECK(old_result.front().objectives.front().value == 1.0);
    CHECK(old_generation.statistics().cache_hits == 1);
    CHECK(old_identity->calls == 0);
}

TEST_CASE("unmanaged evaluators reject escaping snapshots") {
    IdentifiedBatchEvaluator evaluator{
        EvaluationCacheIdentity{"lifetime", "v1"}, "value", 4.0};

    CHECK_THROWS_WITH_AS(evaluator.evaluation_snapshot(),
                         "BatchEvaluator::evaluation_snapshot requires shared ownership",
                         std::logic_error);
}

TEST_CASE("shared evaluator snapshots retain the owner after external reset") {
    auto evaluator = std::make_shared<IdentifiedBatchEvaluator>(
        EvaluationCacheIdentity{"lifetime", "v1"}, "value", 4.0);
    auto snapshot = evaluator->evaluation_snapshot();
    evaluator.reset();

    const auto results = snapshot.evaluate({cv(1.0)}, {});
    REQUIRE(results.size() == 1);
    CHECK(results.front().status == EvaluationStatus::Success);
    CHECK(results.front().objectives.front().value == 4.0);
}

TEST_CASE("active shared snapshots delay final destruction without deadlock") {
    std::barrier entered{2};
    std::barrier release{2};
    std::atomic<bool> destroyed{false};
    auto evaluator = std::make_shared<LifetimeBarrierEvaluator>(entered, release, destroyed);
    auto snapshot = evaluator->evaluation_snapshot();
    std::vector<EvaluationResult> results;
    std::thread worker([snapshot = std::move(snapshot), &results]() mutable {
        results = snapshot.evaluate({cv(1.0)}, {});
    });

    entered.arrive_and_wait();
    evaluator.reset();
    CHECK_FALSE(destroyed.load());
    release.arrive_and_wait();
    worker.join();

    REQUIRE(results.size() == 1);
    CHECK(results.front().objectives.front().value == 7.0);
    CHECK(destroyed.load());
}

TEST_CASE("shared cache and statistics wrapper snapshots retain the whole chain") {
    auto leaf = std::make_shared<IdentifiedBatchEvaluator>(
        EvaluationCacheIdentity{"wrapper-lifetime", "v1"}, "value", 6.0);
    auto statistics = std::make_shared<StatisticsBatchEvaluator>(leaf);
    auto cached = std::make_shared<CachedBatchEvaluator>(
        statistics, std::make_shared<InMemoryEvaluationCache>());
    auto snapshot = cached->evaluation_snapshot();
    cached.reset();
    statistics.reset();
    leaf.reset();

    const auto results = snapshot.evaluate({cv(1.0)}, {});
    REQUIRE(results.size() == 1);
    CHECK(results.front().status == EvaluationStatus::Success);
    CHECK(results.front().objectives.front().value == 6.0);
}

TEST_CASE("cached batch evaluation isolates delegate exceptions and malformed output") {
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    CachedBatchEvaluator throwing{std::make_shared<ThrowingBatchEvaluator>(), cache};
    const auto results = throwing.evaluate_batch({cv(1.0), cv(-1.0), cv(2.0)}, EvaluationContext{4, false});
    REQUIRE(results.size() == 3);
    CHECK(results[0].status == EvaluationStatus::Success);
    CHECK(results[1].status == EvaluationStatus::Failed);
    CHECK(results[1].diagnostics[0].code == "evaluation_batch_output");
    CHECK(results[2].status == EvaluationStatus::Failed);
    CHECK(results[2].diagnostics[0].code == "evaluation_batch_output");

    CachedBatchEvaluator empty{std::make_shared<EmptyBatchEvaluator>(), std::make_shared<InMemoryEvaluationCache>()};
    const auto malformed = empty.evaluate_batch({cv(3.0)}, EvaluationContext{4, false});
    REQUIRE(malformed.size() == 1);
    CHECK(malformed[0].status == EvaluationStatus::Failed);
    CHECK(malformed[0].diagnostics[0].code == "evaluation_batch_output");
}

TEST_CASE("cached batch retries thrown batches per candidate") {
    auto delegate = std::make_shared<CandidateThrowingBatchEvaluator>();
    CachedBatchEvaluator adapter{delegate, std::make_shared<InMemoryEvaluationCache>()};
    const EvaluationContext context{5, false};

    const auto results = adapter.evaluate_batch({cv(1.0), cv(-1.0), cv(2.0)}, context);
    REQUIRE(results.size() == 3);
    CHECK(results[0].status == EvaluationStatus::Success);
    CHECK(results[0].objectives.front().value == 1.0);
    CHECK(results[1].status == EvaluationStatus::Failed);
    CHECK(results[1].diagnostics.front().code == "evaluation_exception");
    CHECK(results[2].status == EvaluationStatus::Success);
    CHECK(results[2].objectives.front().value == 2.0);
    CHECK(delegate->calls == 4);
    CHECK(adapter.statistics().evaluations == 3);
    CHECK(adapter.statistics().successful_evaluations == 2);
    CHECK(adapter.statistics().failed_evaluations == 1);

    const auto cached = adapter.evaluate_batch({cv(1.0), cv(-1.0), cv(2.0)}, context);
    REQUIRE(cached.size() == 3);
    CHECK(cached[0].status == EvaluationStatus::Success);
    CHECK(cached[1].status == EvaluationStatus::Failed);
    CHECK(cached[2].status == EvaluationStatus::Success);
    CHECK(delegate->calls == 4);
    CHECK(adapter.statistics().evaluations == 3);
    CHECK(adapter.statistics().cache_hits == 3);
}

TEST_CASE("nested cached statistics evaluator retries thrown batches per candidate") {
    auto delegate = std::make_shared<CandidateThrowingBatchEvaluator>();
    auto tracked = std::make_shared<StatisticsBatchEvaluator>(delegate);
    CachedBatchEvaluator adapter{tracked, std::make_shared<InMemoryEvaluationCache>()};
    const auto results = adapter.evaluate_batch({cv(1.0), cv(-1.0), cv(2.0)}, EvaluationContext{5, false});
    REQUIRE(results.size() == 3);
    CHECK(results[0].status == EvaluationStatus::Success);
    CHECK(results[1].status == EvaluationStatus::Failed);
    CHECK(results[2].status == EvaluationStatus::Success);
    CHECK(delegate->calls == 4);
}

TEST_CASE("cached evaluator batches misses once and restores input order") {
    auto delegate = std::make_shared<RecordingBatchEvaluator>();
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    const EvaluationContext context{8, false};
    auto cached = EvaluationResult::success();
    cached.objectives.push_back({"value", 10.0, true});
    cache->put(make_cache_key(cv(10.0), context), cached);

    CachedBatchEvaluator adapter{delegate, cache};
    const auto results = adapter.evaluate_batch({cv(10.0), cv(2.0), cv(3.0), cv(10.0)}, context);
    REQUIRE(results.size() == 4);
    CHECK(delegate->calls == 1);
    CHECK(delegate->batches == std::vector<double>{2.0, 3.0});
    CHECK(results[0].objectives.front().value == 10.0);
    CHECK(results[1].objectives.front().value == 2.0);
    CHECK(results[2].objectives.front().value == 3.0);
    CHECK(results[3].objectives.front().value == 10.0);
}

TEST_CASE("unevaluated delegate results are normalized as failures") {
    CachedBatchEvaluator adapter{std::make_shared<UnevaluatedBatchEvaluator>(),
                                 std::make_shared<InMemoryEvaluationCache>()};
    const auto results = adapter.evaluate_batch({cv(1.0)}, EvaluationContext{12, false});
    REQUIRE(results.size() == 1);
    CHECK(results.front().status == EvaluationStatus::Failed);
    CHECK(results.front().diagnostics.front().code == "evaluation_unevaluated");
    CHECK(adapter.statistics().failed_evaluations == 1);
}

TEST_CASE("short delegate batches preserve returned results and fail only missing candidates") {
    StatisticsBatchEvaluator adapter{std::make_shared<ShortBatchEvaluator>()};
    const auto results = adapter.evaluate_batch({cv(4.0), cv(5.0)}, EvaluationContext{13, false});
    REQUIRE(results.size() == 2);
    CHECK(results[0].status == EvaluationStatus::Success);
    CHECK(results[0].objectives.front().value == 4.0);
    CHECK(results[1].status == EvaluationStatus::Failed);
    CHECK(results[1].diagnostics.front().code == "evaluation_batch_output");
    CHECK(adapter.statistics().successful_evaluations == 1);
    CHECK(adapter.statistics().failed_evaluations == 1);
}

TEST_CASE("cached statistics count misses, statuses, cache hits, fallback, and duration") {
    auto cache = std::make_shared<InMemoryEvaluationCache>();
    CachedBatchEvaluator adapter{std::make_shared<SerialBatchEvaluator>(std::make_shared<IncrementingEvaluator>()), cache};
    const EvaluationContext context{9, true};
    const auto first = adapter.evaluate_batch({cv(1.0), cv(-1.0), cv(1.0)}, context);
    REQUIRE(first.size() == 3);
    CHECK(adapter.statistics().evaluations == 2);
    CHECK(adapter.statistics().successful_evaluations == 1);
    CHECK(adapter.statistics().failed_evaluations == 1);
    CHECK(adapter.statistics().cache_hits == 0);
    CHECK(adapter.statistics().fallbacks == 1);
    CHECK(adapter.statistics().seed == 9);
    CHECK(adapter.statistics().elapsed_seconds >= 0.0);

    adapter.evaluate_batch({cv(1.0), cv(-1.0)}, context);
    CHECK(adapter.statistics().evaluations == 2);
    CHECK(adapter.statistics().cache_hits == 2);
    CHECK(adapter.statistics().successful_evaluations == 1);
    CHECK(adapter.statistics().failed_evaluations == 1);
    CHECK(adapter.statistics().fallbacks == 2);
    CHECK(adapter.statistics().elapsed_seconds >= 0.0);
}

TEST_CASE("statistics adapter normalizes malformed batch output") {
    StatisticsBatchEvaluator adapter{std::make_shared<EmptyBatchEvaluator>()};
    const auto results = adapter.evaluate_batch({cv(1.0), cv(2.0)}, EvaluationContext{6, false});
    REQUIRE(results.size() == 2);
    CHECK(results[0].status == EvaluationStatus::Failed);
    CHECK(results[1].status == EvaluationStatus::Failed);
    CHECK(adapter.statistics().evaluations == 2);
    CHECK(adapter.statistics().successful_evaluations == 0);
    CHECK(adapter.statistics().failed_evaluations == 2);
    CHECK(adapter.statistics().elapsed_seconds >= 0.0);
}

TEST_CASE("statistics tracks evaluations, failures, fallback, seed, and duration") {
    auto serial = std::make_shared<SerialBatchEvaluator>(std::make_shared<IncrementingEvaluator>());
    StatisticsBatchEvaluator adapter{serial};
    EvaluationContext context{1234, true};
    const auto results = adapter.evaluate_batch(std::vector<CandidateVariables>{cv(1.0), cv(-1.0)}, context);
    REQUIRE(results.size() == 2);
    CHECK(adapter.statistics().evaluations == 2);
    CHECK(adapter.statistics().failed_evaluations == 1);
    CHECK(adapter.statistics().fallbacks == 1);
    CHECK(adapter.statistics().seed == 1234);
    CHECK(adapter.statistics().elapsed_seconds >= 0.0);
}

TEST_CASE("statistics snapshots retain nested actual fallbacks across cache hits") {
    auto source = std::make_shared<FallbackReportingEvaluator>();
    auto tracked = std::make_shared<StatisticsBatchEvaluator>(source);
    CachedBatchEvaluator cached(tracked, std::make_shared<InMemoryEvaluationCache>());
    const EvaluationContext context{42, false};

    cached.evaluate_batch({cv(1.0)}, context);
    REQUIRE(cached.statistics_snapshot());
    CHECK(cached.statistics_snapshot()->fallbacks == 1);
    CHECK(cached.statistics_snapshot()->cache_hits == 0);

    cached.evaluate_batch({cv(1.0)}, context);
    REQUIRE(cached.statistics_snapshot());
    CHECK(cached.statistics_snapshot()->fallbacks == 1);
    CHECK(cached.statistics_snapshot()->cache_hits == 1);

    cached.evaluate_batch({cv(2.0)}, context);
    REQUIRE(cached.statistics_snapshot());
    CHECK(cached.statistics_snapshot()->fallbacks == 2);
}

TEST_CASE("run statistics collector is context-owned and returns value snapshots") {
    auto collector = std::make_shared<EvaluationStatisticsCollector>();
    EvaluationContext context{17, false, collector};
    collector->add_cache_hits(2);
    collector->add_fallbacks(1);

    const auto first = collector->snapshot();
    collector->add_cache_hits(3);
    const auto second = collector->snapshot();

    CHECK(context.statistics == collector);
    CHECK(first.cache_hits == 2);
    CHECK(first.fallbacks == 1);
    CHECK(second.cache_hits == 5);
    CHECK(first.cache_hits == 2);
}

TEST_CASE("aggregate statistics access returns an independent value snapshot") {
    StatisticsBatchEvaluator adapter{std::make_shared<SerialBatchEvaluator>(
        std::make_shared<IncrementingEvaluator>())};
    adapter.evaluate_batch({cv(1.0)}, EvaluationContext{});
    const auto snapshot = adapter.statistics();
    adapter.evaluate_batch({cv(2.0)}, EvaluationContext{});

    CHECK(snapshot.evaluations == 1);
    CHECK(adapter.statistics().evaluations == 2);
}

TEST_CASE("nested statistics wrappers forward a legacy fallback marker once per call") {
    auto leaf = std::make_shared<SerialBatchEvaluator>(std::make_shared<IncrementingEvaluator>());
    auto inner = std::make_shared<StatisticsBatchEvaluator>(leaf);
    StatisticsBatchEvaluator outer(inner);
    auto collector = std::make_shared<EvaluationStatisticsCollector>();
    EvaluationContext context{19, true, collector};

    outer.evaluate_batch({cv(1.0)}, context);
    CHECK(collector->snapshot().fallbacks == 1);
    outer.evaluate_batch({cv(2.0)}, context);
    CHECK(collector->snapshot().fallbacks == 2);
}

TEST_CASE("collector add preserves the run-owned seed") {
    EvaluationStatisticsCollector collector;
    collector.set_seed(101);
    EvaluationStatistics delta;
    delta.seed = 202;
    delta.cache_hits = 3;

    collector.add(delta);

    const auto snapshot = collector.snapshot();
    CHECK(snapshot.seed == 101);
    CHECK(snapshot.cache_hits == 3);
}

TEST_CASE("same-context concurrent wrapper calls each forward one fallback marker") {
    class BlockingEvaluator final : public BatchEvaluator {
    public:
        explicit BlockingEvaluator(std::barrier<>& barrier) : barrier_(barrier) {}

        std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>& values,
                                                     const EvaluationContext&) override {
            barrier_.arrive_and_wait();
            std::vector<EvaluationResult> results;
            for (const auto& value : values) {
                auto result = EvaluationResult::success();
                result.objectives.push_back({"value", value.values.front(), true});
                results.push_back(std::move(result));
            }
            return results;
        }

    private:
        std::barrier<>& barrier_;
    };

    std::barrier barrier(2);
    auto leaf = std::make_shared<BlockingEvaluator>(barrier);
    auto inner = std::make_shared<StatisticsBatchEvaluator>(leaf);
    StatisticsBatchEvaluator outer(inner);
    auto collector = std::make_shared<EvaluationStatisticsCollector>();
    EvaluationContext context{23, true, collector};

    std::thread first([&] { outer.evaluate_batch({cv(1.0)}, context); });
    std::thread second([&] { outer.evaluate_batch({cv(2.0)}, context); });
    first.join();
    second.join();

    CHECK(collector->snapshot().fallbacks == 2);
}
