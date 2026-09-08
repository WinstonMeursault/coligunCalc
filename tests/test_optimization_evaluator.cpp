#include <doctest/doctest.h>

#include "coilgun/optimization/cache.hpp"
#include "coilgun/optimization/evaluator.hpp"
#include "coilgun/optimization/statistics.hpp"

#include <stdexcept>

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
