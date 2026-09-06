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

class EmptyBatchEvaluator final : public BatchEvaluator {
public:
    std::vector<EvaluationResult> evaluate_batch(const std::vector<CandidateVariables>&,
                                                 const EvaluationContext&) override {
        return {};
    }
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
    CHECK(adapter.statistics().cache_hits == 2);
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
    CHECK(results[1].diagnostics[0].code == "evaluation_exception");
    CHECK(results[2].status == EvaluationStatus::Success);

    CachedBatchEvaluator empty{std::make_shared<EmptyBatchEvaluator>(), std::make_shared<InMemoryEvaluationCache>()};
    const auto malformed = empty.evaluate_batch({cv(3.0)}, EvaluationContext{4, false});
    REQUIRE(malformed.size() == 1);
    CHECK(malformed[0].status == EvaluationStatus::Failed);
    CHECK(malformed[0].diagnostics[0].code == "evaluation_batch_output");
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
    CHECK(adapter.statistics().cache_hits == 1);
    CHECK(adapter.statistics().fallbacks == 1);
    CHECK(adapter.statistics().seed == 9);
    CHECK(adapter.statistics().elapsed_seconds >= 0.0);

    adapter.evaluate_batch({cv(1.0), cv(-1.0)}, context);
    CHECK(adapter.statistics().evaluations == 2);
    CHECK(adapter.statistics().cache_hits == 3);
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
