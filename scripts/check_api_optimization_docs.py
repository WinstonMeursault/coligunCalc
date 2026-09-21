#!/usr/bin/env python3
"""Check the bilingual optimization API contract and real CMake presets."""

from pathlib import Path
import argparse
import json
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
DOCS = {
    "en": (ROOT / "docs/API.md", "### Optimization", "## Physics Constants"),
    "zh": (ROOT / "docs/API_cn.md", "### 优化", "## 物理常量"),
}

# Each entry is a semantic contract, not merely a bag of translated words.
# Marker order is checked for the routing/fallback clauses so reversed policy
# documentation cannot pass accidentally.
SEMANTIC_MARKERS = {
    "problem_contract": {
        "en": ("ProblemSpec", "immutable run contract", "VariableSpec", "RepairPolicy"),
        "zh": ("ProblemSpec", "不可变的运行契约", "VariableSpec", "RepairPolicy"),
    },
    "routing": {
        "en": ("SelectionStrategy::Auto", "one declared objective", "SingleObjective",
               "two or more objectives", "NSGA-II", "NSGA2"),
        "zh": ("SelectionStrategy::Auto", "一个目标", "标量遗传算法", "SingleObjective",
               "两个或更多目标", "NSGA-II", "NSGA2"),
    },
    "pareto_selection": {
        "en": ("Pareto front only", "select_representative", "MaxObjective",
               "LexicographicObjectives", "std::invalid_argument", "std::out_of_range"),
        "zh": ("只返回 Pareto front", "select_representative", "MaxObjective",
               "LexicographicObjectives", "std::invalid_argument", "std::out_of_range"),
    },
    "statistics_ownership": {
        "en": ("EvaluationCacheIdentity", "run-local statistics snapshot",
               "statistics_snapshot", "gpu_fallbacks", "legacy overload",
               "that value wins", "not added twice"),
        "zh": ("EvaluationCacheIdentity", "新的运行期统计快照",
               "statistics_snapshot", "gpu_fallbacks", "旧版重载", "使用该值",
               "不会重复相加"),
    },
    "snapshot_lifetime": {
        "en": ("BatchEvaluationSnapshot", "evaluation_snapshot()", "final lifecycle boundary",
               "std::shared_ptr", "requires shared ownership", "std::logic_error",
               "protected", "make_evaluation_snapshot", "asynchronously"),
        "zh": ("BatchEvaluationSnapshot", "evaluation_snapshot()", "生命周期边界",
               "std::shared_ptr", "需要共享所有权", "std::logic_error",
               "受保护", "make_evaluation_snapshot", "异步调用"),
    },
    "cuda_policy": {
        "en": ("CudaBatchEvaluator", "CudaFallbackPolicy::Strict", "PerCandidateCpu",
               "WholeBatchCpu", "protocol errors", "PeakCurrent", "ExecutionReport::gpu_executed"),
        "zh": ("CudaBatchEvaluator", "CudaFallbackPolicy::Strict", "PerCandidateCpu",
               "WholeBatchCpu", "协议错误", "PeakCurrent", "ExecutionReport::gpu_executed"),
    },
    "installation": {
        "en": ("CPU-only installations", "coilgun::coilgun", "coilgun::coilgun_cuda",
               "not install the CUDA batch adapter"),
        "zh": ("仅 CPU 安装", "coilgun::coilgun", "coilgun::coilgun_cuda", "不安装 CUDA 批量适配器"),
    },
}

COMMON_TOKENS = {
    "ProblemSpec", "VariableSpec", "ObjectiveDefinition", "ConstraintDefinition", "RepairPolicy",
    "SelectionStrategy::Auto", "SingleObjective", "NSGA2", "select_representative",
    "MaxObjective", "MinConstraintViolationMargin", "IdealPointDistance", "WeightedScore",
    "LexicographicObjectives", "CallbackSelector", "EvaluationCacheIdentity", "make_cache_key",
    "statistics_snapshot", "BatchEvaluationSnapshot", "evaluation_snapshot",
    "CoilgunOptimizationProblem", "CudaBatchEvaluator",
    "CudaFallbackPolicy::Strict", "PerCandidateCpu", "WholeBatchCpu", "SimBatch<EulerStepper>",
    "OptimizationLevel::Full", "PeakCurrent", "ExecutionReport::gpu_executed",
    "BackendMode::Fallback", "coilgun::coilgun", "coilgun::coilgun_cuda",
}

REQUIRED_PRESETS = {"cpu-debug", "cpu-release", "cuda-debug", "cuda-release"}


def section(path: Path, start: str, end: str) -> str:
    text = path.read_text(encoding="utf-8")
    try:
        return text.split(start, 1)[1].split(end, 1)[0]
    except IndexError as exc:
        raise SystemExit(f"missing section marker in {path}") from exc


def actual_presets() -> tuple[set[str], set[str], set[str]]:
    data = json.loads((ROOT / "CMakePresets.json").read_text(encoding="utf-8"))
    configure = {item["name"] for item in data.get("configurePresets", [])}
    build = {item["name"] for item in data.get("buildPresets", [])}
    test = {item["name"] for item in data.get("testPresets", [])}
    return configure, build, test


def check(sections: dict[str, str], full_docs: dict[str, str]) -> list[str]:
    failures: list[str] = []
    for language, body in sections.items():
        for token in sorted(COMMON_TOKENS):
            if token not in body:
                failures.append(f"{language}: missing API token {token}")
        for name, markers_by_language in SEMANTIC_MARKERS.items():
            markers = markers_by_language[language]
            missing = [marker for marker in markers if marker not in body]
            if missing:
                failures.append(f"{language}: {name} missing {missing}")

        fences = re.findall(r"```(?:cpp)?\n(.*?)```", body, re.DOTALL)
        if len(fences) != 1:
            failures.append(f"{language}: expected one optimization API example, found {len(fences)}")
        elif not all(token in fences[0] for token in ("ProblemSpec", "GeneticOptimizer", "pareto_front")):
            failures.append(f"{language}: optimization example is incomplete")
        if len(body) < 2500:
            failures.append(f"{language}: optimization section is unexpectedly short")

    # Ordered markers make the policy meaning executable, not just present.
    for language, body in sections.items():
        start = body.find("SelectionStrategy::Auto")
        end_marker = "Multi-objective runs" if language == "en" else "多目标运行"
        end = body.find(end_marker, start)
        routing_body = body[start:end] if start >= 0 and end >= 0 else ""
        expected = ("one declared objective", "SingleObjective", "two or more objectives", "NSGA-II")
        if language == "zh":
            expected = ("一个目标", "标量遗传算法", "SingleObjective", "两个或更多目标", "NSGA-II")
        positions = [routing_body.find(marker) for marker in expected]
        if any(position < 0 for position in positions) or positions != sorted(positions):
            failures.append(f"{language}: Auto routing order is not one-objective → scalar GA → multi-objective → NSGA-II")

    for language, body in sections.items():
        ordered = ("CudaFallbackPolicy::Strict", "PerCandidateCpu", "WholeBatchCpu")
        if language == "zh":
            ordered = ("CudaFallbackPolicy::Strict", "PerCandidateCpu", "WholeBatchCpu")
        positions = [body.rfind(marker) for marker in ordered]
        if any(position < 0 for position in positions) or positions != sorted(positions):
            failures.append(f"{language}: fallback policy order is not Strict → PerCandidateCpu → WholeBatchCpu")

    examples = {
        language: re.findall(r"```(?:cpp)?\n(.*?)```", body, re.DOTALL)[0]
        for language, body in sections.items()
        if re.findall(r"```(?:cpp)?\n(.*?)```", body, re.DOTALL)
    }
    if len(examples) == 2 and examples["en"] != examples["zh"]:
        failures.append("English and Chinese optimization examples differ")

    configure, build, test = actual_presets()
    for kind, available in (("configure", configure), ("build", build), ("test", test)):
        missing = REQUIRED_PRESETS - available
        if missing:
            failures.append(f"CMakePresets.json {kind} presets missing {sorted(missing)}")
    documented = {
        language: set(re.findall(r"`(cpu-(?:debug|release)|cuda-(?:debug|release))`", text))
        for language, text in full_docs.items()
    }
    for language, names in documented.items():
        if REQUIRED_PRESETS - names:
            failures.append(f"{language}: docs missing required presets {sorted(REQUIRED_PRESETS - names)}")
        if names - configure - build - test:
            failures.append(f"{language}: docs mention unknown presets {sorted(names - configure - build - test)}")
    if documented["en"] != documented["zh"]:
        failures.append(f"English and Chinese preset lists differ: {documented}")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true", help="run a negative omission check")
    args = parser.parse_args()
    sections = {language: section(path, start, end)
                for language, (path, start, end) in DOCS.items()}
    full_docs = {language: path.read_text(encoding="utf-8")
                 for language, (path, _, _) in DOCS.items()}
    failures = check(sections, full_docs)
    if args.self_test:
        broken = dict(sections)
        broken["zh"] = broken["zh"].replace("NSGA2", "REMOVED", 1)
        negative = check(broken, full_docs)
        if not any("zh: missing API token NSGA2" in failure for failure in negative):
            print("negative self-test: FAIL")
            return 1
        print("negative self-test: PASS")
    if failures:
        print("API optimization bilingual check: FAIL")
        print("\n".join(f"- {failure}" for failure in failures))
        return 1
    configure, build, test = actual_presets()
    print("API optimization bilingual check: PASS")
    print(f"semantic_contracts={len(SEMANTIC_MARKERS)} common_tokens={len(COMMON_TOKENS)}")
    print("example_blocks=1 per language")
    print(f"actual_presets=configure:{len(configure)} build:{len(build)} test:{len(test)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
