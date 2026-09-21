#!/usr/bin/env python3
"""Validate machine-readable B-T4 GPU workflow evidence."""

import argparse
import json
import math
import re
import sys
import statistics
from typing import Any


def _finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _close(actual: Any, expected: float) -> bool:
    return _finite(actual) and math.isclose(float(actual), expected, rel_tol=1e-9, abs_tol=1e-9)


def validate_document(document: dict[str, Any]) -> list[str]:
    errors: list[str] = []

    def required(mapping: dict[str, Any], key: str, path: str) -> Any:
        if key not in mapping:
            errors.append(f"missing required key {path}.{key}")
            return None
        return mapping[key]

    if not isinstance(document, dict):
        return ["document must be an object"]
    if required(document, "schema_version", "document") != 1:
        errors.append("schema_version must equal 1")
    for key in ("source_revision", "worktree_state", "toolchain", "workload_provenance", "hardware", "workload", "batches", "optimization", "decision"):
        required(document, key, "document")
    source_revision = document.get("source_revision")
    if not isinstance(source_revision, str) or not re.fullmatch(r"[0-9a-fA-F]{7,40}", source_revision):
        errors.append("source_revision must be a non-empty 7-40 character hexadecimal revision")
    if document.get("worktree_state") != "clean":
        errors.append("worktree_state must be clean")

    toolchain = document.get("toolchain")
    if not isinstance(toolchain, dict):
        errors.append("toolchain must be an object")
    else:
        for key in ("compiler", "cmake", "build_type"):
            value = required(toolchain, key, "toolchain")
            if not isinstance(value, str) or not value.strip():
                errors.append(f"toolchain.{key} must be a non-empty string")

    workload_provenance = document.get("workload_provenance")
    if not isinstance(workload_provenance, dict):
        errors.append("workload_provenance must be an object")
    else:
        for key in ("name", "configuration", "command"):
            value = required(workload_provenance, key, "workload_provenance")
            if not isinstance(value, str) or not value.strip():
                errors.append(f"workload_provenance.{key} must be a non-empty string")

    hardware = document.get("hardware")
    if isinstance(hardware, dict):
        for key in ("gpu_name", "driver", "compute_capability", "cuda_toolkit"):
            required(hardware, key, "hardware")
    else:
        errors.append("hardware must be an object")

    workload = document.get("workload")
    repetitions = None
    batch_sizes = None
    if isinstance(workload, dict):
        for key in ("batch_sizes", "warmups", "repetitions", "precision", "seed"):
            required(workload, key, "workload")
        repetitions = workload.get("repetitions")
        batch_sizes = workload.get("batch_sizes")
        if workload.get("precision") != "full":
            errors.append("workload.precision must be full")
        if not isinstance(repetitions, int) or repetitions < 5:
            errors.append("workload.repetitions must be at least five")
        if not isinstance(workload.get("warmups"), int) or workload.get("warmups") < 1:
            errors.append("workload.warmups must be positive")
        if batch_sizes != [1, 8, 32, 128]:
            errors.append("workload.batch_sizes must be [1, 8, 32, 128]")
    else:
        errors.append("workload must be an object")

    batches = document.get("batches")
    if not isinstance(batches, list):
        errors.append("batches must be an array")
        batches = []
    elif batch_sizes is not None and [row.get("batch_size") for row in batches if isinstance(row, dict)] != batch_sizes:
        errors.append("batch count/order does not match workload.batch_sizes")

    timing_keys = ("cpu_samples_ms", "gpu_samples_ms")
    summary_keys = ("cpu_median_ms", "gpu_median_ms", "cpu_per_candidate_ms", "gpu_per_candidate_ms", "cpu_candidates_per_second", "gpu_candidates_per_second", "speedup")
    for index, row in enumerate(batches):
        path = f"batches[{index}]"
        if not isinstance(row, dict):
            errors.append(f"{path} must be an object")
            continue
        for key in ("batch_size", "cpu_rows", "gpu_rows", "representative_rows", "order_ok", "numerical_ok", "gpu_executed", "backend", "requested_backend", "solver", "precision", "requested_precision", "fallback_events", "gpu_report_time_ms", "transfer_time_ms", *timing_keys, *summary_keys):
            required(row, key, path)
        size = row.get("batch_size")
        if not isinstance(size, int) or size <= 0:
            errors.append(f"{path}.batch_size must be positive")
        for count_key in ("cpu_rows", "gpu_rows"):
            if row.get(count_key) != size:
                errors.append(f"{path} row count must equal batch_size")
        if row.get("order_ok") is not True:
            errors.append(f"{path}.order_ok must be true")
        if row.get("numerical_ok") is not True:
            errors.append(f"{path}.numerical_ok must be true")
        if row.get("gpu_executed") is not True:
            errors.append(f"{path}.gpu_executed must be true")
        if row.get("backend") == "fallback" or not isinstance(row.get("backend"), str) or not row.get("backend"):
            errors.append(f"{path}.backend must be a non-fallback backend")
        if row.get("requested_backend") == "fallback" or not isinstance(row.get("requested_backend"), str) or not row.get("requested_backend"):
            errors.append(f"{path}.requested_backend must be a non-fallback backend")
        if not isinstance(row.get("solver"), str) or not row.get("solver") or row.get("solver") == "unknown":
            errors.append(f"{path}.solver must identify a resolved solver")
        if row.get("precision") != "full":
            errors.append(f"{path}.precision must be full")
        if row.get("requested_precision") != "full":
            errors.append(f"{path}.requested_precision must be full")
        if row.get("fallback_events") != 0:
            errors.append(f"{path}.fallback_events must be zero")
        for key in ("gpu_report_time_ms", "transfer_time_ms"):
            if not _finite(row.get(key)) or row[key] < 0:
                errors.append(f"{path}.{key} must be finite and non-negative")
        if not isinstance(row.get("representative_rows"), int) or not 0 < row["representative_rows"] <= size:
            errors.append(f"{path}.representative_rows must be a positive count within batch_size")
        for key in timing_keys:
            samples = row.get(key)
            if not isinstance(samples, list):
                errors.append(f"{path}.{key} must be an array")
                continue
            if isinstance(repetitions, int) and len(samples) != repetitions:
                errors.append(f"{path}.{key} sample count does not match repetitions")
            for sample in samples:
                if not _finite(sample):
                    errors.append(f"{path}.{key} contains non-finite time")
                elif sample < 0:
                    errors.append(f"{path}.{key} contains non-negative time violation")
                elif sample == 0:
                    errors.append(f"{path}.{key} samples must be greater than zero")
        for key in summary_keys:
            value = row.get(key)
            if not _finite(value) or value < 0:
                errors.append(f"{path}.{key} must be finite and non-negative")
        cpu_samples = row.get("cpu_samples_ms")
        gpu_samples = row.get("gpu_samples_ms")
        if (isinstance(cpu_samples, list) and cpu_samples and all(_finite(sample) and sample > 0 for sample in cpu_samples)
                and isinstance(gpu_samples, list) and gpu_samples and all(_finite(sample) and sample > 0 for sample in gpu_samples)):
            derived = {
                "cpu_median_ms": statistics.median(cpu_samples),
                "gpu_median_ms": statistics.median(gpu_samples),
            }
            derived["cpu_per_candidate_ms"] = derived["cpu_median_ms"] / size
            derived["gpu_per_candidate_ms"] = derived["gpu_median_ms"] / size
            derived["cpu_candidates_per_second"] = 1000.0 * size / derived["cpu_median_ms"]
            derived["gpu_candidates_per_second"] = 1000.0 * size / derived["gpu_median_ms"]
            derived["speedup"] = derived["cpu_median_ms"] / derived["gpu_median_ms"]
            for key, expected in derived.items():
                if not _close(row.get(key), expected):
                    errors.append(f"{path}.{key} derived value mismatch")

    optimization = document.get("optimization")
    if not isinstance(optimization, dict):
        errors.append("optimization must be an object")
    else:
        if optimization.get("same_seed") is not True:
            errors.append("optimization.same_seed must be true")
        for key in ("population_size", "max_generations"):
            if not isinstance(optimization.get(key), int) or optimization[key] <= 0:
                errors.append(f"optimization.{key} must be positive")
        for name in ("cpu", "gpu"):
            run = optimization.get(name)
            if not isinstance(run, dict):
                errors.append(f"optimization.{name} must be an object")
                continue
            for key in ("termination", "best_voltage", "best_objective", "reference_objective", "reference_error", "reference_tolerance", "statistics"):
                required(run, key, f"optimization.{name}")
            for key in ("feasible", "reference_valid"):
                if run.get(key) is not True:
                    errors.append(f"optimization.{name}.{key} must be true")
            for key in ("best_voltage", "best_objective", "reference_objective", "reference_error", "reference_tolerance"):
                if not _finite(run.get(key)) or run[key] < 0:
                    errors.append(f"optimization.{name}.{key} must be finite and non-negative")
            if all(_finite(run.get(key)) for key in ("best_objective", "reference_objective", "reference_error", "reference_tolerance")):
                expected_tolerance = 5e-8 + 1e-6 * abs(run["reference_objective"])
                expected_error = abs(run["best_objective"] - run["reference_objective"])
                expected_valid = expected_error <= expected_tolerance
                if not _close(run.get("reference_tolerance"), expected_tolerance):
                    errors.append(f"optimization.{name}.reference_tolerance derived value mismatch")
                if not _close(run.get("reference_error"), expected_error):
                    errors.append(f"optimization.{name}.reference_error derived value mismatch")
                if run.get("reference_valid") is not expected_valid:
                    errors.append(f"optimization.{name}.reference_valid does not match reference error/tolerance")
            stats = run.get("statistics")
            stat_keys = ("seed", "evaluations", "successful_evaluations", "failed_evaluations", "cache_hits", "gpu_fallbacks", "gpu_requested_evaluations", "gpu_executed_evaluations", "gpu_successful_evaluations", "gpu_failed_evaluations", "cpu_fallback_evaluations", "gpu_batches", "gpu_failed_batches", "gpu_transfer_seconds", "gpu_kernel_seconds", "gpu_elapsed_seconds", "generations", "elapsed_seconds")
            if not isinstance(stats, dict):
                errors.append(f"optimization.{name}.statistics must be an object")
                continue
            for key in stat_keys:
                required(stats, key, f"optimization.{name}.statistics")
            for key in ("seed", "evaluations", "successful_evaluations", "failed_evaluations", "cache_hits", "gpu_fallbacks", "gpu_requested_evaluations", "gpu_executed_evaluations", "gpu_successful_evaluations", "gpu_failed_evaluations", "cpu_fallback_evaluations", "gpu_batches", "gpu_failed_batches", "generations"):
                if not isinstance(stats.get(key), int) or stats[key] < 0:
                    errors.append(f"optimization.{name}.statistics.{key} must be a non-negative integer")
            for key in ("gpu_transfer_seconds", "gpu_kernel_seconds", "gpu_elapsed_seconds", "elapsed_seconds"):
                if not _finite(stats.get(key)) or stats[key] < 0:
                    errors.append(f"optimization.{name}.statistics.{key} must be finite and non-negative")
            if _finite(stats.get("evaluations")) and stats.get("successful_evaluations", 0) + stats.get("failed_evaluations", 0) != stats["evaluations"]:
                errors.append(f"optimization.{name}.statistics evaluation count mismatch")
            if _finite(stats.get("gpu_executed_evaluations")) and stats.get("gpu_successful_evaluations", 0) + stats.get("gpu_failed_evaluations", 0) != stats["gpu_executed_evaluations"]:
                errors.append(f"optimization.{name}.statistics GPU evaluation count mismatch")
            if stats.get("gpu_executed_evaluations", 0) > stats.get("gpu_requested_evaluations", 0):
                errors.append(f"optimization.{name}.statistics GPU requested/executed count mismatch")
            if stats.get("gpu_failed_batches", 0) > stats.get("gpu_batches", 0):
                errors.append(f"optimization.{name}.statistics GPU batch count mismatch")

    decision = document.get("decision")
    if not isinstance(decision, dict):
        errors.append("decision must be an object")
    else:
        for key in ("batch32_gpu_faster", "batch128_gpu_faster"):
            if not isinstance(decision.get(key), bool):
                errors.append(f"decision.{key} must be boolean")
        for target in (32, 128):
            key = f"batch{target}_gpu_faster"
            row = next((candidate for candidate in batches
                        if isinstance(candidate, dict) and candidate.get("batch_size") == target), None)
            if row is not None and isinstance(decision.get(key), bool):
                cpu_per_candidate = row.get("cpu_per_candidate_ms")
                gpu_per_candidate = row.get("gpu_per_candidate_ms")
                cpu_throughput = row.get("cpu_candidates_per_second")
                gpu_throughput = row.get("gpu_candidates_per_second")
                if all(_finite(value) for value in (cpu_per_candidate, gpu_per_candidate,
                                                    cpu_throughput, gpu_throughput)):
                    expected = gpu_per_candidate < cpu_per_candidate and gpu_throughput > cpu_throughput
                    if decision[key] is not expected:
                        errors.append(f"decision.{key} does not match measured per-candidate timing/throughput")
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("document", help="JSON evidence file")
    args = parser.parse_args(argv)
    try:
        with open(args.document, encoding="utf-8") as stream:
            document = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        print(f"invalid JSON: {error}", file=sys.stderr)
        return 2
    errors = validate_document(document)
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print(f"valid B-T4 schema v1: {args.document}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
