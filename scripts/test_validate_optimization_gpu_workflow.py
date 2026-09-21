import copy
import importlib.util
import json
import math
import pathlib
import unittest


_MODULE_PATH = pathlib.Path(__file__).with_name("validate_optimization_gpu_workflow.py")
_SPEC = importlib.util.spec_from_file_location("workflow_validator", _MODULE_PATH)
_MODULE = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
_SPEC.loader.exec_module(_MODULE)


def valid_document():
    sample = {
        "schema_version": 1,
        "source_revision": "abc1234",
        "worktree_state": "clean",
        "toolchain": {"compiler": "GNU 16.2.1", "cmake": "4.4.3", "build_type": "release"},
        "workload_provenance": {"name": "optimization-gpu-workflow", "configuration": "fixed-geometry-euler-full", "command": "bench_gpu_optimization_workflow"},
        "hardware": {"gpu_name": "RTX", "driver": "615.71.09", "compute_capability": "12.0", "cuda_toolkit": "13.4"},
        "workload": {"batch_sizes": [1, 8, 32, 128], "warmups": 2, "repetitions": 5, "precision": "full", "seed": 7},
        "batches": [{
            "batch_size": 1,
            "cpu_samples_ms": [2.0, 2.1, 2.2, 2.1, 2.0],
            "gpu_samples_ms": [1.0, 1.1, 1.2, 1.1, 1.0],
            "cpu_median_ms": 2.1,
            "gpu_median_ms": 1.1,
            "cpu_per_candidate_ms": 2.1,
            "gpu_per_candidate_ms": 1.1,
            "cpu_candidates_per_second": 1000.0 / 2.1,
            "gpu_candidates_per_second": 1000.0 / 1.1,
            "speedup": 2.1 / 1.1,
            "cpu_rows": 1,
            "gpu_rows": 1,
            "order_ok": True,
            "numerical_ok": True,
            "gpu_executed": True,
            "backend": "direct",
            "requested_backend": "direct",
            "solver": "cusolver",
            "precision": "full",
            "requested_precision": "full",
            "fallback_events": 0,
            "gpu_report_time_ms": 1.0,
            "transfer_time_ms": 0.1,
        }],
        "optimization": {"same_seed": True, "population_size": 8, "max_generations": 4,
            "cpu": {"termination": "maximum generations", "feasible": True, "best_voltage": 500.0,
                    "best_objective": 2.0, "reference_objective": 2.0, "reference_error": 0.0,
                    "reference_tolerance": 0.00000205, "reference_valid": True,
                    "statistics": {"seed": 7, "evaluations": 8, "successful_evaluations": 8,
                        "failed_evaluations": 0, "cache_hits": 0, "gpu_fallbacks": 0,
                        "gpu_requested_evaluations": 0, "gpu_executed_evaluations": 0,
                        "gpu_successful_evaluations": 0, "gpu_failed_evaluations": 0,
                        "cpu_fallback_evaluations": 0, "gpu_batches": 0, "gpu_failed_batches": 0,
                        "gpu_transfer_seconds": 0.0, "gpu_kernel_seconds": 0.0,
                        "gpu_elapsed_seconds": 0.0, "generations": 4, "elapsed_seconds": 1.0}},
            "gpu": {"termination": "maximum generations", "feasible": True, "best_voltage": 500.0,
                    "best_objective": 2.0, "reference_objective": 2.0, "reference_error": 0.0,
                    "reference_tolerance": 0.00000205, "reference_valid": True,
                    "statistics": {"seed": 7, "evaluations": 8, "successful_evaluations": 8,
                        "failed_evaluations": 0, "cache_hits": 0, "gpu_fallbacks": 0,
                        "gpu_requested_evaluations": 8, "gpu_executed_evaluations": 8,
                        "gpu_successful_evaluations": 8, "gpu_failed_evaluations": 0,
                        "cpu_fallback_evaluations": 0, "gpu_batches": 1, "gpu_failed_batches": 0,
                        "gpu_transfer_seconds": 0.1, "gpu_kernel_seconds": 0.2,
                        "gpu_elapsed_seconds": 0.3, "generations": 4, "elapsed_seconds": 0.4}}},
        "decision": {"batch32_gpu_faster": True, "batch128_gpu_faster": True},
    }
    sample["batches"] = [copy.deepcopy(sample["batches"][0]) for _ in range(4)]
    for size, row in zip((1, 8, 32, 128), sample["batches"]):
        row["batch_size"] = size
        row["cpu_rows"] = size
        row["gpu_rows"] = size
        row["representative_rows"] = 1
        row["cpu_per_candidate_ms"] = row["cpu_median_ms"] / size
        row["gpu_per_candidate_ms"] = row["gpu_median_ms"] / size
        row["cpu_candidates_per_second"] = 1000.0 * size / row["cpu_median_ms"]
        row["gpu_candidates_per_second"] = 1000.0 * size / row["gpu_median_ms"]
        row["speedup"] = row["cpu_median_ms"] / row["gpu_median_ms"]
    return sample


class WorkflowValidatorTests(unittest.TestCase):
    def test_accepts_complete_finite_gpu_workflow(self):
        self.assertEqual(_MODULE.validate_document(valid_document()), [])

    def test_rejects_missing_required_key(self):
        document = valid_document()
        del document["batches"]
        self.assertTrue(any("batches" in error for error in _MODULE.validate_document(document)))

    def test_rejects_missing_or_malformed_provenance(self):
        document = valid_document()
        document["source_revision"] = ""
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("source_revision" in error for error in errors))
        document = valid_document()
        document["toolchain"]["compiler"] = ""
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("toolchain" in error for error in errors))
        document = valid_document()
        del document["workload_provenance"]
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("workload_provenance" in error for error in errors))

    def test_rejects_non_finite_or_negative_time(self):
        document = valid_document()
        document["batches"][0]["gpu_samples_ms"][0] = math.nan
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("finite" in error for error in errors))
        document = valid_document()
        document["batches"][0]["cpu_samples_ms"][0] = -1.0
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("non-negative" in error for error in errors))

    def test_rejects_zero_timing_sample(self):
        document = valid_document()
        document["batches"][0]["gpu_samples_ms"][0] = 0.0
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("greater than zero" in error for error in errors))

    def test_rejects_fallback_or_false_gpu_execution(self):
        document = valid_document()
        document["batches"][0]["backend"] = "fallback"
        document["batches"][0]["gpu_executed"] = False
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("backend" in error for error in errors))
        self.assertTrue(any("gpu_executed" in error for error in errors))

    def test_rejects_inconsistent_counts_and_samples(self):
        document = valid_document()
        document["batches"][0]["gpu_rows"] = 2
        document["batches"][0]["gpu_samples_ms"] = [1.0]
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("count" in error for error in errors))
        self.assertTrue(any("samples" in error for error in errors))

    def test_rejects_missing_backend_solver_precision_and_optimizer_fields(self):
        document = valid_document()
        for key in ("backend", "solver", "precision"):
            mutated = copy.deepcopy(document)
            del mutated["batches"][0][key]
            self.assertTrue(any(key in error for error in _MODULE.validate_document(mutated)))
        mutated = copy.deepcopy(document)
        del mutated["optimization"]["gpu"]["statistics"]
        self.assertTrue(any("statistics" in error for error in _MODULE.validate_document(mutated)))
        mutated = copy.deepcopy(document)
        del mutated["optimization"]["gpu"]["best_objective"]
        self.assertTrue(any("best_objective" in error for error in _MODULE.validate_document(mutated)))

    def test_rejects_derived_timing_mismatch(self):
        document = valid_document()
        for key, value in (("cpu_median_ms", 9.0), ("gpu_per_candidate_ms", 9.0),
                           ("gpu_candidates_per_second", 9.0), ("speedup", 9.0)):
            mutated = copy.deepcopy(document)
            mutated["batches"][0][key] = value
            self.assertTrue(any("derived" in error or key in error
                                for error in _MODULE.validate_document(mutated)))

    def test_rejects_reference_error_or_validity_mismatch(self):
        document = valid_document()
        mutated = copy.deepcopy(document)
        mutated["optimization"]["gpu"]["reference_error"] = 1.0
        self.assertTrue(any("reference" in error for error in _MODULE.validate_document(mutated)))
        mutated = copy.deepcopy(document)
        mutated["optimization"]["gpu"]["reference_valid"] = False
        self.assertTrue(any("reference_valid" in error for error in _MODULE.validate_document(mutated)))

    def test_rejects_optimizer_count_relationship(self):
        document = valid_document()
        document["optimization"]["gpu"]["statistics"]["gpu_successful_evaluations"] = 7
        self.assertTrue(any("count" in error for error in _MODULE.validate_document(document)))

    def test_rejects_inconsistent_derived_decisions(self):
        document = valid_document()
        document["decision"]["batch32_gpu_faster"] = False
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("batch32_gpu_faster" in error for error in errors))

        document = valid_document()
        row = document["batches"][2]
        row["gpu_samples_ms"] = [3.0, 3.1, 3.2, 3.1, 3.0]
        row["gpu_median_ms"] = 3.1
        row["gpu_per_candidate_ms"] = 3.1 / 32.0
        row["gpu_candidates_per_second"] = 1000.0 * 32.0 / 3.1
        row["speedup"] = 2.1 / 3.1
        errors = _MODULE.validate_document(document)
        self.assertTrue(any("batch32_gpu_faster" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
