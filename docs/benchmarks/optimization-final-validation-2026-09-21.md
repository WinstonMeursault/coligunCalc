# Optimization integration validation — 2026-09-21

## Decision

The v4.2.0 optimization workflow passed its final integration gate. The
validated scope includes the generic GA/NSGA-II framework, the coilgun problem
adapter, the fixed-geometry CUDA batch evaluator, public headers and package
exports, bilingual API documentation, and the production benchmark workflow.

This is a dated validation record, not a promise that future CTest discovery
will retain the same test counts. Re-run the repository presets to validate a
new revision.

## Verification snapshot

The final clean validation ran all four supported development presets:

| Preset | Result at the validated revision |
|---|---:|
| `cpu-debug` | 35/35 passed |
| `cpu-release` | 35/35 passed |
| `cuda-debug` | 53/53 passed |
| `cuda-release` | 53/53 passed |

The CUDA configurations included 19 GPU-labelled tests. Independent CUDA
Release executions also passed:

- `test_cuda_batch_evaluator`: 24/24 cases, 244/244 assertions.
- `test_gpu_optimization_baseline`: 1/1 case, 154/154 assertions.

The final read-only whole-branch review approved the implementation with no
Critical, Important, or Minor findings. It specifically checked evaluator
lifetime ownership, cache identity isolation, callback schema validation,
fallback behavior, automatic GA/NSGA-II routing, Pareto result semantics,
installation, documentation, and benchmark contracts.

## Real-GPU evidence

Validation used an NVIDIA GeForce RTX 5080 Laptop GPU (compute capability
12.0), driver 615.71.09, CUDA 13.4.59, GCC 16.2.1, and CMake 4.4.3. Every
accepted production measurement required both
`ExecutionReport::gpu_executed == true` and a resolved backend other than
`BackendMode::Fallback`.

The final CUDA Release benchmark reported direct, Full-precision GPU execution,
stable candidate ordering, CPU numerical parity for supported metrics, and no
fallback events. Per-candidate GPU throughput exceeded CPU throughput at batch
sizes 32 and 128. Batch size 1 remained slower because launch and transfer
overhead could not be amortized.

The GPU optimizer completed 32/32 successful, GPU-executed evaluations in four
batches, with no failed batches, CPU repairs, or fallback events. Its selected
result passed a CPU Reference recheck.

Exact workload parameters, raw timings, schema rules, and reproduction commands
remain in [GPU optimization workflow benchmark](optimization-gpu-workflow-2026-09-17.md).
The diagnostic simulation comparison is recorded in
[CUDA optimization numerical baseline](optimization-cuda-numerical-baseline-2026-09-16.md).

## Public API and packaging

The English and Chinese API documents passed the repository equivalence
checker and its negative self-test. CPU and CUDA Release installations passed
the installed-header split checks. Independent out-of-tree consumers found the
package with `find_package(coilgun CONFIG REQUIRED)`, then built and ran against
only their installed prefixes:

- CPU consumers use `coilgun::coilgun` and the CPU optimization umbrella.
- CUDA consumers use `coilgun::coilgun_cuda` and the CUDA batch evaluator.

Evaluator snapshots require genuine `std::shared_ptr` ownership and retain a
strong owner for saved or concurrent calls. `CudaBatchEvaluator` owns an
immutable copy of the source problem configuration, so snapshots and CPU
fallback do not borrow the source problem. Callback results are validated
against the complete declared objective and constraint schema before caching.
The authoritative contract is in [API.md](../API.md), mirrored by
[API_cn.md](../API_cn.md).

## Durable limitations

- The production CUDA optimizer supports fixed shared geometry, Euler stepping,
  and `OptimizationLevel::Full`; arbitrary-geometry candidates, CUDA RK4,
  thermal optimization metrics, and non-Full candidates are outside its
  contract.
- CUDA `PeakCurrent` remains a diagnostic-only metric. The numerical baseline
  measured a reproducible 5–8% shortfall relative to CPU Full, so it is rejected
  as a production objective or constraint instead of hiding the discrepancy by
  widening tolerances.
- `OptimizationLevel::LookupTable` remains a compatibility value and currently
  follows the same runtime path as `Reference`.
- `BackendMode::Persistent` is not a supported execution path and resolves to
  fallback.
- Additional algorithm families such as differential evolution, CMA-ES,
  Bayesian optimization, and MOEA/D were not part of this release.

## Revalidation

Use the repository presets rather than copying the historical counts above:

```sh
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug --output-on-failure

cmake --preset cuda-release
cmake --build --preset cuda-release
ctest --preset cuda-release --output-on-failure
```

For benchmark reproduction and validation, follow the commands in the linked
workflow benchmark document.
