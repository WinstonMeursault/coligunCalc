# A-T5 — Document optimization public API

## Scope and outcome

The optimization API manuals now describe the actual headers and implementation
on the B-T3 base in both languages. The sections cover `ProblemSpec`, variable /
objective / constraint definitions and repair, automatic scalar-GA versus
NSGA-II routing, Pareto-only multi-objective results and explicit selectors,
validation/error behavior, termination and seed semantics, cache identity and
statistics ownership, the `CoilgunOptimizationProblem` adapter, and the
fixed-geometry `CudaBatchEvaluator` support/fallback matrix. CPU-only versus
CUDA umbrella, link, and install behavior is stated explicitly.

The original plan now uses the real `cpu-debug`, `cpu-release`, `cuda-debug`, and
`cuda-release` names. Its implementation addendum distinguishes the original
design from the implemented fixed-geometry-first CUDA scope. The 2026-09-08
benchmark keeps all raw measurements, but labels the no-CUDA-evaluator text as
historical and points current throughput work to B-T4.

## Executable validation

`tests/test_optimization_integration.cpp` now asserts fixed-seed deterministic
equality between two independent runs, hard-constraint feasibility, a finite
Full terminal velocity within the named current post-RNG baseline envelope,
strict Full / Reference distinction, and the existing
`5e-8 + 1e-6 * abs(reference_velocity)` bound. The current Full baseline is
`0.0096453039804834419 m/s`; its envelope uses the same tight absolute/relative
terms (`5e-8 + 1e-6 * abs(baseline)`) as the Reference check. Platform-brittle
exact variable/objective literals were removed. The historical
`0.009764939958 m/s` value remains unchanged in the 2026-09-08 benchmark; the
intentional evolution is attributable to the post-RNG stream implementation,
not a relaxed numerical tolerance. The test also checks seed,
submitted/success/failure counts, cache/fallback ownership, GPU counters, budget
skips, generation count, and non-negative finite run time.

TDD evidence:

- RED: the new tight-envelope test intentionally used the historical value and
  failed (`0.0096453` versus `0.009764939958`; error `0.000119636` versus
  tolerance `5.97649e-08`).
- GREEN: the baseline was recorded from the current post-RNG run, exact
  cross-run equality was retained, and both Full and Reference checks use the
  named `5e-8 + 1e-6 * abs(value)` formula.

Review-fix TDD also covered packaging and documentation gates:

- RED: the new CPU install verifier found
  `include/coilgun/optimization/cuda_batch_evaluator.hpp` in a CPU-only install.
- GREEN: `CMakeLists.txt` excludes that header from the generic optimization
  install and installs it only when `COILGUN_ENABLE_CUDA` is enabled. The
  checker was strengthened with semantic marker maps, ordered routing/fallback
  clauses, identical executable examples, and JSON-backed preset validation;
  its negative omission self-test passes.

`OptimizationStatistics::gpu_fallbacks` is documented with its compatibility
alias rule: a non-zero explicit collector `gpu_fallbacks` wins; otherwise the
legacy collector `fallbacks` value is used, with no double counting.

The install split is verified by `scripts/verify_install_headers.py`: CPU-only
installs omit `cuda_batch_evaluator.hpp`, while CUDA installs include it and
export `coilgun::coilgun_cuda`.

## Bilingual equivalence and verification

The maintained checker `scripts/check_api_optimization_docs.py` checks explicit
semantic contract marker maps in both optimization sections, policy-marker
ordering, matching executable examples, and preset names loaded from the real
`CMakePresets.json`. Its negative omission self-test also proves that removing
one marker fails the gate. Current output:

```text
negative self-test: PASS
API optimization bilingual check: PASS
semantic_contracts=6 common_tokens=30
example_blocks=1 per language
actual_presets=configure:6 build:6 test:12
```

Verification commands and results:

```text
cmake --preset cpu-debug && cmake --build --preset cpu-debug -j2
  configure succeeded; build completed with no work remaining
cpu_prefix=$(mktemp -d) && cmake --install build/cpu-debug --prefix "$cpu_prefix"
python3 scripts/verify_install_headers.py "$cpu_prefix"
  install header split: PASS (CPU-only)
cmake --build --preset cuda-debug --target coilgun_cuda -j2
cuda_prefix=$(mktemp -d) && cmake --install build/cuda-debug --prefix "$cuda_prefix"
python3 scripts/verify_install_headers.py "$cuda_prefix" --cuda
  install header split: PASS (CUDA); coilgun::coilgun_cuda export: PASS
ctest --preset cpu-debug -R '^test_optimization_integration$' --output-on-failure
  1/1 passed
ctest --preset cpu-debug -R '^(test_optimization_.*|test_coilgun_optimization)$' --output-on-failure
  14/14 passed
cmake --list-presets && ctest --list-presets
  real configure/build and test preset names listed
git diff --check
  passed
```

The current preset inventories were checked with `cmake --list-presets` and
`ctest --list-presets`. The configure/build presets are `cpu-debug`,
`cpu-release`, `cuda-debug`, `cuda-release`, `cpu-release-library`, and
`cuda-release-library`; matching CTest presets exist for the four debug/release
builds, plus the quick/parallel/slow/integration filters.

## Deferred work and limits

The durable ledgers record the B-T1 `PeakCurrent` discrepancy, immutable or
revalidated cache identity, arbitrary-geometry CUDA candidates, thermal CUDA
optimization metrics, and future DE/CMA-ES/Bayesian/MOEA-D families. This task
does not run or claim B-T4 benchmark success and does not claim whole-branch
completion.
