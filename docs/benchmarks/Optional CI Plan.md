# Optional CI Performance Feedback Plan

This is a feedback plan, not a release gate. It proposes how CI can expose
correctness and performance observations without making machine-dependent
timings block normal development.

## Jobs

| Job | Runner | Command/preset | Frequency | Required result |
|---|---|---|---|---|
| CPU quick | Linux CPU | `ctest --preset cpu-quick --output-on-failure` | Every PR | Tests pass |
| CPU parallel | Linux CPU | `ctest --preset cpu-parallel --output-on-failure` | Every PR or merge queue | Tests pass; wall time is a report |
| CPU slow/integration | Linux CPU | `ctest --preset cpu-integration --output-on-failure` | Merge queue/nightly | Tests pass |
| CUDA quick | Pinned CUDA runner | `ctest --preset cuda-quick --output-on-failure` | CUDA PRs, nightly | Tests pass; check actual GPU fields when device exists |
| CUDA integration | Pinned CUDA runner | `ctest --preset cuda-integration --output-on-failure` | Nightly or CUDA changes | Tests pass; check GPU/fallback semantics |
| GPU benchmark | Same pinned CUDA runner | Build and run `bench_gpu_engine` | Scheduled/manual | Artifact is complete and finite; timing is advisory |

Existing default `cpu-*` and `cuda-*` test presets remain the complete
regression entry points. Library-only builds can use
`cmake --build --preset cpu-release-library` or
`cmake --build --preset cuda-release-library` when tests are not needed.

## Concurrency and Artifacts

- Allow at most one GPU job per physical device. Keep CTest's `gpu`
  `RESOURCE_LOCK` and serialize the benchmark job with integration.
- Limit CPU parallel jobs to the preset's bounded `jobs=4`; do not infer a
  runtime optimization from CI test wall time.
- Upload complete CTest output, benchmark stdout, commit, preset, compiler,
  GPU, driver, CUDA runtime, and runner image metadata.
- Keep raw benchmark artifacts by commit and hardware identity. Do not replace
  an older measurement with a newer machine's result.

## Validation Policy

The job fails for build/test failure, non-finite benchmark fields, malformed
headers, missing required metadata, or a fallback row incorrectly marked as
`gpu_executed=yes` or given a numeric GPU speedup. A requested GPU backend is
not proof of GPU execution.

Timing comparisons are advisory. Flag a regression only when the same workload,
preset, hardware identity, and protocol have at least three repeats. Compare
medians and dispersion; do not fail a pull request for scheduling noise. Tune
an alert threshold only after history accumulates, and never turn fallback
latency into GPU speedup.

## Benchmark Artifact Checklist

1. Confirm the expected schema header before parsing Markdown rows.
2. Confirm the expected phase set and repeat count.
3. Confirm all measured timing fields are finite.
4. Partition rows by `gpu_executed` before computing speedup.
5. Report raw samples, median, dispersion, fallback reasons, and Graph rebuild totals.
6. Link the artifact to the exact commit and hardware metadata.

This plan is intentionally non-blocking for runtime optimization. Changes to
the protocol or alert thresholds require a new dated schema record and an
explicit before/after comparison.
