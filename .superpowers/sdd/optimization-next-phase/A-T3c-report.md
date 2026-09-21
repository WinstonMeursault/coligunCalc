# A-T3c — Unify optimization termination reasons

## Scope

Unified optimization termination reporting around `TerminationReason`. Evaluation
budgets now report `MaxEvaluations` directly, the legacy
`GeneticTerminationReason` name is a type alias, and
`genetic_termination_reason` returns the structured reason without inspecting
diagnostic text. Unsupported NSGA-II no-improvement termination now fails fast
according to when its strategy can be known.

## TDD RED evidence

Tests were added before production edits. The focused executables were rebuilt
against the old implementation and failed with the expected behavior:

```text
test_optimization_single: 2 failed, 10 passed
CHECK( budget_result.termination.reason != TerminationReason::MaxGenerations )
  values: CHECK( 1 != 1 )
CHECK( result.termination.reason != TerminationReason::MaxGenerations )
  values: CHECK( 1 != 1 )

test_optimization_routing: 2 failed, 5 passed
explicit NSGA-II: reason was MaxGenerations, evaluator calls 3, candidates 36
Auto NSGA-II: reason was MaxGenerations, evaluator calls 3, candidates 36,
  generations 3, and a non-empty Pareto front
```

The RED failures demonstrated that budgets were encoded as
`MaxGenerations`, while NSGA-II stagnation was ignored and evolution
continued.

## Implementation and compatibility

- Appended `TerminationReason::MaxEvaluations` so existing enumerator values
  remain unchanged.
- Replaced the duplicate `GeneticTerminationReason` enum with
  `using GeneticTerminationReason = TerminationReason`.
- Kept `genetic_termination_reason(const OptimizationTermination&)` as a
  source-compatible helper that directly returns `termination.reason`.
- Added structured `to_string` coverage for all unified values, including
  `Cancelled` and `MaxEvaluations`.
- Explicit `SelectionStrategy::NSGA2` with a nonzero
  `max_no_improvement_generations` is rejected before evaluator work.
- `SelectionStrategy::Auto` is rejected immediately after the first objective
  schema resolves to NSGA-II, before NSGA-II selection or evolution.
- Single-objective stagnation behavior and ordinary NSGA-II routing remain
  unchanged.

## GREEN evidence

Focused tests:

```text
build/cpu-debug/tests/test_optimization_types
6 test cases | 28 assertions passed

build/cpu-debug/tests/test_optimization_single
12 test cases | 64 assertions passed

build/cpu-debug/tests/test_optimization_routing
7 test cases | 39 assertions passed
```

The fail-fast tests now prove the evaluator-call boundaries exactly:

- Explicit NSGA-II: `0` evaluator calls and `0` candidates evaluated.
- Auto multiobjective NSGA-II: exactly `1` evaluator call over the initial
  `12` candidates, `1` generation recorded, no Pareto evolution.

## Verification

- `cmake --build --preset cpu-debug -j2` — passed.
- `ctest --preset cpu-debug -R 'optimization|coilgun_optimization' --output-on-failure`
  — 12/12 tests passed.
- `git diff --check` — clean.

No separate project lint or typecheck command is configured in this worktree.

## Deferred concerns

No multiobjective convergence metric was added. The existing generic
termination message remains diagnostic only; callers should use the structured
`reason` field (or the compatibility helper) for control flow.
