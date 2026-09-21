# A-T4b — Stabilize optimization problem contract

## RED/GREEN evidence

The new `test_optimization_problem_spec` target was added before the production
implementation. The first build failed as expected because `ProblemSpec` and
the problem-only optimizer constructors did not exist. For this review-fix
cycle, three regressions were added before the fixes: a const dual problem
initially dispatched its legacy non-const batch method, spec Auto→NSGA-II
stagnation initially repaired and evaluated one generation, and the scaled
Penalty comparison initially selected no valid declared result. Each failed in
the expected way before production changes and passed after them.

The focused target now passes all eight contract tests, including exact
zero-call assertions for the pre-initialization configuration rejection.

## Contract and compatibility

`ProblemSpec` owns one `VariableSchema`, ordered objective definitions, ordered
constraint definitions, and a callable `RepairPolicy`. Its accessors return only
const references. Construction validates all definitions and unique IDs;
explicitly empty repair policies are rejected. `ProblemSpec::repair` invokes the
custom policy and then runs the owned schema repair, preserving dimensionality,
bounds, integer, and enum invariants. The default policy delegates to
`VariableSchema::repair`.

`OptimizationProblem` retains its default constructor for legacy subclasses and
adds an owning spec constructor. The new `GeneticOptimizer(problem, ...)` path
rejects a legacy problem immediately with `std::invalid_argument`, while all
schema-first constructors and `optimize_single_objective` remain available.
The new `(ProblemSpec, BatchEvaluator)` and shared-pointer constructors support
batch callers without a duplicate schema.

Spec-aware runs use the declared schema for initialization and genetic operators,
apply repair to initial candidates and offspring, freeze objective and
constraint metadata, route from the declared objective count, and use declared
direction/scale definitions in scalar and NSGA-II comparisons. Successful
evaluator metadata mismatches terminate structurally with
`TerminationReason::ConfigurationError`; legacy schema-first runs retain their
dynamic schema discovery.

## Coilgun and batch-path evidence

`CoilgunOptimizationProblem` now publishes an owned spec containing the
`muzzle_velocity` maximizing objective at scale `1.0` and its configured metric
constraints in order. `schema()` reads from that authoritative spec. The
optional `BatchEvaluator::evaluate_batch_const` hook is safe and unavailable by
default; const problem-only optimization therefore serializes through
`OptimizationProblem::evaluate` unless a type opts into the hook. The const
dual regression observes zero legacy batch calls and exactly one serial call
per candidate. Non-const dual problems still dispatch their legacy batch
override, and Coilgun's const-native adapter opts into the hook while retaining
its existing const and non-const overloads; its GPU callback regression remains
green.

Spec-aware scalar runs construct a run-local `FeasibilityComparator` from the
caller's strategy/penalty weight and the declared single objective, so scale is
used consistently for best-candidate selection, elites, tournaments, and the
incumbent comparison. A deterministic soft-penalty case proves scale `100`
reverses the winner selected with implicit scale `1`.

When a spec's Auto route resolves to NSGA-II, a nonzero stagnation limit now
returns `ConfigurationError` before population initialization/repair or
evaluator work. Legacy schema-first Auto routing retains its one-batch
objective discovery behavior.

## Verification

Commands run successfully:

```text
cmake --preset cpu-debug
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug -R 'test_optimization|test_coilgun_optimization' --output-on-failure
ctest --preset cpu-debug --output-on-failure
git diff --check
```

The focused optimization/Coilgun selection ran 13/13 tests, and the complete
CPU preset ran 34/34 tests with zero failures.

## Review-fix regression

The base `BatchEvaluator` no longer exposes a direct const `evaluate_batch`
overload, avoiding a silent empty-vector result when the optional
`evaluate_batch_const` hook is unavailable. A compile-time concept proves that
the base has no direct const batch call while
`CoilgunOptimizationProblem` retains its concrete const overload; that
overload is no longer marked as a base override. `ProblemBatchEvaluator` uses
only the optional hook on the const path and preserves serial fallback when it
returns `nullopt`.

## Deferred concerns

Public bilingual API manuals and install/export documentation remain deferred
to the later documentation task as requested.
