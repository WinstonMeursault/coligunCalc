# A-T4a — Harden optimization result selectors

## Baseline

The characterization selector suite was green before production edits:

```text
ctest --preset cpu-debug -R '^test_optimization_selectors$' --output-on-failure
1/1 tests passed
```

The focused tests were strengthened before the refactor with coverage for:

- first-candidate tie behavior across every built-in selector;
- copied `OptimizationResult` values;
- callback and custom-subclass compatibility; and
- preservation of the complete Pareto front after selection.

Those strengthened tests also passed against the characterization implementation (1/1).

## Structural RED/GREEN evidence

Before the production change, the required structural check found address-derived
indices in `src/optimization/selectors.cpp` and returned exit code 0:

```text
rg -n -U '&(lhs|rhs)|pareto_front\.data\(\)|lhs\s*-|rhs\s*-' src/optimization/selectors.cpp
140: const auto index = static_cast<std::size_t>(&lhs - result.pareto_front.data());
141: const auto other = static_cast<std::size_t>(&rhs - result.pareto_front.data());
158: const auto left_index = static_cast<std::size_t>(&lhs - result.pareto_front.data());
159: const auto right_index = static_cast<std::size_t>(&rhs - result.pareto_front.data());
structural_red_exit=0
```

After the refactor, the same check returned no matches and passed:

```text
structural_green=PASS
```

The internal `choose` helper now passes explicit `std::size_t` row indices to
every built-in selector comparator. Normalized score matrices are indexed only
with those explicit indices; no address-to-index recovery or pointer subtraction
remains in built-in selector logic.

## Behavioral verification

```text
ctest --preset cpu-debug -R '^test_optimization_selectors$' --output-on-failure
1/1 tests passed

ctest --preset cpu-debug -R 'optimization' --output-on-failure
12/12 tests passed

git diff --check
passed
```

The focused suite confirms stable first-candidate ties, normalized selector
results, copied-result selection, custom callback/subclass behavior, and
non-mutation of the Pareto front.

## Compatibility

No public selector headers or signatures changed. `RepresentativeSelector::select`
still returns `std::optional<Candidate>`, `OptimizationResult::select_representative`
is unchanged, and `CallbackSelector` plus existing custom subclasses continue to
use the same API and behavior.

## Changed files

- `src/optimization/selectors.cpp` — switched the internal built-in selector
  protocol from candidate references to explicit row indices.
- `tests/test_optimization_selectors.cpp` — added tie, copied-result,
  compatibility, and non-mutation coverage.
- `.superpowers/sdd/optimization-next-phase/A-T4a-report.md` — this report.

## Deferred concerns

None within A-T4a scope. Public API exposure and documentation synchronization
remain assigned to later tasks.
