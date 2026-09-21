# A-T2b — Strengthen optimization regression tests

## Scope

Strengthened the single-objective direction regression, unknown non-`std::exception`
batch isolation regression, and coilgun adapter behavior coverage for independent
`CoilTurns` and `TriggerValue` bindings. Changes are limited to
`tests/test_optimization_single.cpp`, `tests/test_coilgun_optimization.cpp`, and
this report. No production change remains in the final diff.

## Test-effectiveness mutation evidence

Each mutation was applied temporarily, rebuilt, run against its focused test case,
and restored before the next mutation.

### Single-objective direction

Mutation: in `src/optimization/comparator.cpp`, forced the objective definition
used by `FeasibilityComparator` to `maximize == false`, making every objective
comparison minimize.

Command:

```text
cmake --build --preset cpu-debug --target test_optimization_single -j2 && \
build/cpu-debug/tests/test_optimization_single \
  --test-case='single objective optimizer honors maximize and minimize directions'
```

Expected RED result (exit code 1):

```text
CHECK( max_value == *std::max_element(...) ) is NOT correct!
values: CHECK( 0.497917 == 8.35125 )
CHECK( max_value > min_value ) is NOT correct!
values: CHECK( 0.497917 >  0.497917 )
test cases: 1 | 0 passed | 1 failed
```

The same seeded candidate stream is now compared against both deterministic
extrema, so a direction inversion cannot pass.

### Unknown singleton exception isolation

Mutation: changed `UnknownExceptionEvaluator` so its singleton retry never threw.
The batch still throws the non-`std::exception` value `42`, but all singleton
retries then incorrectly succeed.

Command:

```text
cmake --build --preset cpu-debug --target test_optimization_single -j2 && \
build/cpu-debug/tests/test_optimization_single \
  --test-case='optimizer isolates unknown batch exceptions and preserves successful siblings'
```

Expected RED result (exit code 1):

```text
CHECK( successful_evaluations == population_size - 1 ) is NOT correct!
values: CHECK( 4 == 3 )
CHECK( failed_evaluations == 1 ) is NOT correct!
values: CHECK( 0 == 1 )
CHECK( singleton_calls == population_size ) is NOT correct!
values: CHECK( 0 == 4 )
REQUIRE( successful_values.size() == population_size - 1 ) is NOT correct!
values: REQUIRE( 4 == 3 )
test cases: 1 | 0 passed | 1 failed
```

The real test forces the first singleton retry to throw and checks exact
evaluation, success, and failure counts, `MaxGenerations` termination, retry
count, and usability of the remaining successful siblings through the selected
best result.

### `CoilTurns` binding

Mutation: made the `CoilTurns` assignment in the adapter's coil-spec decoding a
no-op.

Command:

```text
cmake --build --preset cpu-debug --target test_coilgun_optimization -j2 && \
build/cpu-debug/tests/test_coilgun_optimization \
  --test-case='coil turns binding changes the physical objective'
```

Expected RED result (exit code 1):

```text
CHECK( std::abs(low_velocity - high_velocity) > 0.0 ) is NOT correct!
values: CHECK( 0 > 0 )
test cases: 1 | 0 passed | 1 failed
```

### `TriggerValue` binding

Mutation: made the adapter's effective trigger-value assignment a no-op.

Command:

```text
cmake --build --preset cpu-debug --target test_coilgun_optimization -j2 && \
build/cpu-debug/tests/test_coilgun_optimization \
  --test-case='trigger value binding changes the physical objective'
```

Expected RED result (exit code 1):

```text
CHECK( std::abs(immediate_velocity - delayed_velocity) > 0.0 ) is NOT correct!
values: CHECK( 0 > 0 )
test cases: 1 | 0 passed | 1 failed
```

Both adapter tests keep all other configuration and candidate values fixed,
require successful finite objectives, and assert a nonzero physical effect.

## Final GREEN verification

Focused executables:

```text
build/cpu-debug/tests/test_optimization_single
```

```text
test cases: 11 | 11 passed | 0 failed
assertions: 62 | 62 passed | 0 failed
```

```text
build/cpu-debug/tests/test_coilgun_optimization
```

```text
test cases: 13 | 13 passed | 0 failed
assertions: 59 | 59 passed | 0 failed
```

Full optimization-named CTest set:

```text
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug -R 'optimization|coilgun_optimization' --output-on-failure
```

```text
100% tests passed out of 12
Total Test time (real) = 0.31 sec
```

## Concerns

None. All temporary mutations were restored, and the final source diff contains
only the two relevant test files and this report.
