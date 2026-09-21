# A-T3a — Validate thermal optimization constraints

## Scope

Added focused adapter tests for maximum-temperature constraint configuration and
added constructor-time validation in `CoilgunOptimizationProblem`. The
validation rejects every `MaximumTemperature` constraint when
`Config::enable_thermal` is false, while allowing the constraint when thermal
simulation is enabled and preserving non-temperature constraints in the
disabled-thermal mode.

## TDD RED evidence

The new focused tests were added before the production change and the focused
test executable was rebuilt and run:

```text
cmake --build --preset cpu-debug --target test_coilgun_optimization -j2
ctest --preset cpu-debug -R '^test_coilgun_optimization$' --output-on-failure
```

The expected constructor checks failed because no validation existed:

```text
2 failed, 15 passed
CHECK_THROWS_WITH_AS(...) did NOT throw at all!
```

The failures were specifically the disabled-thermal maximum-temperature tests;
the enabled-thermal and non-temperature acceptance tests passed in this RED
run.

## GREEN implementation

The constructor now checks each constraint entry and throws:

```text
maximum-temperature constraints require thermal simulation to be enabled
```

The check runs during configuration construction, before candidate evaluation,
and does not alter thermal configuration or reject other metrics.

## Final GREEN verification

Focused adapter test:

```text
ctest --preset cpu-debug -R '^test_coilgun_optimization$' --output-on-failure
1/1 Test #27: test_coilgun_optimization ........   Passed
100% tests passed out of 1
```

Optimization-named CTest set:

```text
ctest --preset cpu-debug -R optimization --output-on-failure
100% tests passed out of 12
Total Test time (real) = 0.34 sec
```

No project lint or typecheck command is present in the worktree.

## Concerns

None. The final diff is limited to the adapter constructor, focused adapter
tests, and this report.
