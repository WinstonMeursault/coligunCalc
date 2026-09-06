# OPT-T9 Report: Add coilgun optimization adapter

## Scope

Added `CoilgunOptimizationProblem` as the physics-facing adapter between
generic optimization variables and the existing CPU `MultiStageSim` API.
Bindings cover coil geometry/turns/position, excitation voltage/capacitance,
trigger values, and armature position/velocity. The adapter reports terminal
velocity as the default maximize objective and exposes maximum temperature,
peak current, peak voltage, efficiency, and energy loss as named metadata and
optional metric constraints.

The adapter implements both the context-aware `Evaluator` and `BatchEvaluator`
interfaces, so it can be passed directly to `GeneticOptimizer` while preserving
the injectable GPU callback and CPU fallback path. Peak voltage includes the
absolute initial excitation voltages as well as recorded history.

## TDD Evidence

### RED

Extended `tests/test_coilgun_optimization.cpp` with regression coverage for
direct `GeneticOptimizer` construction and initial peak voltage. Before the fix,
the callback count remained zero and a zero-step simulation reported 0 V.

### GREEN

`CoilgunOptimizationProblem` now provides the required non-const
`BatchEvaluator` override plus a const forwarding overload. `GeneticOptimizer`
has a dedicated overload that selects the batch path without base-class
overload ambiguity.

## Verification

```text
cmake --build --preset cpu-debug --target test_coilgun_optimization -j2
./build/cpu-debug/tests/test_coilgun_optimization
```

Result: 8 test cases passed, 25 assertions passed.

```text
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug --output-on-failure
```

Result: 30/30 CPU tests passed, including the adapter suite.
