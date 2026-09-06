# OPT-T9 Report: Add coilgun optimization adapter

## Scope

Added `CoilgunOptimizationProblem` as the physics-facing adapter between
generic optimization variables and the existing CPU `MultiStageSim` API.
Bindings cover coil geometry/turns/position, excitation voltage/capacitance,
trigger values, and armature position/velocity. The adapter reports terminal
velocity as the default maximize objective and exposes maximum temperature,
peak current, peak voltage, efficiency, and energy loss as named metadata and
optional metric constraints.

The adapter also implements the context-aware `Evaluator` interface and has an
injectable GPU batch callback. A callback that throws or returns malformed
batch output falls back to isolated CPU evaluations; non-finite successful GPU
values are converted to per-candidate failures.

## TDD Evidence

### RED

Added `tests/test_coilgun_optimization.cpp` and registered it in
`tests/CMakeLists.txt`. The focused target initially failed because
`coilgun/optimization/coilgun_problem.hpp` did not exist.

### GREEN

Implemented the adapter and registered `src/optimization/coilgun_problem.cpp`
in the library target. The focused test now covers variable decoding and
metrics, malformed candidates, GPU callback fallback, and non-finite GPU
results.

## Verification

```text
cmake --build --preset cpu-debug --target test_coilgun_optimization -j2
./build/cpu-debug/tests/test_coilgun_optimization
```

Result: 3 test cases passed, 16 assertions passed.

```text
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug --output-on-failure
```

Result: 30/30 CPU tests passed, including the new adapter suite.
