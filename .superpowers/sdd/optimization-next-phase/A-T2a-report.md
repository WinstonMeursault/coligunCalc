# A-T2a — Advance the optimization random stream

## Scope

`Population::initialize` now consumes the caller-owned `RandomContext` by reference. Production and test call sites pass lvalue contexts; no by-value overload was added.

## TDD evidence

### RED

Added `population initialization advances the caller random stream` to `tests/test_optimization_operators.cpp` before changing production code. Against the original by-value signature, the focused test failed as intended:

```text
CHECK( next_draw == reference.uniform() ) is NOT correct!
values: CHECK( 0.412857 == 0.942841 )
CHECK( next_draw != fresh.uniform() ) is NOT correct!
values: CHECK( 0.412857 != 0.412857 )
```

The failure demonstrates that initialization consumed a copy and left the caller at the fresh first draw.

### GREEN

Changed the declaration and definition to accept `RandomContext&`. The focused test then passed, proving that the caller's next draw matches an independently advanced reference stream and differs from a fresh stream. The reproducibility test also uses two independent lvalue contexts with the same seed and checks their subsequent draws remain identical.

## Verification

- `cmake --build --preset cpu-debug --target test_optimization_operators -j2`
- `./build/cpu-debug/tests/test_optimization_operators` — 6 test cases, 73 assertions passed.
- `cmake --build --preset cpu-debug -j2`
- `ctest --preset cpu-debug -R 'optimization|coilgun_optimization' --output-on-failure` — 12/12 tests passed.

## Concerns

No known concerns within the requested scope. The CTest project does not define an `optimization` label; the optimization-named tests (including `test_coilgun_optimization` and integration) were selected by name.
