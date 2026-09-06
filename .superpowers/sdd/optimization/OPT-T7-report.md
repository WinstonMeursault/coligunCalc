# OPT-T7 Report

## RED

Added `tests/test_optimization_nsga2.cpp` and registered the focused target.
The first CPU debug build failed because `coilgun/optimization/nsga2.hpp` did
not exist.

## GREEN

Implemented physics-independent NSGA-II ranking and selection in
`nsga2.hpp`/`nsga2.cpp`:

- constraint-domination gives feasible candidates precedence and orders
  infeasible candidates by hard normalized violation;
- fixed objective vectors of two or more objectives support maximize/minimize
  directions and optional objective definitions/scales;
- stable non-dominated fronts and candidate-indexed ranks are returned;
- crowding distance handles boundary points, repeated values, and degenerate
  objective ranges without NaNs;
- parent and offspring populations are merged and truncated by rank then
  descending crowding distance with stable tie order.

## Verification

`cmake --preset cpu-debug`

`cmake --build --preset cpu-debug --target test_optimization_nsga2`

`./build/cpu-debug/tests/test_optimization_nsga2`

Result: 6 test cases passed, 29 assertions passed.

`ctest --preset cpu-debug -R 'test_optimization_(nsga2|types|constraints|variables|operators|evaluator)' --output-on-failure`

Result: 6/6 focused optimization tests passed.
