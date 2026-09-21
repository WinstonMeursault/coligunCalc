# Optimization execution ledger

Source plan: `docs/superpowers/plans/2026-09-06-optimization-algorithm-execution-plan.md`

This ledger records the durable OPT-T1…OPT-T12 history from the existing reports
and commits. It is separate from the next-phase A/B ledger; B-T4 benchmarking
and whole-branch completion are intentionally not claimed here.

| Task | Status | Evidence |
|---|---|---|
| OPT-T1 domain types | complete | `0acf6a3 Add optimization domain types`; subsequent type reports |
| OPT-T2 variables/repair | complete | `1a6e585 Implement variable encoding and repair`, `3481491 Fix optimization variable validation`, `OPT-T2-report.md` |
| OPT-T3 objectives/constraints | complete | `af1fb2d Implement objectives and constraints`, `591f504 Align optimization constraint semantics`, `OPT-T3-report.md` |
| OPT-T4 batch evaluator/cache/statistics | complete | `0314937 Add batch evaluation interfaces`, `OPT-T4-report.md`, A-T4c ownership follow-up |
| OPT-T5 genetic operators | complete | `2e018c0 Implement genetic operators`, `OPT-T5-report.md` |
| OPT-T6 scalar GA | complete | `e59bf76 Add single objective optimization`, `OPT-T6-report.md` |
| OPT-T7 NSGA-II | complete | `6de12f7 Add NSGA-II selection`, `OPT-T7-report.md` |
| OPT-T8 automatic routing | complete | `OPT-T8-report.md`, later A-T3c termination/configuration fixes |
| OPT-T9 coilgun adapter | complete | `131f549 feat(optimization): implement coilgun optimization framework`, `OPT-T9-report.md` |
| OPT-T10 Pareto selectors | complete | `f69e81c Harden optimization result selectors`, `OPT-T10-report.md` |
| OPT-T11 public API/integration | complete | `include/coilgun/coilgun.hpp`, `CMakeLists.txt` install/export rules, bilingual API docs, public API smoke test, and the `Document optimization public API` task commit |
| OPT-T12 workflow validation | complete | `OPT-T12-report.md`, `tests/test_optimization_integration.cpp`, historical benchmark record |

The current implementation additionally includes the post-plan fixed-geometry
CUDA batch evaluator and run-local GPU metrics from B-T2/B-T3. The historical
2026-09-08 benchmark remains raw evidence of its measurement revision; its
“no concrete CUDA evaluator” statement is historical and does not describe the
current adapter.

## Deferred and follow-up work

- Resolve the B-T1 `PeakCurrent` discrepancy before using that metric as CUDA
  production fitness; no numerical tolerance or physics was changed here.
- Extend CUDA evaluation to arbitrary-geometry candidates and to thermal
  optimization metrics after their numerical contracts are established.
- B-T4 throughput benchmarking has now been rerun in W5 with fresh real-device
  proof and preserved historical baselines; batch 32/128 gates passed. The
  full W5 validation report is in
  `.superpowers/sdd/optimization-next-phase/final-validation-report.md`.
- Future algorithm families remain out of scope: DE, CMA-ES, Bayesian
  optimization, and MOEA-D (plus distributed evaluation and Pareto tooling).
