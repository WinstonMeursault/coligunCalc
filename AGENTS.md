# AGENTS.md — coligunCalc

## What this is
Multi-stage synchronous induction coilgun simulation using the current filament method. GPLv3.

## State: C++ physics foundation complete
The Python prototype (`basics.py`, `core.py`, `__init__.py`) has been deleted from the working tree but lives in git history (`git show 563c1d1:basics.py`). The C++ rewrite of the physics foundation layer is complete on branch `refactor/cpp` — all 13 planned commits implemented, 18 test suites active (16 passing, 1 GPU precision, 1 GPU WIP).

## Project structure
```
include/coilgun/   — Public headers (core types, physics, components)
src/               — Library implementation (static lib coilgun)
tests/             — Unit and integration tests (doctest, 9 suites)
tools/             — T(q,p) lookup table generator
docs/              — API.md / API_cn.md (C++ API), NumericalModel.md (physics)
.references/       — PDF papers; gitignored, local-only
```

## Build
```sh
cmake --preset ninja-debug && cmake --build --preset ninja-debug && ctest --preset debug
```
Dependencies (Boost.Math, Eigen, doctest) are fetched automatically via CMake FetchContent. Requires C++17 compiler and CMake ≥ 3.20.

## Key documents
- `docs/NumericalModel.md` — canonical physics reference: filament discretization, elliptic integral kernel, circuit equations, forward Euler solver
- `docs/API.md` — C++ function/class API (English)
- `docs/API_cn.md` — C++ function/class API (Chinese)

## Documentation sync rule
**API.md and API_cn.md must always be kept in sync.** Any change to one must be mirrored in the other. Same applies to README.md ↔ README_cn.md. Both language versions must contain equivalent content.

## MinerU OCR pipeline for .references/
When OCR-ing PDFs in `.references/` via local MinerU API (`http://127.0.0.1:8000`):
- **Backend**: Always `hybrid-engine` (CUDA), not `pipeline`
- **Output**: one folder per PDF under `.references/<basename>/`:
  - `md_content` → `<basename>.md`
  - `middle_json` → `middle_json.json`
  - `content_list` → `content_list.json`
  - OCR images → `images/`
- **No root pollution**: move any auto-created `output/` into the target ref folder and delete from root
