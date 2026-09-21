# Repository Guidelines

## Project Structure & Module Organization

This repository is a C++20 library for multi-stage induction coilgun simulation. Public headers live in `include/coilgun/`; implementations under `src/` mirror the `components/`, `physics/`, `simulation/`, and `optimization/` namespaces. CUDA sources are isolated in `src/cuda/`. Add doctest suites to `tests/`, developer utilities to `tools/`, and validation scripts to `scripts/`. User-facing references belong in `docs/`; keep `README.md` synchronized with `README_cn.md`, and `docs/API.md` with `docs/API_cn.md`. Treat `.references/` as local research material rather than production source.

## Build, Test, and Development Commands

- `cmake --preset cpu-debug` configures a CPU debug build in `build/cpu-debug/` and fetches missing dependencies.
- `cmake --build --preset cpu-debug` builds the library and tests.
- `ctest --preset cpu-debug` runs the complete CPU test suite.
- `ctest --preset cpu-debug -L fast` runs the quick test subset.
- `cmake --preset cuda-debug && cmake --build --preset cuda-debug` enables the optional CUDA backend; use `ctest --preset cuda-debug -L gpu` for GPU-labelled tests.
- Use `cpu-release` or `cuda-release` for optimized builds. CMake 3.20+, Ninja, OpenMP, and a C++20 compiler are required; CUDA presets also require the CUDA Toolkit.

## Coding Style & Naming Conventions

Follow the existing four-space C++ indentation and nearby brace/layout conventions. Use `snake_case` for files, functions, and variables; `PascalCase` for classes and structs; and descriptive namespaces under `coilgun`. Keep public declarations in `include/coilgun/` and implementation details in `src/`. Prefer standard-library facilities and RAII, and avoid introducing a dependency when an existing Eigen or Boost facility suffices. No repository-wide formatter is configured, so keep diffs consistent with surrounding code.

## Testing Guidelines

Tests use doctest and CTest. Name new files `tests/test_<feature>.cpp`, register them in `tests/CMakeLists.txt`, and cover normal, boundary, and numerical-tolerance cases. Apply the existing `fast`, `validation`, `slow`/`integration`, or `gpu` labels appropriately. GPU tests must skip cleanly when no compatible device is available. No fixed coverage threshold is enforced; every behavior change should include a focused regression test.

## Commit & Pull Request Guidelines

Recent history follows Conventional Commit-style subjects such as `feat(optimization): ...`, `perf(cuda): ...`, and `chore(release): ...`. Keep each commit focused and use an imperative, concise subject. Pull requests should explain intent and numerical or performance impact, list verification commands, link relevant issues, and update both language versions of affected documentation. Include benchmark evidence for performance changes and explicitly note CPU/CUDA behavior differences.
