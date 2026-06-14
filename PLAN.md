# Refactor Plan

This plan prioritizes stabilization, build hygiene, and clearer boundaries between the C++ core and Python bindings. The goal is to improve maintainability without changing the public Python API unless explicitly planned and tested.

## Guiding Principles

- Preserve existing behavior first; refactor behind tests.
- Keep the C++ core as the source of reusable algorithms.
- Keep Python bindings thin: adapt inputs, call core code, expose results.
- Prefer small, reviewable changes over broad rewrites.
- Modernize build and packaging incrementally.

## Phase 1: Correctness And Packaging Drift

### Fix Known Defects

- Done: fix `include/dionysus/matrix-filtration.h`: the constructor used `m_->size()` even though `m_` is stored by value.
- Done: fix `bindings/python/persistence.cpp`: reduced-matrix pickle restore checked for a one-element tuple but reads four tuple elements.
- Done: fix `CMakeLists.txt`: the `if(debug)` block referenced an undefined option.

### Fix Runtime Metadata

- Done: add `sortedcontainers` to `pyproject.toml`; `bindings/python/dionysus/zigzag.py` imports `SortedSet` directly.
- Done: review classifiers in `pyproject.toml` so they match `requires-python >=3.10` and the wheel build script.
- Done: update README/docs that still imply pybind11 is vendored, since the current build uses `find_package(pybind11 CONFIG REQUIRED)`.

### Clean Stale Files

- Done: remove or archive legacy packaging and CI files if they are no longer used:
  - `setup.py.bak` and `setup.cfg.bak` were ignored local backup files and were removed from the working tree.
  - `.travis.yml`
- Done: expand `.gitignore` for local/generated artifacts:
  - `.venv/`
  - `.mypy_cache/`
  - `.pytest_cache/`
  - existing build/doc/cache outputs as needed

### Phase 1 Progress

Completed in the first stabilization batch:

- Fixed the known `MatrixFiltration` assertion and reduced-matrix pickle restore bugs.
- Added the explicit CMake `debug` option for the existing `backward.cpp` block.
- Added the missing `sortedcontainers` runtime dependency and lockfile entry.
- Added common local cache directories to `.gitignore`.
- Converted `tests/test-issue72.py` from top-level execution into pytest tests.
- Added reduced-matrix pickle roundtrip coverage in `tests/test_matrix_u.py`.
- Made `tests/test_issue39.py` path-safe, but continue to exclude it from routine verification because it is known to hang.

Completed in the metadata cleanup batch:

- Added Python 3.10-3.14 classifiers to match `requires-python` and `build-wheels.sh`.
- Updated README and documentation wording so pybind11 is no longer described as vendored.
- Removed stale pybind11 vendoring from `peru.yaml`; pybind11 is provided by the build environment.
- Removed stale Travis CI configuration.
- Removed ignored local `setup.py.bak` and `setup.cfg.bak` files from the working tree.

## Phase 2: Modernize CMake

### Add A Public C++ Target

Done: create a target-scoped interface library for the header-only C++ core:

```cmake
add_library(dionysus INTERFACE)
target_include_directories(dionysus INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_compile_features(dionysus INTERFACE cxx_std_14)
```

### Replace Global CMake State

- Done: replace global `include_directories(...)` with `target_include_directories(...)` in the project CMake files.
- Done: replace global `add_definitions(...)` with `target_compile_definitions(...)` on the `dionysus` interface target.
- Done: replace shared `${libraries}` plumbing with explicit target links for project examples and Python bindings.
- Done: update these files first:
  - `CMakeLists.txt`
  - `bindings/python/CMakeLists.txt`
  - `examples/CMakeLists.txt`
  - `examples/*/CMakeLists.txt`

### Phase 2 Progress

Completed in the first CMake modernization batch:

- Added the header-only `dionysus` interface target with public include directories and `cxx_std_14` compile features.
- Moved `DEBUG`, `COUNTERS`, `TRACE`, `DIONYSUS_ZIGZAG_DEBUG`, and `BACKWARD_HAS_DW=1` from global definitions onto `dionysus`.
- Linked `_dionysus` privately against `dionysus` and scoped its binding-local and Hera include directories to the extension target.
- Added a `dionysus_examples` interface target for shared example includes and links.
- Updated example executables to link `dionysus_examples` instead of the shared `${libraries}` variable.

### Installation Boundary

- Done: install/export the headers and the `dionysus` interface target for C++ consumers.
- Done: keep Python packaging through `scikit-build-core` working throughout the transition.

Completed in the install/export batch:

- Added standard CMake package config generation for `DionysusConfig.cmake` and `DionysusConfigVersion.cmake`.
- Installed `include/dionysus` and exported `Dionysus::dionysus` from `DionysusTargets.cmake`.
- Used `GNUInstallDirs` for install destinations and kept the existing Python package install layout unchanged.
- Verified a downstream CMake consumer can `find_package(Dionysus CONFIG REQUIRED)` and link `Dionysus::dionysus` from a temporary install prefix.

## Phase 3: Thin The Python Binding Layer

Several binding files contain reusable algorithmic code. Move that code into the C++ core or shared internal helpers, then leave pybind files as wrappers.

### Extract Boundary Matrix Construction

- Done: move boundary/coboundary matrix construction from `bindings/python/boundary.cpp` into a reusable C++ helper.
- Reuse the helper from vineyard construction and persistence-related paths where applicable.
- Done: keep Python functions `boundary()` and `coboundary()` unchanged.

Completed in the boundary extraction batch:

- Added `include/dionysus/boundary-matrix.h` with reusable `make_boundary_matrix_filtration()` and `make_coboundary_matrix_filtration()` helpers.
- Reduced `bindings/python/boundary.cpp` to thin wrappers that call the reusable helpers for `Filtration` and `MultiFiltration`.
- Preserved the existing default prime-3 signed coefficient behavior.
- Added Python regression coverage for boundary and coboundary construction from both `Filtration` and `MultiFiltration`.

### Extract Vineyard Linear Homotopy

- Move algorithmic logic and result structs from `bindings/python/vineyard.cpp` into core headers or a small C++ source module.
- Keep `bindings/python/vineyard.cpp` focused on:
  - pybind class/function exports
  - Python object ownership
  - conversion between Python containers and C++ chains
- Preserve existing vineyard tests while moving code.

### Extract Freudenthal Construction

- Move Freudenthal triangulation generation out of `bindings/python/freudenthal.cpp`.
- Keep NumPy dtype/shape handling in the binding layer.
- Avoid `const_cast` by constructing simplices with final data where possible.

### Extract Zigzag Helpers

- Move cone construction and zigzag diagram classification from `bindings/python/zigzag-persistence.cpp` into reusable C++ helpers.
- Keep Python callback and object exposure in the binding layer.
- Keep the Python `dionysus.fast_zigzag` wrapper behavior unchanged.

## Phase 4: Consolidate Filtration Containers

The filtration containers share substantial logic:

- `include/dionysus/filtration.h`
- `include/dionysus/multi-filtration.h`
- `include/dionysus/linked-multi-filtration.h`

Refactor toward shared internal storage policies for:

- unique vs non-unique cell lookup
- optional linked index
- checked vs unchecked lookup
- index update behavior after sort/rearrange

Do this after tests cover the public behavior of all three containers.

## Phase 5: Improve Mutation And Numeric Boundaries

### Remove Unsafe Data Mutation Patterns

- Avoid `const_cast` in `bindings/python/rips.cpp` and `bindings/python/freudenthal.cpp`.
- Prefer constructing simplices with final filtration values before insertion.
- If post-insertion data updates are truly required, add a narrow, explicit API for updating simplex data safely.

### Clarify Numeric Policy

- Python-facing simplex data currently uses `float` via `bindings/python/simplex.h`.
- Diagrams use the same type via `bindings/python/diagram.h`.
- Decide whether to keep `float` for compatibility or migrate to `double` for precision.
- If migrating, add tests and document the compatibility impact.

## Phase 6: Test Improvements

### Immediate Test Hygiene

- Rename/fix `tests/test-issue72.py`; it currently has top-level execution and no assertions.
- Fix `tests/test_issue39.py` to load data relative to the test file rather than the current working directory.
- Add pytest configuration to `pyproject.toml`.
- Add `tests/conftest.py` fixtures for common filtrations, primes, and reduction methods.

### Broaden Python Coverage

Add tests for:

- `Simplex`
- `Filtration`
- `MultiFiltration`
- `LinkedMultiFiltration`
- `fill_rips`
- `fill_freudenthal`
- cohomology persistence
- omnifield persistence
- zigzag helpers
- reduced-matrix pickling
- package import from an installed wheel

### Add C++ Tests

- Add CTest integration.
- Add small C++ tests for core headers and reduction algorithms.
- Use these tests to protect refactors that move algorithmic code out of bindings.

## Suggested Execution Order

1. Fix correctness bugs and dependency metadata.
2. Add or repair tests for the fixed behavior.
3. Modernize CMake target boundaries.
4. Extract boundary/coboundary matrix construction.
5. Extract vineyard linear homotopy from the binding file.
6. Extract Freudenthal and zigzag helpers.
7. Consolidate filtration container duplication.
8. Revisit numeric precision and public API cleanup.

## Highest-Value First Refactor

The highest-value structural refactor is extracting algorithm code out of these binding files while preserving the Python API:

- `bindings/python/vineyard.cpp`
- `bindings/python/freudenthal.cpp`
- `bindings/python/zigzag-persistence.cpp`
- `bindings/python/boundary.cpp`

This should reduce binding complexity, make the C++ core easier to test directly, and prevent future algorithms from being trapped behind Python-specific code.
