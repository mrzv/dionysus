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
- Done: reuse shared boundary construction from vineyard construction paths where applicable.
- Done: keep Python functions `boundary()` and `coboundary()` unchanged.

Completed in the boundary extraction batch:

- Added `include/dionysus/boundary-matrix.h` with reusable `make_boundary_matrix_filtration()` and `make_coboundary_matrix_filtration()` helpers.
- Reduced `bindings/python/boundary.cpp` to thin wrappers that call the reusable helpers for `Filtration` and `MultiFiltration`.
- Preserved the existing default prime-3 signed coefficient behavior.
- Added Python regression coverage for boundary and coboundary construction from both `Filtration` and `MultiFiltration`.

Completed in the vineyard boundary-chain reuse batch:

- Added reusable `make_boundary_chains()` construction to `include/dionysus/boundary-matrix.h`.
- Reduced `VineyardV`, `VineyardU`, and `Vineyard(...)` filtration constructors to call the shared helper instead of a binding-local boundary builder.
- Kept invalid filtration order validation behavior through standard C++ exceptions mapped by pybind.

### Extract Vineyard Linear Homotopy

- Done: move algorithmic logic and result structs from `bindings/python/vineyard.cpp` into core headers or a small C++ source module.
- Keep `bindings/python/vineyard.cpp` focused on:
  - pybind class/function exports
  - Python object ownership
  - conversion between Python containers and C++ chains
- Preserve existing vineyard tests while moving code.

Completed in the vineyard homotopy preparation batch:

- Added `include/dionysus/vineyard-linear-homotopy.h` with reusable endpoint validation, stable-order construction, value interpolation, slope calculation, and stable boundary-chain remapping helpers.
- Reduced `bindings/python/vineyard.cpp` by delegating linear-homotopy data preparation to the reusable helper.
- Kept event scheduling, result structs, and pybind-owned vineyard result conversion in the binding layer for follow-up vineyard extraction batches.

Completed in the vineyard scheduler extraction batch:

- Consolidated vineyard linear-homotopy tolerance on `dionysus::vineyard_linear_homotopy_epsilon`.
- Moved crossing candidate ordering, adjacent inversion checks, crossing-time calculation, event queue construction, neighborhood updates, and stale candidate popping into `include/dionysus/vineyard-linear-homotopy.h`.
- Reduced `bindings/python/vineyard.cpp` to call the reusable scheduler helpers while keeping event recording and Python result ownership in the binding layer.

Completed in the vineyard result extraction batch:

- Moved linear-homotopy event, segment, vine, active-vine, and closed-vine structs into `include/dionysus/vineyard-linear-homotopy.h`.
- Moved feature discovery, segment open/close/reopen bookkeeping, persistence-point matching, transposition validation, and transposition event recording into reusable helpers.
- Kept the Python-visible `VineyardLinearHomotopyResult` wrapper in `bindings/python/vineyard.cpp` because it owns the pybind-managed vineyard object.

### Extract Freudenthal Construction

- Done: move Freudenthal triangulation generation out of `bindings/python/freudenthal.cpp`.
- Done: keep NumPy dtype/shape handling in the binding layer.
- Done: avoid `const_cast` by constructing simplices with final data where possible.

Completed in the Freudenthal extraction batch:

- Added `include/dionysus/freudenthal.h` with a reusable `fill_freudenthal()` helper that accepts shape, element strides, and a data pointer.
- Reduced `bindings/python/freudenthal.cpp` to dtype dispatch and NumPy shape/stride adaptation.
- Preserved lower-star and upper-star behavior while assigning simplex data at construction time instead of mutating through `const_cast`.
- Added Python regression coverage for 1D lower-star, 1D upper-star, 2D lower-star, and unsupported dtype errors.

### Extract Zigzag Helpers

- Done: move cone construction from `bindings/python/zigzag-persistence.cpp` into a reusable C++ helper.
- Done: move zigzag diagram classification from `bindings/python/zigzag-persistence.cpp` into reusable C++ helpers.
- Keep Python callback and object exposure in the binding layer.
- Done: keep the Python `dionysus.fast_zigzag` wrapper behavior unchanged.

Completed in the zigzag cone extraction batch:

- Added `include/dionysus/zigzag-cone.h` with reusable `make_zigzag_cone()` construction.
- Reduced the C++ `fast_zigzag()` binding implementation to a thin wrapper that supplies the existing base and cone ordering comparators.
- Added Python regression coverage for cone ordering and the existing matrix-V issue72 path.

Completed in the zigzag diagram extraction batch:

- Moved `init_zigzag_diagrams()` classification into the reusable `include/dionysus/zigzag-cone.h` helper header.
- Kept pybind overload registration and Python object exposure in `bindings/python/zigzag-persistence.cpp`.
- Added Python regression coverage for classified closed-closed, closed-open diagonal, and open-closed intervals.

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

Completed in the filtration coverage batch:

- Added Python regression coverage for `Filtration`, `MultiFiltration`, and `LinkedMultiFiltration` construction, sorting, rearranging, indexing, and containment behavior.
- Documented current duplicate handling behavior: `Filtration` keeps unique simplices, while `MultiFiltration` and `LinkedMultiFiltration` preserve duplicate simplex vertex sets.
- Captured the current duplicate lookup contract where `index(simplex, bound)` returns the duplicate at or before the search bound and simplex `data` is not part of lookup identity.
- Added boundary-face lookup coverage for duplicate-preserving filtrations to ensure simplex boundary searches select the latest matching duplicate before the coface.

Completed in the first filtration consolidation batch:

- Added `include/dionysus/multi-filtration-common.h` for shared duplicate-preserving filtration internals.
- Reused shared indexed-cell ordering in `MultiFiltration` and `LinkedMultiFiltration`.
- Reused a shared cell projection functor for order iterators in both duplicate-preserving filtration containers.
- Reused shared rearrangement reference and reverse-index helpers for duplicate-preserving filtration reordering.
- Reused shared old-to-new index mapping and update traversal helpers while keeping plain and linked index mutations container-specific.
- Kept storage, linked-index maintenance, and public lookup behavior unchanged for later, more invasive consolidation steps.

## Phase 5: Improve Mutation And Numeric Boundaries

### Remove Unsafe Data Mutation Patterns

- Done: avoid `const_cast` in `bindings/python/rips.cpp` and `bindings/python/freudenthal.cpp`.
- Done: prefer constructing simplices with final filtration values before insertion.
- If post-insertion data updates are truly required, add a narrow, explicit API for updating simplex data safely.

Completed in the Rips mutation cleanup batch:

- Updated `bindings/python/rips.cpp` so generated Rips simplices receive their final filtration value before insertion.
- Preserved point-cloud radius checks with squared distances while storing public simplex data as unsquared distances.
- Added Python regression coverage for point-cloud Rips, explicit-distance Rips, and unsupported dtype errors.

### Clarify Numeric Policy

- Python-facing simplex data currently uses `float` via `bindings/python/simplex.h`.
- Diagrams use the same type via `bindings/python/diagram.h`.
- Done: keep `float` for compatibility and document the current precision policy.
- Done: add tests that pin simplex data and diagram coordinates to the current single-precision behavior.

Completed in the numeric policy batch:

- Documented in `doc/API.rst` that Python-facing simplex data and diagram coordinates use C++ `float` precision.
- Added regression coverage for float32 rounding of simplex data and persistence diagram coordinates.

## Phase 6: Test Improvements

### Immediate Test Hygiene

- Done: rename/fix `tests/test-issue72.py`; it currently has top-level execution and no assertions.
- Done: fix `tests/test_issue39.py` to load data relative to the test file rather than the current working directory.
- Done: add pytest configuration to `pyproject.toml`.
- Done: add `tests/conftest.py` fixtures for common filtrations, primes, and reduction methods.

Completed in the test hygiene batch:

- Renamed `tests/test-issue72.py` to `tests/test_issue72.py` so it matches standard pytest discovery patterns.
- Added `tool.pytest.ini_options.testpaths = ["tests"]` to `pyproject.toml`.
- Marked `tests/test_issue39.py` skipped with a documented reason because the Wasserstein regression is still known to hang.

Completed in the shared fixture batch:

- Added `tests/conftest.py` with shared fixtures for the common prime, matrix reduction methods, and triangle filtration cells.
- Migrated boundary, matrix-U, and vineyard linear-homotopy tests to use the shared fixtures where that removed local duplication.

### Broaden Python Coverage

Add tests for:

- Done: `Simplex`
- `Filtration`
- `MultiFiltration`
- `LinkedMultiFiltration`
- `fill_rips`
- `fill_freudenthal`
- Done: cohomology persistence
- Done: omnifield persistence
- zigzag helpers
- reduced-matrix pickling
- Done: package import from an installed wheel

Completed in the installed-wheel import batch:

- Added an opt-in packaging regression that builds a wheel, installs it into an isolated virtualenv, and imports `dionysus` from outside the source tree.
- Registered the `packaging` pytest marker and kept the installed-wheel test skipped by default unless `DIONYSUS_RUN_PACKAGING_TESTS=1` is set.

Completed in the persistence variant coverage batch:

- Added cohomology persistence coverage that compares pairs and diagrams against ordinary homology persistence on a small filtration.
- Added cocycle access coverage for alive cohomology representatives when `keep_cocycles=True`.
- Added omnifield persistence coverage for basic prime-specific columns and diagram initialization.

Completed in the simplex coverage batch:

- Added focused `Simplex` coverage for construction, sorted vertex iteration, sequence access, containment, boundary order, joining, mutable data, identity semantics, hashing, and index errors.

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
