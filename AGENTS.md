# AGENTS.md

## Project Shape
- Dionysus is a C++14 computational topology library with pybind11 Python bindings; the C++ API is mostly header-only under `include/dionysus/`.
- Python packaging uses `scikit-build-core` from `pyproject.toml`; `bindings/python/CMakeLists.txt` builds the `_dionysus` extension and installs pure Python files from `bindings/python/dionysus/`.
- `src/` only contains optional `backward.cpp` for CMake `debug=ON`; do not look there for normal implementation code.
- `ext/hera` is a pinned vendored subset for bottleneck/Wasserstein distances; refresh it through `peru.yaml` rather than ad-hoc edits.

## Commands
- Default Python tests: `uv run --group test pytest -q`.
- Focused Python test: `uv run --group test pytest -q tests/test_simplex.py::test_simplex_construction_sorts_vertices_and_exposes_sequence_api`.
- Opt-in wheel smoke test: `DIONYSUS_RUN_PACKAGING_TESTS=1 uv run --group test pytest -q -m packaging`.
- Fast C++ check for header-only changes: `cmake -S . -B /tmp/dionysus-cmake -Dbuild_python_bindings=OFF -Dbuild_examples=OFF && cmake --build /tmp/dionysus-cmake && ctest --test-dir /tmp/dionysus-cmake --output-on-failure`.
- Full CMake defaults build examples and Python bindings: `cmake -S . -B /tmp/dionysus-cmake-full && cmake --build /tmp/dionysus-cmake-full`.
- Docs doctests: `uv run --group doc --group test sphinx-build -M doctest doc /tmp/dionysus-doc`; non-interactive matplotlib warnings are expected.
- `./build-wheels.sh` is a release helper: it deletes `.venv`, builds CPython 3.10-3.14 wheels locally and via `orb`, then runs `auditwheel`; do not use it for routine verification.

## Test Notes
- Default pytest currently has two expected skips: `tests/test_issue39.py` is marked skipped because the Wasserstein regression hangs, and packaging tests skip unless `DIONYSUS_RUN_PACKAGING_TESTS=1`.
- Python tests import `dionysus`; when using a manual CMake build instead of `uv`, set `PYTHONPATH=<cmake-build>/bindings/python` as in `README.rst` and `.build.yml`.
- C++ tests are registered through CTest from `tests/cpp`, with target/test name `dionysus-core-tests`.

## Gotchas
- There is no repo-level linter, formatter, or typecheck config; do not invent a lint step unless you add the config.
- The project version is duplicated in `pyproject.toml` and root `CMakeLists.txt`; keep both in sync.
- `build/`, `dist/`, `.venv/`, `doc/_build/`, and `dionysus.egg-info/` are ignored generated artifacts; prefer temp build dirs so CMake and scikit-build outputs do not mix.
