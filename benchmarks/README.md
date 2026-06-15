# Refactor Performance Benchmarks

`refactor_performance.py` builds selected commits in isolated git worktrees and virtual environments, then runs the same deterministic Python workloads against each installed package.

Default comparison:

```bash
python benchmarks/refactor_performance.py run --json-output /tmp/dionysus-refactor-results.json
```

The default commits are `acd84eea`, `001f9523`, and `447a64a0`. The report prints median/min/IQR timing triples in milliseconds and ratio columns normalized to the first successful commit in the requested order. With the defaults, ratio columns are `acd84eea/acd84eea`, `001f9523/acd84eea`, and `447a64a0/acd84eea`; values above `1.0x` are slower than `acd84eea`.

Useful options:

```bash
python benchmarks/refactor_performance.py list
python benchmarks/refactor_performance.py run --size quick --repeats 3 --min-time 0.01
python benchmarks/refactor_performance.py run --bench 'homology_*' --bench '*duplicate*'
python benchmarks/refactor_performance.py run --rebuild
```

For Python-level profiling, add `--cprofile` with a benchmark name, shell-style pattern, or `all`:

```bash
python benchmarks/refactor_performance.py run --size quick --cprofile homology_matrix_v
```

`cProfile` will show C++ extension calls as Python leaves. For C++ attribution, first build the target environments and print worker commands:

```bash
python benchmarks/refactor_performance.py run --build-only --print-worker-commands --commits 447a64a0
```

Then wrap the printed worker command with a native sampler, for example on macOS:

```bash
xcrun xctrace record --template 'Time Profiler' --launch -- /path/to/venv/bin/python /path/to/benchmarks/refactor_performance.py worker --bench homology_matrix_v --repeats 1 --warmups 0 --min-time 10
```

The build uses CMake `Release` and adds `-g` by default, which keeps optimized code while preserving symbols for native profilers. Use `--no-debug-symbols` to omit debug symbols.
