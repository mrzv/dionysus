#!/usr/bin/env python3
"""Benchmark Python-facing Dionysus workloads across git commits.

The script has two modes:

* ``run`` creates one detached worktree and one virtual environment per commit,
  installs that commit, then invokes the same worker workload in each venv.
* ``worker`` imports the installed ``dionysus`` package in the current Python
  environment and emits JSON timings. This mode is also useful under native
  profilers such as xctrace/Instruments or perf.
"""

from __future__ import annotations

import argparse
import cProfile
import fnmatch
import gc
import io
import json
import math
import os
import platform
import pstats
import random
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence


DEFAULT_COMMITS = ("acd84eea", "001f9523", "447a64a0")

BENCHMARK_DESCRIPTIONS = (
    ("freudenthal_2d_build", "build a 2D lower-star Freudenthal filtration"),
    ("freudenthal_3d_build", "build a 3D lower-star Freudenthal filtration"),
    ("rips_points_build", "build a Rips filtration from point coordinates"),
    ("rips_distances_build", "build a Rips filtration from lower-triangular distances"),
    ("boundary_freudenthal_2d", "materialize a boundary matrix from a Freudenthal filtration"),
    ("coboundary_rips_points", "materialize a coboundary matrix from a Rips filtration"),
    ("homology_clearing", "run homology_persistence(method='clearing')"),
    ("homology_column", "run homology_persistence(method='column')"),
    ("homology_matrix_v", "run homology_persistence(method='matrix_v')"),
    ("homology_matrix_u", "run homology_persistence(method='matrix_u')"),
    ("fast_zigzag_build", "build the fast zigzag cone"),
    ("zigzag_homology_matrix_v", "run matrix_v persistence on the fast zigzag cone"),
    ("vineyard_linear_homotopy_matrix_v", "run vineyard_linear_homotopy(method='matrix_v')"),
    ("vineyard_linear_homotopy_matrix_u", "run vineyard_linear_homotopy(method='matrix_u')"),
    ("multi_filtration_sort_duplicates", "construct and sort duplicate-heavy MultiFiltration"),
    ("multi_filtration_index_duplicates", "perform duplicate boundary lookups in MultiFiltration"),
    ("linked_multi_filtration_sort_duplicates", "construct and sort duplicate-heavy LinkedMultiFiltration"),
    ("linked_multi_filtration_index_duplicates", "perform duplicate boundary lookups in LinkedMultiFiltration"),
)

BENCHMARK_ORDER = tuple(name for name, _ in BENCHMARK_DESCRIPTIONS)

SIZE_PRESETS: dict[str, dict[str, Any]] = {
    "quick": {
        "freudenthal_2d_shape": (18, 18),
        "freudenthal_3d_shape": (6, 6, 6),
        "rips_points": 24,
        "rips_radius": 0.34,
        "rips_skeleton": 2,
        "zigzag_vertices": 12,
        "vineyard_points": 9,
        "duplicate_vertices": 24,
        "duplicate_repeats": 6,
        "duplicate_queries": 300,
    },
    "default": {
        "freudenthal_2d_shape": (64, 64),
        "freudenthal_3d_shape": (14, 14, 14),
        "rips_points": 80,
        "rips_radius": 0.32,
        "rips_skeleton": 2,
        "zigzag_vertices": 48,
        "vineyard_points": 22,
        "duplicate_vertices": 160,
        "duplicate_repeats": 14,
        "duplicate_queries": 6000,
    },
    "large": {
        "freudenthal_2d_shape": (110, 110),
        "freudenthal_3d_shape": (20, 20, 16),
        "rips_points": 130,
        "rips_radius": 0.30,
        "rips_skeleton": 2,
        "zigzag_vertices": 80,
        "vineyard_points": 30,
        "duplicate_vertices": 300,
        "duplicate_repeats": 18,
        "duplicate_queries": 15000,
    },
}

BLACKHOLE = 0


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    description: str
    run: Callable[[], Any]


@dataclass(frozen=True)
class Target:
    ref: str
    sha: str
    short_sha: str
    label: str
    path_label: str
    worktree: Path
    venv: Path

    @property
    def python(self) -> Path:
        if os.name == "nt":
            return self.venv / "Scripts" / "python.exe"
        return self.venv / "bin" / "python"


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "list":
        print_benchmark_list()
        return 0
    if args.command == "worker":
        return worker_main(args)
    return controller_main(args)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build selected Dionysus commits and benchmark Python-facing workloads."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="build commits and run benchmarks")
    add_shared_benchmark_args(run_parser)
    run_parser.add_argument(
        "--commits",
        nargs="+",
        default=list(DEFAULT_COMMITS),
        help="git refs to compare, in display order",
    )
    run_parser.add_argument(
        "--repo",
        type=Path,
        default=Path.cwd(),
        help="source repository containing the target commits",
    )
    run_parser.add_argument(
        "--work-root",
        type=Path,
        default=Path(tempfile.gettempdir()) / "dionysus-refactor-performance",
        help="directory for generated worktrees, venvs, profiles, and optional output",
    )
    run_parser.add_argument(
        "--python",
        type=Path,
        default=Path(sys.executable),
        help="Python executable used to create per-commit virtual environments",
    )
    run_parser.add_argument(
        "--rebuild",
        action="store_true",
        help="discard and recreate existing per-commit virtual environments",
    )
    run_parser.add_argument(
        "--skip-install",
        action="store_true",
        help="reuse existing virtual environments without installing packages",
    )
    run_parser.add_argument(
        "--build-only",
        action="store_true",
        help="prepare worktrees and venvs, then print worker commands without running benchmarks",
    )
    run_parser.add_argument(
        "--no-debug-symbols",
        dest="debug_symbols",
        action="store_false",
        help="do not add -g to C/C++ compiler flags during Release builds",
    )
    run_parser.add_argument(
        "--pip-install-arg",
        action="append",
        default=[],
        help="extra argument appended to each pip install command; repeatable",
    )
    run_parser.add_argument(
        "--json-output",
        type=Path,
        help="write full machine-readable results to this path",
    )
    run_parser.add_argument(
        "--keep-going",
        action="store_true",
        help="continue with later commits if one commit fails",
    )
    run_parser.add_argument(
        "--print-worker-commands",
        action="store_true",
        help="print per-commit worker commands for manual native profiling",
    )
    run_parser.set_defaults(command="run", debug_symbols=True)

    worker_parser = subparsers.add_parser(
        "worker",
        help="run workloads against the dionysus package installed in this Python environment",
    )
    add_shared_benchmark_args(worker_parser)
    worker_parser.add_argument(
        "--profile-output-dir",
        type=Path,
        help="directory for cProfile .pstats and text summaries",
    )
    worker_parser.set_defaults(command="worker")

    list_parser = subparsers.add_parser("list", help="list available benchmarks")
    list_parser.set_defaults(command="list")

    return parser


def add_shared_benchmark_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--size",
        choices=tuple(SIZE_PRESETS),
        default="default",
        help="input-size preset",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=7,
        help="timed samples per benchmark",
    )
    parser.add_argument(
        "--warmups",
        type=int,
        default=1,
        help="untimed warmup iterations per benchmark",
    )
    parser.add_argument(
        "--min-time",
        type=float,
        default=0.05,
        help="minimum seconds per timed sample; fast benchmarks are looped",
    )
    parser.add_argument(
        "--bench",
        action="append",
        default=[],
        help="benchmark name or shell-style pattern; repeatable; defaults to all",
    )
    parser.add_argument(
        "--cprofile",
        action="append",
        default=[],
        help="also cProfile benchmarks matching this name or pattern; repeatable; use 'all' for selected benchmarks",
    )
    parser.add_argument(
        "--profile-iterations",
        type=int,
        default=1,
        help="iterations per selected cProfile run",
    )


def print_benchmark_list() -> None:
    width = max(len(name) for name, _ in BENCHMARK_DESCRIPTIONS)
    for name, description in BENCHMARK_DESCRIPTIONS:
        print(f"{name:<{width}}  {description}")


def controller_main(args: argparse.Namespace) -> int:
    repo = resolve_repo(args.repo)
    work_root = args.work_root.expanduser().resolve()
    args.work_root = work_root
    script = Path(__file__).resolve()
    targets = make_targets(repo, work_root, args.commits)

    work_root.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []

    for target in targets:
        try:
            print(f"Preparing {target.label} ({target.short_sha})", file=sys.stderr)
            ensure_worktree(repo, target)
            if not args.skip_install:
                ensure_venv(target, args)
            elif not target.python.exists():
                raise RuntimeError(f"--skip-install requested, but {target.python} does not exist")

            worker_command = build_worker_command(script, target, args)
            if args.print_worker_commands or args.build_only:
                print(shell_join(worker_command))
            if args.build_only:
                continue

            result = run_worker(target, worker_command, work_root, args)
            results.append(result)
        except Exception as exc:  # noqa: BLE001 - controller should report commit-level failures.
            failures.append({"ref": target.ref, "sha": target.sha, "error": str(exc)})
            print(f"FAILED {target.label}: {exc}", file=sys.stderr)
            if not args.keep_going:
                raise

    if args.build_only:
        return 1 if failures else 0

    payload = build_report(repo, work_root, targets, results, failures, args)
    if args.json_output:
        output_path = args.json_output.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"Wrote JSON results to {output_path}", file=sys.stderr)

    print_human_report(payload)
    return 1 if failures else 0


def resolve_repo(path: Path) -> Path:
    completed = run_checked(
        ["git", "-C", str(path), "rev-parse", "--show-toplevel"],
        capture=True,
    )
    return Path(completed.stdout.strip()).resolve()


def make_targets(repo: Path, work_root: Path, refs: Sequence[str]) -> list[Target]:
    targets: list[Target] = []
    for index, ref in enumerate(refs):
        sha = git_stdout(repo, "rev-parse", ref)
        short_sha = git_stdout(repo, "rev-parse", "--short=8", sha)
        label = ref
        path_label = f"{index:02d}-{short_sha}-{sanitize_path_part(ref)}"
        targets.append(
            Target(
                ref=ref,
                sha=sha,
                short_sha=short_sha,
                label=label,
                path_label=path_label,
                worktree=work_root / "worktrees" / path_label,
                venv=work_root / "venvs" / path_label,
            )
        )
    return targets


def ensure_worktree(repo: Path, target: Target) -> None:
    target.worktree.parent.mkdir(parents=True, exist_ok=True)
    if target.worktree.exists():
        existing = git_stdout(target.worktree, "rev-parse", "HEAD")
        if existing != target.sha:
            raise RuntimeError(
                f"existing worktree {target.worktree} is at {existing[:12]}, expected {target.sha[:12]}"
            )
        return

    run_checked(
        ["git", "-C", str(repo), "worktree", "add", "--detach", str(target.worktree), target.sha]
    )


def ensure_venv(target: Target, args: argparse.Namespace) -> None:
    target.venv.parent.mkdir(parents=True, exist_ok=True)
    marker = target.venv / ".dionysus-refactor-performance.json"
    if target.python.exists() and marker.exists() and not args.rebuild:
        try:
            metadata = json.loads(marker.read_text())
        except json.JSONDecodeError:
            metadata = {}
        if (
            metadata.get("sha") == target.sha
            and metadata.get("debug_symbols") == args.debug_symbols
            and metadata.get("pip_install_arg") == args.pip_install_arg
        ):
            return

    if target.venv.exists():
        shutil.rmtree(target.venv)

    run_checked([str(args.python), "-m", "venv", str(target.venv)])
    run_checked([str(target.python), "-m", "pip", "install", "--upgrade", "pip"])

    env = build_install_env(args.debug_symbols)
    install_command = [
        str(target.python),
        "-m",
        "pip",
        "install",
        "--verbose",
        "--config-settings=cmake.build-type=Release",
        *args.pip_install_arg,
        str(target.worktree),
    ]
    run_checked(install_command, env=env)
    marker.write_text(
        json.dumps(
            {
                "sha": target.sha,
                "ref": target.ref,
                "debug_symbols": args.debug_symbols,
                "pip_install_arg": args.pip_install_arg,
                "installed_at": datetime.now(timezone.utc).isoformat(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def build_install_env(debug_symbols: bool) -> dict[str, str]:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["PYTHONNOUSERSITE"] = "1"
    env.setdefault("CMAKE_BUILD_TYPE", "Release")
    env.setdefault("CMAKE_BUILD_PARALLEL_LEVEL", str(os.cpu_count() or 2))
    if debug_symbols:
        env["CFLAGS"] = append_flag(env.get("CFLAGS", ""), "-g")
        env["CXXFLAGS"] = append_flag(env.get("CXXFLAGS", ""), "-g")
    return env


def append_flag(value: str, flag: str) -> str:
    parts = value.split()
    if flag not in parts:
        parts.append(flag)
    return " ".join(parts)


def build_worker_command(script: Path, target: Target, args: argparse.Namespace) -> list[str]:
    command = [
        str(target.python),
        str(script),
        "worker",
        "--size",
        args.size,
        "--repeats",
        str(args.repeats),
        "--warmups",
        str(args.warmups),
        "--min-time",
        str(args.min_time),
        "--profile-iterations",
        str(args.profile_iterations),
    ]
    for pattern in args.bench:
        command.extend(["--bench", pattern])
    for pattern in args.cprofile:
        command.extend(["--cprofile", pattern])
    if args.cprofile:
        command.extend(["--profile-output-dir", str(args.work_root / "profiles" / target.path_label)])
    return command


def run_worker(
    target: Target,
    command: Sequence[str],
    work_root: Path,
    args: argparse.Namespace,
) -> dict[str, Any]:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["PYTHONNOUSERSITE"] = "1"
    completed = run_checked(command, cwd=work_root, env=env, capture=True)
    try:
        worker_payload = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            f"worker for {target.label} did not emit JSON: {exc}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        ) from exc

    worker_payload.update(
        {
            "ref": target.ref,
            "sha": target.sha,
            "short_sha": target.short_sha,
            "label": target.label,
            "worktree": str(target.worktree),
            "venv": str(target.venv),
            "worker_command": shell_join(command),
        }
    )
    if completed.stderr:
        print(completed.stderr, file=sys.stderr, end="")
    return worker_payload


def worker_main(args: argparse.Namespace) -> int:
    validate_worker_args(args)
    print(f"Preparing {args.size} worker inputs", file=sys.stderr)

    import dionysus as d  # noqa: PLC0415 - imported only in worker mode.
    import numpy as np  # noqa: PLC0415 - imported only in worker mode.

    cases = make_benchmark_cases(d, np, args.size)
    selected = filter_cases(cases, args.bench)
    if not selected:
        raise SystemExit("no benchmarks matched --bench filters")

    results = []
    for case in selected:
        print(f"Running {case.name}", file=sys.stderr)
        try:
            results.append(run_timing(case, args.repeats, args.warmups, args.min_time))
        except Exception as exc:  # noqa: BLE001 - keep later benchmarks usable.
            print(f"FAILED {case.name}: {exc}", file=sys.stderr)
            results.append(
                {
                    "name": case.name,
                    "description": case.description,
                    "error": str(exc),
                }
            )

    profile_outputs: list[dict[str, str | int]] = []
    profile_cases = select_profile_cases(selected, args.cprofile)
    if profile_cases:
        if args.profile_output_dir is None:
            raise SystemExit("--cprofile requires --profile-output-dir in worker mode")
        args.profile_output_dir.mkdir(parents=True, exist_ok=True)
        profile_outputs = run_cprofiles(
            profile_cases,
            args.profile_output_dir,
            args.profile_iterations,
        )

    payload = {
        "schema": 1,
        "mode": "worker",
        "size": args.size,
        "size_config": json_ready(SIZE_PRESETS[args.size]),
        "repeats": args.repeats,
        "warmups": args.warmups,
        "min_time_seconds": args.min_time,
        "python_executable": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "dionysus_module": getattr(d, "__file__", None),
        "dionysus_version": getattr(d, "__version__", None),
        "benchmarks": results,
        "cprofile_outputs": profile_outputs,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


def validate_worker_args(args: argparse.Namespace) -> None:
    if args.repeats < 1:
        raise SystemExit("--repeats must be at least 1")
    if args.warmups < 0:
        raise SystemExit("--warmups must be non-negative")
    if args.min_time < 0:
        raise SystemExit("--min-time must be non-negative")
    if args.profile_iterations < 1:
        raise SystemExit("--profile-iterations must be at least 1")


def make_benchmark_cases(d: Any, np: Any, size: str) -> list[BenchmarkCase]:
    config = SIZE_PRESETS[size]

    freudenthal_2d = deterministic_array(np, config["freudenthal_2d_shape"], offset=0.17).astype(
        np.float32
    )
    freudenthal_3d = deterministic_array(np, config["freudenthal_3d_shape"], offset=0.41).astype(
        np.float64
    )
    points = deterministic_points(np, config["rips_points"])
    distances = lower_triangular_distances(np, points)

    rips_skeleton = config["rips_skeleton"]
    rips_radius = config["rips_radius"]
    freudenthal_2d_filtration = d.fill_freudenthal(freudenthal_2d)
    rips_points_filtration = d.fill_rips(points, rips_skeleton, rips_radius)

    zigzag_filtration = make_zigzag_filtration(d, config["zigzag_vertices"])
    zigzag_times = make_zigzag_times(zigzag_filtration)
    zigzag_cone = d.fast_zigzag(zigzag_filtration, zigzag_times)

    vineyard_filtration = d.Filtration(complete_simplex_skeleton(config["vineyard_points"], 2))
    vineyard_values0 = [float(index) for index in range(len(vineyard_filtration))]
    vineyard_values1 = values_for_dimension_shuffle(vineyard_filtration)

    duplicate_cells = make_duplicate_cells(
        d,
        config["duplicate_vertices"],
        config["duplicate_repeats"],
    )
    duplicate_queries = make_duplicate_queries(d, duplicate_cells, config["duplicate_queries"])
    multi_index_filtration = d.MultiFiltration(duplicate_cells)
    multi_index_filtration.sort()
    linked_index_filtration = d.LinkedMultiFiltration(
        [(cell, index) for index, cell in enumerate(duplicate_cells)]
    )
    linked_index_filtration.sort()

    def multi_sort() -> Any:
        filtration = d.MultiFiltration(duplicate_cells)
        filtration.sort()
        return filtration

    def linked_multi_sort() -> Any:
        filtration = d.LinkedMultiFiltration(
            [(cell, index) for index, cell in enumerate(duplicate_cells)]
        )
        filtration.sort()
        return filtration

    def multi_index() -> int:
        total = 0
        for face, bound in duplicate_queries:
            total += multi_index_filtration.index(face, bound)
        return total

    def linked_multi_index() -> int:
        total = 0
        for face, bound in duplicate_queries:
            total += linked_index_filtration.index(face, bound)
        return total

    return [
        BenchmarkCase(
            "freudenthal_2d_build",
            description_for("freudenthal_2d_build"),
            lambda: d.fill_freudenthal(freudenthal_2d),
        ),
        BenchmarkCase(
            "freudenthal_3d_build",
            description_for("freudenthal_3d_build"),
            lambda: d.fill_freudenthal(freudenthal_3d),
        ),
        BenchmarkCase(
            "rips_points_build",
            description_for("rips_points_build"),
            lambda: d.fill_rips(points, rips_skeleton, rips_radius),
        ),
        BenchmarkCase(
            "rips_distances_build",
            description_for("rips_distances_build"),
            lambda: d.fill_rips(distances, rips_skeleton, rips_radius),
        ),
        BenchmarkCase(
            "boundary_freudenthal_2d",
            description_for("boundary_freudenthal_2d"),
            lambda: d.boundary(freudenthal_2d_filtration),
        ),
        BenchmarkCase(
            "coboundary_rips_points",
            description_for("coboundary_rips_points"),
            lambda: d.coboundary(rips_points_filtration),
        ),
        BenchmarkCase(
            "homology_clearing",
            description_for("homology_clearing"),
            lambda: d.homology_persistence(rips_points_filtration, prime=2, method="clearing"),
        ),
        BenchmarkCase(
            "homology_column",
            description_for("homology_column"),
            lambda: d.homology_persistence(rips_points_filtration, prime=2, method="column"),
        ),
        BenchmarkCase(
            "homology_matrix_v",
            description_for("homology_matrix_v"),
            lambda: d.homology_persistence(rips_points_filtration, prime=2, method="matrix_v"),
        ),
        BenchmarkCase(
            "homology_matrix_u",
            description_for("homology_matrix_u"),
            lambda: d.homology_persistence(rips_points_filtration, prime=2, method="matrix_u"),
        ),
        BenchmarkCase(
            "fast_zigzag_build",
            description_for("fast_zigzag_build"),
            lambda: d.fast_zigzag(zigzag_filtration, zigzag_times),
        ),
        BenchmarkCase(
            "zigzag_homology_matrix_v",
            description_for("zigzag_homology_matrix_v"),
            lambda: d.homology_persistence(zigzag_cone, prime=2, method="matrix_v"),
        ),
        BenchmarkCase(
            "vineyard_linear_homotopy_matrix_v",
            description_for("vineyard_linear_homotopy_matrix_v"),
            lambda: d.vineyard_linear_homotopy(
                vineyard_filtration,
                vineyard_values0,
                vineyard_values1,
                field=d.Zp(5),
                method="matrix_v",
            ),
        ),
        BenchmarkCase(
            "vineyard_linear_homotopy_matrix_u",
            description_for("vineyard_linear_homotopy_matrix_u"),
            lambda: d.vineyard_linear_homotopy(
                vineyard_filtration,
                vineyard_values0,
                vineyard_values1,
                field=d.Zp(5),
                method="matrix_u",
            ),
        ),
        BenchmarkCase(
            "multi_filtration_sort_duplicates",
            description_for("multi_filtration_sort_duplicates"),
            multi_sort,
        ),
        BenchmarkCase(
            "multi_filtration_index_duplicates",
            description_for("multi_filtration_index_duplicates"),
            multi_index,
        ),
        BenchmarkCase(
            "linked_multi_filtration_sort_duplicates",
            description_for("linked_multi_filtration_sort_duplicates"),
            linked_multi_sort,
        ),
        BenchmarkCase(
            "linked_multi_filtration_index_duplicates",
            description_for("linked_multi_filtration_index_duplicates"),
            linked_multi_index,
        ),
    ]


def description_for(name: str) -> str:
    return dict(BENCHMARK_DESCRIPTIONS)[name]


def deterministic_array(np: Any, shape: Sequence[int], offset: float) -> Any:
    values = np.arange(math.prod(shape), dtype=np.float64).reshape(tuple(shape))
    return np.sin(values * 0.013 + offset) + np.cos(values * 0.021 + 2.0 * offset)


def deterministic_points(np: Any, count: int) -> Any:
    indices = np.arange(count, dtype=np.float64)
    x = np.mod(indices * 0.6180339887498949, 1.0)
    y = np.mod(indices * 0.4142135623730950, 1.0)
    z = np.mod(indices * 0.7320508075688772, 1.0)
    return np.column_stack((x, y, z)).astype(np.float64)


def lower_triangular_distances(np: Any, points: Any) -> Any:
    count = len(points)
    distances = np.empty(count * (count - 1) // 2, dtype=np.float64)
    cursor = 0
    for i in range(1, count):
        diff = points[i] - points[:i]
        values = np.sqrt(np.sum(diff * diff, axis=1))
        distances[cursor : cursor + i] = values
        cursor += i
    return distances


def complete_simplex_skeleton(points: int, skeleton: int) -> list[tuple[int, ...]]:
    return [
        simplex
        for dimension in range(skeleton + 1)
        for simplex in combinations(range(points), dimension + 1)
    ]


def make_zigzag_filtration(d: Any, vertices: int) -> Any:
    edges = {
        tuple(sorted((vertex, (vertex + offset) % vertices)))
        for vertex in range(vertices)
        for offset in (1, 2)
    }
    cells = [(vertex,) for vertex in range(vertices)] + sorted(edges)
    return d.Filtration(cells)


def make_zigzag_times(filtration: Any) -> list[list[float]]:
    horizon = float(len(filtration) + 12)
    times: list[list[float]] = []
    for simplex in filtration:
        if simplex.dimension() == 0:
            times.append([0.0])
            continue
        start = 1.0 + float(sum(simplex) % 11)
        times.append([start, min(horizon - 1.0, start + 9.0)])
    return times


def values_for_dimension_shuffle(filtration: Any) -> list[float]:
    rng = random.Random(20260614)
    groups: dict[int, list[int]] = {}
    for index, simplex in enumerate(filtration):
        groups.setdefault(simplex.dimension(), []).append(index)

    order: list[int] = []
    for dimension in sorted(groups):
        indices = groups[dimension]
        rng.shuffle(indices)
        order.extend(indices)

    values = [0.0] * len(order)
    for value, cell in enumerate(order):
        values[cell] = float(value)
    return values


def make_duplicate_cells(d: Any, vertices: int, repeats: int) -> list[Any]:
    cells = [
        d.Simplex([vertex], float(repeat))
        for repeat in range(repeats)
        for vertex in range(vertices)
    ]
    cells.extend(
        d.Simplex([vertex, vertex + 1], float(repeats + repeat))
        for repeat in range(repeats)
        for vertex in range(vertices - 1)
    )
    return cells


def make_duplicate_queries(d: Any, cells: Sequence[Any], count: int) -> list[tuple[Any, int]]:
    filtration = d.MultiFiltration(cells)
    filtration.sort()
    edges = [index for index, simplex in enumerate(filtration) if simplex.dimension() == 1]
    if not edges:
        return []

    queries: list[tuple[Any, int]] = []
    for query_index in range(count):
        edge_index = edges[query_index % len(edges)]
        boundary = list(filtration[edge_index].boundary())
        face = boundary[query_index % len(boundary)]
        queries.append((face, edge_index))
    return queries


def filter_cases(cases: Sequence[BenchmarkCase], patterns: Sequence[str]) -> list[BenchmarkCase]:
    expanded = expand_patterns(patterns)
    if not expanded:
        return list(cases)
    return [case for case in cases if any(fnmatch.fnmatchcase(case.name, pattern) for pattern in expanded)]


def expand_patterns(patterns: Sequence[str]) -> list[str]:
    expanded: list[str] = []
    for pattern in patterns:
        expanded.extend(part.strip() for part in pattern.split(",") if part.strip())
    return expanded


def select_profile_cases(
    selected: Sequence[BenchmarkCase],
    patterns: Sequence[str],
) -> list[BenchmarkCase]:
    expanded = expand_patterns(patterns)
    if not expanded:
        return []
    if "all" in expanded:
        return list(selected)
    return [case for case in selected if any(fnmatch.fnmatchcase(case.name, pattern) for pattern in expanded)]


def run_timing(
    case: BenchmarkCase,
    repeats: int,
    warmups: int,
    min_time: float,
) -> dict[str, Any]:
    global BLACKHOLE

    for _ in range(warmups):
        BLACKHOLE ^= consume(case.run())

    samples: list[float] = []
    loop_counts: list[int] = []
    for _ in range(repeats):
        gc.collect()
        was_enabled = gc.isenabled()
        gc.disable()
        count = 0
        started = time.perf_counter()
        try:
            while True:
                BLACKHOLE ^= consume(case.run())
                count += 1
                elapsed = time.perf_counter() - started
                if count >= 1 and elapsed >= min_time:
                    break
        finally:
            if was_enabled:
                gc.enable()
        samples.append(elapsed / count)
        loop_counts.append(count)

    return {
        "name": case.name,
        "description": case.description,
        "seconds": samples,
        "loops": loop_counts,
        "median": statistics.median(samples),
        "min": min(samples),
        "stdev": statistics.stdev(samples) if len(samples) > 1 else 0.0,
        "iqr": iqr(samples),
        "blackhole": BLACKHOLE,
    }


def consume(value: Any) -> int:
    if value is None:
        return 0
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value & 0xFFFF
    if isinstance(value, float):
        return int(value) & 0xFFFF
    if isinstance(value, (list, tuple)):
        total = len(value)
        for item in value[:4]:
            total ^= consume(item)
        return total

    total = 1
    for attr in ("events", "vines", "final_order"):
        if hasattr(value, attr):
            try:
                total ^= len(getattr(value, attr))
            except TypeError:
                pass
    try:
        total ^= len(value)
    except TypeError:
        pass
    return total


def iqr(values: Sequence[float]) -> float:
    ordered = sorted(values)
    if len(ordered) < 2:
        return 0.0
    return percentile(ordered, 0.75) - percentile(ordered, 0.25)


def percentile(ordered: Sequence[float], fraction: float) -> float:
    if len(ordered) == 1:
        return ordered[0]
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def run_cprofiles(
    cases: Sequence[BenchmarkCase],
    output_dir: Path,
    iterations: int,
) -> list[dict[str, str | int]]:
    outputs: list[dict[str, str | int]] = []
    for case in cases:
        print(f"cProfile {case.name}", file=sys.stderr)
        profiler = cProfile.Profile()
        try:
            profiler.enable()
            for _ in range(iterations):
                consume(case.run())
        except Exception as exc:  # noqa: BLE001 - profile failures should be reportable.
            profiler.disable()
            print(f"FAILED cProfile {case.name}: {exc}", file=sys.stderr)
            outputs.append({"benchmark": case.name, "error": str(exc), "iterations": iterations})
            continue
        profiler.disable()

        stats_path = output_dir / f"{case.name}.pstats"
        text_path = output_dir / f"{case.name}.txt"
        profiler.dump_stats(str(stats_path))
        stream = io.StringIO()
        pstats.Stats(profiler, stream=stream).sort_stats("cumulative").print_stats(80)
        text_path.write_text(stream.getvalue())
        outputs.append({"benchmark": case.name, "pstats": str(stats_path), "text": str(text_path)})
    return outputs


def build_report(
    repo: Path,
    work_root: Path,
    targets: Sequence[Target],
    results: Sequence[dict[str, Any]],
    failures: Sequence[dict[str, str]],
    args: argparse.Namespace,
) -> dict[str, Any]:
    result_targets = {result["sha"] for result in results}
    active_targets = [target for target in targets if target.sha in result_targets]
    comparisons = build_comparisons(active_targets, results)
    return {
        "schema": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_repo": str(repo),
        "work_root": str(work_root),
        "commits": [
            {
                "ref": target.ref,
                "sha": target.sha,
                "short_sha": target.short_sha,
                "label": target.label,
                "worktree": str(target.worktree),
                "venv": str(target.venv),
            }
            for target in active_targets
        ],
        "benchmark_config": {
            "size": args.size,
            "size_config": json_ready(SIZE_PRESETS[args.size]),
            "repeats": args.repeats,
            "warmups": args.warmups,
            "min_time_seconds": args.min_time,
            "bench_filters": args.bench,
            "cprofile_filters": args.cprofile,
            "profile_iterations": args.profile_iterations,
        },
        "results": list(results),
        "comparisons": comparisons,
        "failures": list(failures),
    }


def build_comparisons(
    targets: Sequence[Target],
    results: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    if not targets:
        return []

    result_by_sha = {result["sha"]: result for result in results}
    benchmark_names = ordered_benchmark_names(results)
    baseline = targets[0]
    comparisons: list[dict[str, Any]] = []

    for name in benchmark_names:
        ratios: dict[str, float | None] = {}
        baseline_median = benchmark_stat(result_by_sha.get(baseline.sha), name, "median")
        for target in targets:
            target_median = benchmark_stat(result_by_sha.get(target.sha), name, "median")
            ratio_name = f"{target.label}/{baseline.label}"
            ratios[ratio_name] = (
                target_median / baseline_median
                if baseline_median and target_median is not None
                else None
            )
        comparisons.append(
            {"benchmark": name, "baseline": baseline.label, "ratios": ratios}
        )
    return comparisons


def ordered_benchmark_names(results: Sequence[dict[str, Any]]) -> list[str]:
    names = {
        benchmark["name"]
        for result in results
        for benchmark in result.get("benchmarks", [])
    }
    ordered = [name for name in BENCHMARK_ORDER if name in names]
    ordered.extend(sorted(names.difference(ordered)))
    return ordered


def benchmark_stat(result: dict[str, Any] | None, name: str, stat: str) -> float | None:
    if result is None:
        return None
    for benchmark in result.get("benchmarks", []):
        if benchmark.get("name") == name:
            value = benchmark.get(stat)
            return float(value) if value is not None else None
    return None


def print_human_report(payload: dict[str, Any]) -> None:
    commits = payload["commits"]
    if not commits:
        print("No successful benchmark results.")
        return

    result_by_sha = {result["sha"]: result for result in payload["results"]}
    ratio_by_benchmark = {comparison["benchmark"]: comparison["ratios"] for comparison in payload["comparisons"]}
    benchmark_names = ordered_benchmark_names(payload["results"])
    ratio_names = list(payload["comparisons"][0]["ratios"].keys()) if payload["comparisons"] else []

    headers = ["benchmark"]
    headers.extend(f"{commit['label']} med/min/iqr ms" for commit in commits)
    headers.extend(ratio_names)

    rows: list[list[str]] = []
    for name in benchmark_names:
        row = [name]
        for commit in commits:
            stats = benchmark_stats(result_by_sha[commit["sha"]], name)
            row.append(format_stats_ms(stats))
        ratios = ratio_by_benchmark.get(name, {})
        row.extend(format_ratio(ratios.get(ratio_name)) for ratio_name in ratio_names)
        rows.append(row)

    baseline_label = commits[0]["label"]
    print(
        "Timings are per operation. "
        f"Ratios are commit median / {baseline_label} median; >1 means the commit is slower."
    )
    print_table(headers, rows)
    if payload.get("failures"):
        print("\nFailures:")
        for failure in payload["failures"]:
            print(f"{failure['ref']} ({failure['sha'][:8]}): {failure['error']}")

    errors = benchmark_errors(payload["results"])
    if errors:
        print("\nBenchmark errors:")
        for label, name, error in errors:
            print(f"{label} {name}: {error}")


def benchmark_stats(result: dict[str, Any], name: str) -> dict[str, float] | None:
    for benchmark in result.get("benchmarks", []):
        if benchmark.get("name") == name:
            if "median" not in benchmark:
                return None
            return {
                "median": float(benchmark["median"]),
                "min": float(benchmark["min"]),
                "iqr": float(benchmark["iqr"]),
                "stdev": float(benchmark["stdev"]),
            }
    return None


def benchmark_errors(results: Sequence[dict[str, Any]]) -> list[tuple[str, str, str]]:
    errors: list[tuple[str, str, str]] = []
    for result in results:
        label = result.get("label", result.get("python_executable", "worker"))
        for benchmark in result.get("benchmarks", []):
            if "error" in benchmark:
                errors.append((label, benchmark.get("name", "unknown"), benchmark["error"]))
    return errors


def format_stats_ms(stats: dict[str, float] | None) -> str:
    if stats is None:
        return "-"
    return "/".join(format_ms(stats[key]) for key in ("median", "min", "iqr"))


def format_ms(seconds: float) -> str:
    return f"{seconds * 1000.0:.4g}"


def format_ratio(value: float | None) -> str:
    return "-" if value is None else f"{value:.3g}x"


def print_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    widths = [len(header) for header in headers]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))

    def format_row(row: Sequence[str]) -> str:
        return " | ".join(cell.ljust(widths[index]) for index, cell in enumerate(row))

    print(format_row(headers))
    print("-+-".join("-" * width for width in widths))
    for row in rows:
        print(format_row(row))


def json_ready(value: Any) -> Any:
    if isinstance(value, tuple):
        return [json_ready(item) for item in value]
    if isinstance(value, list):
        return [json_ready(item) for item in value]
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    return value


def git_stdout(repo: Path, *args: str) -> str:
    completed = run_checked(["git", "-C", str(repo), *args], capture=True)
    return completed.stdout.strip()


def run_checked(
    command: Sequence[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    capture: bool = False,
) -> subprocess.CompletedProcess[str]:
    if not capture:
        print(f"+ {shell_join(command)}", file=sys.stderr)
    completed = subprocess.run(
        [str(part) for part in command],
        cwd=str(cwd) if cwd else None,
        env=env,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
        check=False,
    )
    if completed.returncode != 0:
        details = ""
        if capture:
            details = f"\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        raise RuntimeError(f"command failed ({completed.returncode}): {shell_join(command)}{details}")
    return completed


def shell_join(command: Sequence[str | Path]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def sanitize_path_part(value: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._-")
    return sanitized or "commit"


if __name__ == "__main__":
    raise SystemExit(main())
