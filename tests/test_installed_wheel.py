import os
import subprocess
import sys
from pathlib import Path

import pytest


pytestmark = pytest.mark.packaging


def run_command(command, **kwargs):
    subprocess.run(command, check=True, **kwargs)


def test_installed_wheel_imports_from_virtualenv(tmp_path):
    if os.environ.get("DIONYSUS_RUN_PACKAGING_TESTS") != "1":
        pytest.skip("set DIONYSUS_RUN_PACKAGING_TESTS=1 to build and install the wheel")

    root = Path(__file__).resolve().parents[1]
    wheelhouse = tmp_path / "wheelhouse"
    venv = tmp_path / "venv"

    run_command(["uv", "build", "--wheel", "--out-dir", str(wheelhouse)], cwd=root)
    wheels = sorted(wheelhouse.glob("dionysus-*.whl"))
    assert len(wheels) == 1

    run_command([sys.executable, "-m", "venv", str(venv)])
    python = venv / "bin" / "python"
    run_command([str(python), "-m", "pip", "install", str(wheels[0])], cwd=tmp_path)
    run_command(
        [
            str(python),
            "-c",
            "import dionysus as d; f = d.Filtration([[0], [1], [0, 1]]); assert len(f) == 3",
        ],
        cwd=tmp_path,
        env={key: value for key, value in os.environ.items() if key != "PYTHONPATH"},
    )
