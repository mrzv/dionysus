from pathlib import Path
import subprocess
import sys

import pytest


def test_issue39():
    data_dir = Path(__file__).parent / "data" / "issue39"
    script = """
import sys
from pathlib import Path

import dionysus as d
import numpy as np

data_dir = Path(sys.argv[1])
dgm1 = d.Diagram(np.loadtxt(data_dir / "dgm1.txt", delimiter=","))
dgm2 = d.Diagram(np.loadtxt(data_dir / "dgm2.txt", delimiter=","))
print(d.wasserstein_distance(dgm1, dgm2, q=5))
"""
    result = subprocess.run(
        [sys.executable, "-c", script, str(data_dir)],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )

    # Hera guarantees delta-relative accuracy; use the exact 1-D reference value.
    assert float(result.stdout) == pytest.approx(0.009632668295414, rel=0.01)
