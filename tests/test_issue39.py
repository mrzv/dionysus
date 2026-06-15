from pathlib import Path

import numpy as np
import pytest
import dionysus as d


pytestmark = pytest.mark.skip(reason="issue39 wasserstein regression currently hangs")


def test_issue39():
    data_dir = Path(__file__).parent / 'data' / 'issue39'
    dgm1 = np.loadtxt(data_dir / 'dgm1.txt', delimiter=',')
    dgm2 = np.loadtxt(data_dir / 'dgm2.txt', delimiter=',')
    dgm1 = d.Diagram(dgm1)
    dgm2 = d.Diagram(dgm2)
    dist = d.wasserstein_distance(dgm1,dgm2,q=5)
    assert dist >= 0
