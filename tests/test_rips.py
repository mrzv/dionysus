import numpy as np
import pytest

import dionysus as d


def simplex_summary(filtration):
    return [(list(simplex), simplex.data) for simplex in filtration]


def test_fill_rips_pairwise_points_store_unsquared_distances():
    points = np.array([[0.0], [3.0], [4.0]], dtype=np.float32)

    filtration = d.fill_rips(points, 1, 4.0)

    assert simplex_summary(filtration) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([1, 2], 1.0),
        ([0, 1], 3.0),
        ([0, 2], 4.0),
    ]


def test_fill_rips_explicit_distances_use_input_values():
    # Lower-triangular distances for pairs (1, 0), (2, 0), (2, 1).
    distances = np.array([3.0, 4.0, 1.0], dtype=np.float64)

    filtration = d.fill_rips(distances, 1, 4.0)

    assert simplex_summary(filtration) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([1, 2], 1.0),
        ([0, 1], 3.0),
        ([0, 2], 4.0),
    ]


def test_fill_rips_rejects_unknown_dtype():
    with pytest.raises(RuntimeError, match="Unknown array dtype"):
        d.fill_rips(np.array([[0], [1]], dtype=np.int64), 1, 1.0)
