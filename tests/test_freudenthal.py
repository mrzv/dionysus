import numpy as np
import pytest

import dionysus as d


def simplices(filtration):
    return [(list(simplex), simplex.data) for simplex in filtration]


def test_fill_freudenthal_1d_lower_star():
    f = d.fill_freudenthal(np.array([2.0, 1.0, 3.0], dtype=np.float32))

    assert simplices(f) == [
        ([1], 1.0),
        ([0], 2.0),
        ([0, 1], 2.0),
        ([2], 3.0),
        ([1, 2], 3.0),
    ]


def test_fill_freudenthal_1d_upper_star():
    f = d.fill_freudenthal(np.array([2.0, 1.0, 3.0], dtype=np.float64), reverse=True)

    assert simplices(f) == [
        ([2], 3.0),
        ([0], 2.0),
        ([1], 1.0),
        ([0, 1], 1.0),
        ([1, 2], 1.0),
    ]


def test_fill_freudenthal_2d_lower_star():
    f = d.fill_freudenthal(np.array([[0.0, 2.0], [1.0, 3.0]], dtype=np.float32))

    assert simplices(f) == [
        ([0], 0.0),
        ([2], 1.0),
        ([0, 2], 1.0),
        ([1], 2.0),
        ([0, 1], 2.0),
        ([3], 3.0),
        ([0, 3], 3.0),
        ([1, 3], 3.0),
        ([2, 3], 3.0),
        ([0, 1, 3], 3.0),
        ([0, 2, 3], 3.0),
    ]


def test_fill_freudenthal_rejects_unknown_dtype():
    try:
        d.fill_freudenthal(np.array([1, 2, 3], dtype=np.int64))
    except RuntimeError as e:
        assert "Unknown array dtype" in str(e)
    else:
        raise AssertionError("expected RuntimeError")


def test_fill_freudenthal_rejects_negative_strides():
    data = np.array([1.0, 2.0, 3.0], dtype=np.float64)[::-1]

    with pytest.raises(RuntimeError, match="negative array strides"):
        d.fill_freudenthal(data)
