import math

import dionysus as d


def simplices(filtration):
    return [(list(simplex), simplex.data) for simplex in filtration]


def test_fast_zigzag_builds_ordered_cone():
    f = d.Filtration([[0], [1], [0, 1]])
    times = [[0], [1], [1, 3]]

    cone = d.fast_zigzag(f, times)
    result = simplices(cone)

    assert result == [
        ([-1], -math.inf),
        ([0], 0.0),
        ([1], 1.0),
        ([0, 1], 1.0),
        ([-1, 0], math.inf),
        ([-1, 1], math.inf),
        ([-1, 0, 1], 3.0),
    ]


def test_fast_zigzag_matrix_v_regression():
    f = d.Filtration([[0], [1], [0, 1], [2], [1, 2]])
    times = [[0], [1], [1, 3, 3], [2], [2]]

    cone = d.fast_zigzag(f, times)
    r, v = d.homology_persistence(cone, method='matrix_v')

    assert len(cone) > 0
    assert len(r) == len(cone)
    assert len(v) == len(cone)
