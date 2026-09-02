import math

import dionysus as d


def test_bottleneck_distance_returns_scalar_by_default():
    dgm1 = d.Diagram([(1, 2), (3, 4), (1.0, float("inf"))])
    dgm2 = d.Diagram([(0, 2), (3, 5), (2.0, float("inf"))])

    distance = d.bottleneck_distance(dgm1, dgm2)

    assert isinstance(distance, float)
    assert distance == 1


def test_bottleneck_distance_returns_longest_point_to_point_edge():
    dgm1 = d.Diagram([(0, 10)])
    dgm2 = d.Diagram([(2, 12)])

    distance, edge = d.bottleneck_distance(
        dgm1, dgm2, delta=0, compute_longest_edge=True
    )

    assert distance == 2
    assert edge == (0, 0)


def test_bottleneck_distance_normalizes_diagonal_edge_indices():
    point = d.Diagram([(0, 4)])
    empty = d.Diagram([])

    assert d.bottleneck_distance(
        point, empty, delta=0, compute_longest_edge=True
    ) == (2, (0, -1))
    assert d.bottleneck_distance(
        empty, point, delta=0, compute_longest_edge=True
    ) == (2, (-1, 0))


def test_bottleneck_distance_uses_sentinel_edge_for_zero_and_infinity():
    point = d.Diagram([(0, 4)])
    essential_point = d.Diagram([(1, float("inf"))])
    empty = d.Diagram([])

    assert d.bottleneck_distance(
        point, point, delta=0, compute_longest_edge=True
    ) == (0, (-1, -1))

    distance, edge = d.bottleneck_distance(
        essential_point, empty, delta=0, compute_longest_edge=True
    )
    assert math.isinf(distance)
    assert edge == (-1, -1)


def test_bottleneck_distance_returns_longest_edge_in_approximate_mode():
    dgm1 = d.Diagram([(0, 10)])
    dgm2 = d.Diagram([(2, 12)])

    distance, edge = d.bottleneck_distance(
        dgm1, dgm2, delta=0.01, compute_longest_edge=True
    )

    assert 2 <= distance <= 2.02
    assert edge == (0, 0)


def test_bottleneck_distance_selects_edge_from_dominating_component():
    for delta, essential_birth in ((0, 1.999), (0.01, 1.99)):
        dgm1 = d.Diagram(
            [(0, float("inf")), (0, 4)]
        )
        dgm2 = d.Diagram([(essential_birth, float("inf"))])

        distance, edge = d.bottleneck_distance(
            dgm1, dgm2, delta=delta, compute_longest_edge=True
        )

        assert distance == 2
        assert edge == (1, -1)
