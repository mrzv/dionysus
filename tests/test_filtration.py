import pytest

import dionysus as d


def simplex_summary(filtration):
    return [(list(simplex), simplex.data) for simplex in filtration]


def test_filtration_construction_sort_index_contains_and_rearrange():
    filtration = d.Filtration([
        ([0, 1], 2.0),
        ([0], 0.0),
        ([1], 1.0),
    ])

    assert len(filtration) == 3
    assert d.Simplex([0]) in filtration
    assert d.Simplex([2]) not in filtration

    filtration.sort()
    assert simplex_summary(filtration) == [
        ([0], 0.0),
        ([1], 1.0),
        ([0, 1], 2.0),
    ]
    assert filtration.index(d.Simplex([0])) == 0
    assert filtration.index(d.Simplex([1])) == 1
    assert filtration.index(d.Simplex([0, 1])) == 2

    filtration.rearrange([2, 0, 1])
    assert [list(simplex) for simplex in filtration] == [[0, 1], [0], [1]]
    assert filtration.index(d.Simplex([0, 1])) == 0
    assert filtration.index(d.Simplex([0])) == 1
    assert filtration.index(d.Simplex([1])) == 2


def test_filtration_keeps_unique_simplices():
    filtration = d.Filtration([d.Simplex([0], 0.0), d.Simplex([0], 1.0)])

    assert len(filtration) == 1
    assert simplex_summary(filtration) == [([0], 0.0)]


def test_multi_filtration_preserves_duplicate_simplices():
    filtration = d.MultiFiltration([d.Simplex([0], 0.0), d.Simplex([0], 1.0)])

    assert len(filtration) == 2
    assert simplex_summary(filtration) == [([0], 0.0), ([0], 1.0)]

    filtration.sort(reverse=True)
    assert simplex_summary(filtration) == [([0], 1.0), ([0], 0.0)]
    assert filtration.index(d.Simplex([0], 1.0), 0) == 0
    assert filtration.index(d.Simplex([0], 0.0), 1) == 1


def test_multi_filtration_index_respects_search_bound():
    filtration = d.MultiFiltration([d.Simplex([0], 0.0), d.Simplex([0], 1.0), d.Simplex([0], 2.0)])

    assert filtration.index(d.Simplex([0], 99.0), 2) == 2
    assert filtration.index(d.Simplex([0], 99.0), 1) == 1
    assert filtration.index(d.Simplex([0], 99.0), 0) == 0

    with pytest.raises(RuntimeError, match="Trying to access non-existent cell"):
        filtration.index(d.Simplex([1]), 3)


def test_multi_filtration_checked_lookup_rejects_missing_predecessor():
    filtration = d.MultiFiltration([d.Simplex([1], 0.0)])

    with pytest.raises(RuntimeError, match="Trying to access non-existent cell"):
        filtration.index(d.Simplex([0]), 0)

    with pytest.raises(RuntimeError, match="Trying to access non-existent cell"):
        d.MultiFiltration().index(d.Simplex([0]))


def test_multi_filtration_boundary_lookup_uses_latest_duplicate_before_cell():
    filtration = d.MultiFiltration([
        d.Simplex([0], 0.0),
        d.Simplex([0], 1.0),
        d.Simplex([1], 0.0),
        d.Simplex([0, 1], 2.0),
    ])

    edge_index = 3
    assert [filtration.index(face, edge_index) for face in filtration[edge_index].boundary()] == [2, 1]


def test_linked_multi_filtration_sort_rearrange_and_index():
    filtration = d.LinkedMultiFiltration([
        ([0], 2.0, 0),
        ([1], 1.0, 1),
        ([0, 1], 3.0, 2),
    ])

    assert simplex_summary(filtration) == [
        ([0], 2.0),
        ([1], 1.0),
        ([0, 1], 3.0),
    ]

    filtration.sort()
    assert simplex_summary(filtration) == [
        ([1], 1.0),
        ([0], 2.0),
        ([0, 1], 3.0),
    ]
    assert filtration.index(d.Simplex([1], 1.0)) == 0
    assert filtration.index(d.Simplex([0], 2.0)) == 1
    assert filtration.index(d.Simplex([0, 1], 3.0)) == 2

    filtration.rearrange([1, 0, 2])
    assert simplex_summary(filtration) == [
        ([0], 2.0),
        ([1], 1.0),
        ([0, 1], 3.0),
    ]
    assert filtration.index(d.Simplex([0], 2.0)) == 0
    assert filtration.index(d.Simplex([1], 1.0)) == 1
    assert filtration.index(d.Simplex([0, 1], 3.0)) == 2


def test_linked_multi_filtration_preserves_duplicate_simplices():
    filtration = d.LinkedMultiFiltration([
        (d.Simplex([0], 0.0), 0),
        (d.Simplex([0], 1.0), 1),
    ])

    assert len(filtration) == 2
    assert simplex_summary(filtration) == [([0], 0.0), ([0], 1.0)]

    filtration.sort(reverse=True)
    assert simplex_summary(filtration) == [([0], 1.0), ([0], 0.0)]


def test_linked_multi_filtration_default_index_uses_duplicate_lookup():
    filtration = d.LinkedMultiFiltration([
        (d.Simplex([0], 0.0), 0),
        (d.Simplex([0], 1.0), 1),
    ])

    assert filtration.index(d.Simplex([0], 99.0)) == 1


def test_linked_multi_filtration_falls_back_when_link_does_not_match():
    filtration = d.LinkedMultiFiltration([
        (d.Simplex([0], 0.0), 0),
        (d.Simplex([1], 1.0), 1),
    ])

    assert filtration.index(d.Simplex([0], 99.0), 1) == 0


def test_linked_multi_filtration_checked_lookup_rejects_missing_predecessor():
    filtration = d.LinkedMultiFiltration([(d.Simplex([1], 0.0), 0)])

    with pytest.raises(RuntimeError, match="Trying to access non-existent cell"):
        filtration.index(d.Simplex([0]), 0)

    with pytest.raises(RuntimeError, match="Trying to access non-existent cell"):
        d.LinkedMultiFiltration().index(d.Simplex([0]))


def test_linked_multi_filtration_boundary_lookup_uses_latest_duplicate_before_cell():
    filtration = d.LinkedMultiFiltration([
        (d.Simplex([0], 0.0), 0),
        (d.Simplex([0], 1.0), 1),
        (d.Simplex([1], 0.0), 2),
        (d.Simplex([0, 1], 2.0), 3),
    ])

    edge_index = 3
    assert [filtration.index(face, edge_index) for face in filtration[edge_index].boundary()] == [2, 1]
