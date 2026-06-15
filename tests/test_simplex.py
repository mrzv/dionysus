import pytest

import dionysus as d


def simplex_summary(simplex):
    return list(simplex), simplex.data, repr(simplex)


def test_simplex_construction_sorts_vertices_and_exposes_sequence_api():
    simplex = d.Simplex([2, 0, 1], 7.0)

    assert list(simplex) == [0, 1, 2]
    assert len(simplex) == 3
    assert simplex.dimension() == 2
    assert simplex[0] == 0
    assert simplex[2] == 2
    assert 1 in simplex
    assert 3 not in simplex
    assert repr(simplex) == "<0,1,2> 7"


def test_empty_simplex_has_minus_one_dimension():
    simplex = d.Simplex()

    assert list(simplex) == []
    assert len(simplex) == 0
    assert simplex.dimension() == -1
    assert repr(simplex) == "<> 0"


def test_simplex_boundary_order_and_default_face_data():
    simplex = d.Simplex([0, 1, 2], 7.0)

    faces = list(simplex.boundary())

    assert [simplex_summary(face) for face in faces] == [
        ([1, 2], 0.0, "<1,2> 0"),
        ([0, 2], 0.0, "<0,2> 0"),
        ([0, 1], 0.0, "<0,1> 0"),
    ]


def test_simplex_join_returns_new_sorted_simplex_without_mutating_original():
    simplex = d.Simplex([1], 2.0)

    joined = simplex.join(0)

    assert simplex_summary(simplex) == ([1], 2.0, "<1> 2")
    assert simplex_summary(joined) == ([0, 1], 2.0, "<0,1> 2")


def test_simplex_data_is_mutable_but_not_part_of_identity():
    first = d.Simplex([0], 0.0)
    second = d.Simplex([0], 1.0)

    assert first == second
    assert hash(first) == hash(second)
    assert not first < second

    first.data = 3.0
    assert first.data == 3.0
    assert first == second


def test_simplex_index_errors_propagate():
    with pytest.raises(IndexError):
        d.Simplex([0])[1]
