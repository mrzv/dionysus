import math

import dionysus as d
import pytest


PRIME = 5


def simplex_key(simplex):
    return (simplex.dimension(), tuple(simplex))


def stable_maps(filtration, values0):
    stable_to_input = sorted(
        range(len(filtration)),
        key=lambda index: (values0[index], simplex_key(filtration[index])),
    )
    input_to_stable = {
        input_index: stable for stable, input_index in enumerate(stable_to_input)
    }
    return stable_to_input, input_to_stable


def endpoint_order(filtration, values0, values1):
    stable_to_input, _ = stable_maps(filtration, values0)
    return sorted(
        range(len(filtration)),
        key=lambda stable: (
            values1[stable_to_input[stable]],
            simplex_key(filtration[stable_to_input[stable]]),
        ),
    )


def endpoint_matrix_filtration(filtration, values0, values1):
    stable_to_input, input_to_stable = stable_maps(filtration, values0)
    current_order = endpoint_order(filtration, values0, values1)
    position = {cell: p for p, cell in enumerate(current_order)}
    matrix = d.ReducedMatrix(d.Zp(PRIME), len(filtration))

    for column_position, stable in enumerate(current_order):
        simplex = filtration[stable_to_input[stable]]
        column = []
        for boundary_index, face in enumerate(simplex.boundary()):
            coefficient = 1 if boundary_index % 2 == 0 else PRIME - 1
            face_stable = input_to_stable[filtration.index(face)]
            column.append((coefficient, position[face_stable]))
        matrix[column_position] = column

    dimensions = [filtration[stable_to_input[cell]].dimension() for cell in current_order]
    values = [values1[stable_to_input[cell]] for cell in current_order]
    return d.MatrixFiltration(matrix, dimensions, values), current_order


def recompute_endpoint_pairs(filtration, values0, values1, method):
    matrix_filtration, current_order = endpoint_matrix_filtration(filtration, values0, values1)
    result = d.homology_persistence(matrix_filtration, prime=PRIME, method=method)
    reduced = result[0] if isinstance(result, tuple) else result

    pairs = [reduced.unpaired] * len(current_order)
    for position, stable in enumerate(current_order):
        pair = reduced.pair(position)
        if pair != reduced.unpaired:
            pairs[stable] = current_order[pair]
    return pairs, current_order


def assert_final_pairs_match_recomputation(filtration, values0, values1, method="matrix_v"):
    result = d.vineyard_linear_homotopy(
        filtration, values0, values1, field=d.Zp(PRIME), method=method
    )
    expected_pairs, expected_order = recompute_endpoint_pairs(
        filtration, values0, values1, method
    )

    assert result.final_order == expected_order
    assert [result.vineyard.pair(i) for i in range(len(filtration))] == expected_pairs


def test_linear_homotopy_records_single_crossing_and_vines():
    filtration = d.Filtration([[0], [1]])

    result = d.vineyard_linear_homotopy(
        filtration, [0.0, 1.0], [1.0, 0.0], field=d.Zp(PRIME)
    )

    assert isinstance(result.vineyard, d.VineyardV)
    assert len(result.events) == 1
    assert result.events[0].time == pytest.approx(0.5)
    assert (result.events[0].first, result.events[0].second) == (0, 1)
    assert result.events[0].first_pair_before == result.vineyard.unpaired
    assert result.events[0].second_pair_before == result.vineyard.unpaired
    assert result.events[0].first_pair_after == result.vineyard.unpaired
    assert result.events[0].second_pair_after == result.vineyard.unpaired
    assert result.final_order == [1, 0]

    segments = [segment for vine in result.vines for segment in vine.segments]
    assert len(segments) == 4
    assert {segment.event1 for segment in segments if segment.t1 == pytest.approx(0.5)} == {0}
    assert {segment.event0 for segment in segments if segment.t0 == pytest.approx(0.5)} == {0}
    assert all(math.isinf(segment.death0) and math.isinf(segment.death1) for segment in segments)


def test_linear_homotopy_processes_simultaneous_degenerate_block():
    filtration = d.Filtration([[0], [1], [2]])

    result = d.vineyard_linear_homotopy(
        filtration, [0.0, 1.0, 2.0], [2.0, 1.0, 0.0], field=d.Zp(PRIME)
    )

    assert [event.time for event in result.events] == [pytest.approx(0.5)] * 3
    assert [(event.first, event.second) for event in result.events] == [
        (1, 2),
        (0, 2),
        (0, 1),
    ]
    assert result.final_order == [2, 1, 0]


def test_linear_homotopy_processes_disjoint_simultaneous_crossings():
    filtration = d.Filtration([[0], [1], [2], [3]])

    result = d.vineyard_linear_homotopy(
        filtration, [0.0, 1.0, 10.0, 11.0], [1.0, 0.0, 11.0, 10.0], field=d.Zp(PRIME)
    )

    assert [event.time for event in result.events] == [pytest.approx(0.5)] * 2
    assert {(event.first, event.second) for event in result.events} == {(0, 1), (2, 3)}
    assert result.final_order == [1, 0, 3, 2]


def test_linear_homotopy_processes_initial_equal_value_inversions():
    filtration = d.Filtration([[0], [1], [2]])

    result = d.vineyard_linear_homotopy(
        filtration, [0.0, 0.0, 0.0], [2.0, 1.0, 0.0], field=d.Zp(PRIME)
    )

    assert [event.time for event in result.events] == [pytest.approx(0.0)] * 3
    assert result.final_order == [2, 1, 0]


def test_linear_homotopy_preserves_filtration_order_with_dimension_ties():
    filtration = d.Filtration([[0], [1], [0, 1]])

    result = d.vineyard_linear_homotopy(
        filtration, [0.0, 1.0, 1.0], [1.0, 0.0, 1.0], field=d.Zp(PRIME)
    )

    assert len(result.events) == 1
    assert result.final_order == [1, 0, 2]
    assert result.vineyard.position(0) < result.vineyard.position(2)
    assert result.vineyard.position(1) < result.vineyard.position(2)


def test_linear_homotopy_supports_matrix_u():
    filtration = d.Filtration([[0], [1], [0, 1]])

    result = d.vineyard_linear_homotopy(
        filtration,
        [0.0, 1.0, 1.0],
        [1.0, 0.0, 1.0],
        field=d.Zp(PRIME),
        method="matrix_u",
    )

    assert isinstance(result.vineyard, d.VineyardU)
    assert len(result.events) == 1
    assert result.final_order == [1, 0, 2]
    assert result.vineyard.position(0) < result.vineyard.position(2)
    assert result.vineyard.position(1) < result.vineyard.position(2)


def test_linear_homotopy_final_pairs_match_matrix_v_recomputation():
    filtration = d.Filtration([[0], [1], [2], [0, 1], [1, 2], [0, 2]])
    values0 = [0.0, 1.0, 2.0, 2.0, 3.0, 3.0]
    values1 = [2.0, 0.0, 1.0, 2.0, 1.0, 3.0]

    assert_final_pairs_match_recomputation(filtration, values0, values1)


def test_linear_homotopy_final_pairs_match_matrix_u_recomputation():
    filtration = d.Filtration([[0], [1], [2], [0, 1], [1, 2], [0, 2]])
    values0 = [0.0, 1.0, 2.0, 2.0, 3.0, 3.0]
    values1 = [2.0, 0.0, 1.0, 2.0, 1.0, 3.0]

    assert_final_pairs_match_recomputation(
        filtration, values0, values1, method="matrix_u"
    )


def test_linear_homotopy_rejects_non_filtration_endpoint_values():
    filtration = d.Filtration([[0], [1], [0, 1]])

    with pytest.raises(ValueError, match="not a filtration"):
        d.vineyard_linear_homotopy(
            filtration, [0.0, 0.0, 1.0], [0.0, 1.0, 0.5], field=d.Zp(PRIME)
        )


def test_linear_homotopy_rejects_unknown_method():
    filtration = d.Filtration([[0]])

    with pytest.raises(ValueError, match="unknown vineyard method"):
        d.vineyard_linear_homotopy(filtration, [0.0], [1.0], method="unknown")
