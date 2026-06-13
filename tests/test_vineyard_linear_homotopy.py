import math
import random
from itertools import combinations

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


def normalize(column):
    return [
        (coefficient % PRIME, index)
        for coefficient, index in sorted(column, key=lambda entry: entry[1])
        if coefficient % PRIME
    ]


def recompute_endpoint_matrix_v(filtration, values0, values1):
    matrix_filtration, current_order = endpoint_matrix_filtration(filtration, values0, values1)
    reduced, chains = d.homology_persistence(
        matrix_filtration, prime=PRIME, method="matrix_v"
    )

    stable_reduced = [[] for _ in current_order]
    stable_chains = [[] for _ in current_order]
    stable_pairs = [reduced.unpaired for _ in current_order]

    for position, stable in enumerate(current_order):
        stable_reduced[stable] = normalize(
            (entry.element, current_order[entry.index]) for entry in reduced[position]
        )
        stable_chains[stable] = normalize(
            (entry.element, current_order[entry.index]) for entry in chains[position]
        )

        pair = reduced.pair(position)
        if pair != reduced.unpaired:
            stable_pairs[stable] = current_order[pair]

    return stable_reduced, stable_chains, stable_pairs, current_order


def stable_boundary_columns(filtration, values0):
    stable_to_input, input_to_stable = stable_maps(filtration, values0)
    columns = [[] for _ in stable_to_input]
    for stable, input_index in enumerate(stable_to_input):
        simplex = filtration[input_index]
        column = []
        for boundary_index, face in enumerate(simplex.boundary()):
            coefficient = 1 if boundary_index % 2 == 0 else PRIME - 1
            column.append((coefficient, input_to_stable[filtration.index(face)]))
        columns[stable] = normalize(column)
    return columns


def multiply_reduced_by_trails(reduced, trails):
    columns = [{} for _ in range(len(reduced))]
    for row, trail in enumerate(trails):
        for coefficient, column_index in trail:
            column = columns[column_index]
            for reduced_coefficient, reduced_index in reduced[row]:
                column[reduced_index] = (
                    column.get(reduced_index, 0) + coefficient * reduced_coefficient
                ) % PRIME

    return [normalize((coefficient, index) for index, coefficient in column.items()) for column in columns]


def assert_final_pairs_match_recomputation(filtration, values0, values1, method="matrix_v"):
    result = d.vineyard_linear_homotopy(
        filtration, values0, values1, field=d.Zp(PRIME), method=method
    )
    expected_pairs, expected_order = recompute_endpoint_pairs(
        filtration, values0, values1, method
    )

    assert result.final_order == expected_order
    assert [result.vineyard.pair(i) for i in range(len(filtration))] == expected_pairs


def complete_simplex_skeleton(points, skeleton):
    return [
        simplex
        for dimension in range(skeleton + 1)
        for simplex in combinations(range(points), dimension + 1)
    ]


def shuffled_dimension_filtration(simplices, rng):
    order = []
    max_dimension = max(len(simplex) for simplex in simplices) - 1
    for dimension in range(max_dimension + 1):
        cells = [
            i for i, simplex in enumerate(simplices) if len(simplex) == dimension + 1
        ]
        rng.shuffle(cells)
        order.extend(cells)

    return [simplices[i] for i in order]


def random_adjacent_filtration_order(size, dimensions, rng, swaps):
    order = list(range(size))
    positions_by_dimension = {
        dimension: [i for i, d in enumerate(dimensions) if d == dimension]
        for dimension in set(dimensions)
    }

    for _ in range(swaps):
        dimension = rng.choice(list(positions_by_dimension))
        positions = positions_by_dimension[dimension]
        position = rng.choice(positions[:-1])
        order[position], order[position + 1] = order[position + 1], order[position]

    return order


def values_for_order(order):
    values = [0.0] * len(order)
    for value, cell in enumerate(order):
        values[cell] = float(value)
    return values


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
    assert len(result.vines) == 2
    assert [len(vine.segments) for vine in result.vines] == [2, 2]
    assert len(segments) == 4
    assert {segment.event1 for segment in segments if segment.t1 == pytest.approx(0.5)} == {0}
    assert {segment.event0 for segment in segments if segment.t0 == pytest.approx(0.5)} == {0}
    assert all(math.isinf(segment.death0) and math.isinf(segment.death1) for segment in segments)


def test_linear_homotopy_continues_vines_through_pairing_switch():
    filtration = d.Filtration([[0], [1], [0, 1]])

    result = d.vineyard_linear_homotopy(
        filtration, [0.0, 1.0, 2.0], [1.0, 0.0, 2.0], field=d.Zp(PRIME)
    )

    assert len(result.events) == 1
    assert len(result.vines) == 2
    assert [len(vine.segments) for vine in result.vines] == [2, 2]

    essential, finite = result.vines
    assert essential.segments[0].birth_cell == 0
    assert essential.segments[0].death_cell == result.vineyard.unpaired
    assert essential.segments[1].birth_cell == 1
    assert essential.segments[1].death_cell == result.vineyard.unpaired

    assert finite.segments[0].birth_cell == 1
    assert finite.segments[0].death_cell == 2
    assert finite.segments[1].birth_cell == 0
    assert finite.segments[1].death_cell == 2


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


def test_linear_homotopy_complete_simplex_2_skeleton_on_50_points():
    rng = random.Random(20260612)
    simplices = shuffled_dimension_filtration(
        complete_simplex_skeleton(points=50, skeleton=2), rng
    )
    filtration = d.Filtration(simplices)
    values0 = [float(i) for i in range(len(filtration))]
    dimensions = [simplex.dimension() for simplex in filtration]
    final_order = random_adjacent_filtration_order(
        len(filtration), dimensions, rng, swaps=500
    )
    values1 = values_for_order(final_order)

    vineyard = d.Vineyard(filtration, field=d.Zp(PRIME))
    expected_start_pairs, expected_start_order = recompute_endpoint_pairs(
        filtration, values0, values0, method="matrix_v"
    )

    assert expected_start_order == list(range(len(filtration)))
    assert [vineyard.pair(i) for i in range(len(filtration))] == expected_start_pairs

    result = d.vineyard_linear_homotopy(filtration, values0, values1, field=d.Zp(PRIME))
    expected_reduced, expected_chains, expected_final_pairs, expected_final_order = (
        recompute_endpoint_matrix_v(filtration, values0, values1)
    )
    vineyard_reduced = [result.vineyard.reduced_column(i) for i in range(len(filtration))]
    vineyard_chains = [result.vineyard.chain(i) for i in range(len(filtration))]

    assert result.final_order == expected_final_order
    assert [result.vineyard.pair(i) for i in range(len(filtration))] == expected_final_pairs
    assert vineyard_reduced != expected_reduced
    assert vineyard_chains != expected_chains

    lazy_result = d.vineyard_linear_homotopy(
        filtration, values0, values1, field=d.Zp(PRIME), lazy=True
    )
    lazy_reduced = [
        lazy_result.vineyard.reduced_column(i) for i in range(len(filtration))
    ]
    lazy_chains = [lazy_result.vineyard.chain(i) for i in range(len(filtration))]

    assert lazy_result.final_order == expected_final_order
    assert [lazy_result.vineyard.pair(i) for i in range(len(filtration))] == expected_final_pairs
    assert lazy_reduced == expected_reduced
    assert lazy_chains == expected_chains


def test_linear_homotopy_lazy_matrix_u_preserves_ru_decomposition():
    filtration = d.Filtration([[0], [1], [2], [0, 1], [1, 2], [0, 2]])
    values0 = [0.0, 1.0, 2.0, 2.0, 3.0, 3.0]
    values1 = [2.0, 0.0, 1.0, 2.0, 1.0, 3.0]

    result = d.vineyard_linear_homotopy(
        filtration,
        values0,
        values1,
        field=d.Zp(PRIME),
        method="matrix_u",
        lazy=True,
    )

    reduced = [result.vineyard.reduced_column(i) for i in range(len(filtration))]
    trails = [result.vineyard.trail(i) for i in range(len(filtration))]

    assert multiply_reduced_by_trails(reduced, trails) == stable_boundary_columns(
        filtration, values0
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
