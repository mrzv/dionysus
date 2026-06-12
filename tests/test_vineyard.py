import dionysus as d
import pytest


PRIME = 5
BOUNDARY = [
    [],
    [],
    [(1, 0), (1, 1)],
    [(1, 0), (2, 1)],
]

TRIANGLE_BOUNDARY = [
    [],
    [],
    [],
    [(4, 0), (1, 1)],
    [(4, 1), (1, 2)],
    [(4, 0), (1, 2)],
    [(1, 3), (1, 4), (4, 5)],
]


def normalize(column):
    return [
        (coefficient % PRIME, index)
        for coefficient, index in sorted(column, key=lambda entry: entry[1])
        if coefficient % PRIME
    ]


def order(vineyard):
    return [vineyard.cell_at(i) for i in range(len(vineyard))]


def columns_from_filtration(filtration):
    columns = []
    for simplex in filtration:
        column = []
        for boundary_index, face in enumerate(simplex.boundary()):
            coefficient = 1 if boundary_index % 2 == 0 else PRIME - 1
            column.append((coefficient, filtration.index(face)))
        columns.append(sorted(column, key=lambda entry: entry[1]))
    return columns


def matrix_filtration(columns, current_order):
    position = {cell: p for p, cell in enumerate(current_order)}
    matrix = d.ReducedMatrix(d.Zp(PRIME), len(columns))

    for column_position, cell in enumerate(current_order):
        matrix[column_position] = [
            (coefficient, position[row]) for coefficient, row in columns[cell]
        ]

    return d.MatrixFiltration(matrix, [0] * len(columns), list(range(len(columns))))


def recompute_matrix_v(columns, current_order):
    reduced, chains = d.homology_persistence(
        matrix_filtration(columns, current_order), prime=PRIME, method="matrix_v"
    )

    stable_reduced = [[] for _ in current_order]
    stable_chains = [[] for _ in current_order]
    stable_pairs = [reduced.unpaired for _ in current_order]

    for column_position, cell in enumerate(current_order):
        stable_reduced[cell] = normalize(
            (entry.element, current_order[entry.index]) for entry in reduced[column_position]
        )
        stable_chains[cell] = normalize(
            (entry.element, current_order[entry.index]) for entry in chains[column_position]
        )

        pair = reduced.pair(column_position)
        if pair != reduced.unpaired:
            stable_pairs[cell] = current_order[pair]

    return stable_reduced, stable_chains, stable_pairs


def recompute_matrix_u(columns, current_order):
    reduced, trails = d.homology_persistence(
        matrix_filtration(columns, current_order), prime=PRIME, method="matrix_u"
    )

    stable_reduced = [[] for _ in current_order]
    stable_trails = [[] for _ in current_order]
    stable_pairs = [reduced.unpaired for _ in current_order]

    for position, cell in enumerate(current_order):
        stable_reduced[cell] = normalize(
            (entry.element, current_order[entry.index]) for entry in reduced[position]
        )
        stable_trails[cell] = normalize(
            (entry.element, current_order[entry.index]) for entry in trails[position]
        )

        pair = reduced.pair(position)
        if pair != reduced.unpaired:
            stable_pairs[cell] = current_order[pair]

    return stable_reduced, stable_trails, stable_pairs


def vineyard_state(vineyard):
    return (
        [vineyard.reduced_column(i) for i in range(len(vineyard))],
        [vineyard.chain(i) for i in range(len(vineyard))],
        [vineyard.pair(i) for i in range(len(vineyard))],
    )


def vineyard_u_state(vineyard):
    return (
        [vineyard.reduced_column(i) for i in range(len(vineyard))],
        [vineyard.trail(i) for i in range(len(vineyard))],
        [vineyard.pair(i) for i in range(len(vineyard))],
    )


def pairs(vineyard):
    return [vineyard.pair(i) for i in range(len(vineyard))]


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


def assert_low_pivot_caches_match_columns(vineyard):
    current_order = order(vineyard)
    position = {cell: p for p, cell in enumerate(current_order)}
    expected_lows = [vineyard.unpaired] * len(vineyard)
    expected_pivots = [vineyard.unpaired] * len(vineyard)

    for column in range(len(vineyard)):
        reduced_column = vineyard.reduced_column(column)
        if reduced_column:
            low = max(reduced_column, key=lambda entry: position[entry[1]])[1]
            expected_lows[column] = low
            expected_pivots[low] = column

    assert [vineyard.low(i) for i in range(len(vineyard))] == expected_lows
    assert [vineyard.pivot(i) for i in range(len(vineyard))] == expected_pivots


def assert_matches_recomputation(vineyard):
    expected = recompute_matrix_v(BOUNDARY, order(vineyard))

    assert vineyard_state(vineyard) == expected
    assert_low_pivot_caches_match_columns(vineyard)


def assert_u_matches_recomputation(vineyard):
    expected = recompute_matrix_u(BOUNDARY, order(vineyard))

    assert vineyard_u_state(vineyard) == expected
    assert_low_pivot_caches_match_columns(vineyard)


def assert_reconstructs_boundary(vineyard, columns):
    reduced = [vineyard.reduced_column(i) for i in range(len(vineyard))]
    trails = [vineyard.trail(i) for i in range(len(vineyard))]

    assert multiply_reduced_by_trails(reduced, trails) == [normalize(column) for column in columns]


def assert_pairs_match_recomputation(vineyard, columns):
    _, _, expected_pairs = recompute_matrix_v(columns, order(vineyard))

    assert pairs(vineyard) == expected_pairs
    assert_low_pivot_caches_match_columns(vineyard)


def transpose_and_assert_pair_locality(vineyard, position):
    before = pairs(vineyard)
    a = vineyard.cell_at(position)
    b = vineyard.cell_at(position + 1)
    local = {a, b, before[a], before[b]}
    local.discard(vineyard.unpaired)

    vineyard.transpose_position(position)

    after = pairs(vineyard)
    for cell in range(len(vineyard)):
        if cell not in local:
            assert after[cell] == before[cell]


def test_vineyard_initializes_from_matrix_v_reduction():
    vineyard = d.VineyardV(BOUNDARY, field=d.Zp(PRIME))

    assert_matches_recomputation(vineyard)


def test_vineyard_transpose_repairs_duplicate_low():
    vineyard = d.VineyardV(BOUNDARY, field=d.Zp(PRIME))

    assert vineyard.transpose_position(0) == (0, 1)

    assert order(vineyard) == [1, 0, 2, 3]
    assert_matches_recomputation(vineyard)


def test_vineyard_handles_multiple_adjacent_transpositions():
    vineyard = d.VineyardV(BOUNDARY, field=d.Zp(PRIME))

    for position in [0, 0, 2, 2]:
        vineyard.transpose_position(position)
        assert_matches_recomputation(vineyard)


def test_vineyard_handles_triangle_vertex_and_edge_transpositions():
    vineyard = d.VineyardV(TRIANGLE_BOUNDARY, field=d.Zp(PRIME))

    assert_pairs_match_recomputation(vineyard, TRIANGLE_BOUNDARY)
    for position in [0, 1, 0, 3, 4, 3]:
        transpose_and_assert_pair_locality(vineyard, position)
        assert_pairs_match_recomputation(vineyard, TRIANGLE_BOUNDARY)


def test_vineyard_u_initializes_from_matrix_u_reduction():
    vineyard = d.VineyardU(BOUNDARY, field=d.Zp(PRIME))

    assert_u_matches_recomputation(vineyard)
    assert_reconstructs_boundary(vineyard, BOUNDARY)


def test_vineyard_u_handles_multiple_adjacent_transpositions():
    vineyard = d.VineyardU(BOUNDARY, field=d.Zp(PRIME))

    for position in [0, 0, 2, 2]:
        vineyard.transpose_position(position)
        assert_u_matches_recomputation(vineyard)
        assert_reconstructs_boundary(vineyard, BOUNDARY)


def test_vineyard_u_handles_triangle_vertex_and_edge_transpositions():
    vineyard = d.VineyardU(TRIANGLE_BOUNDARY, field=d.Zp(PRIME))

    assert_pairs_match_recomputation(vineyard, TRIANGLE_BOUNDARY)
    assert_reconstructs_boundary(vineyard, TRIANGLE_BOUNDARY)
    for position in [0, 1, 0, 3, 4, 3]:
        transpose_and_assert_pair_locality(vineyard, position)
        assert_pairs_match_recomputation(vineyard, TRIANGLE_BOUNDARY)
        assert_reconstructs_boundary(vineyard, TRIANGLE_BOUNDARY)


def test_vineyard_factory_defaults_to_matrix_v():
    vineyard = d.Vineyard(BOUNDARY, field=d.Zp(PRIME))

    assert isinstance(vineyard, d.VineyardV)
    assert_matches_recomputation(vineyard)


def test_vineyard_factory_supports_matrix_u():
    vineyard = d.Vineyard(BOUNDARY, field=d.Zp(PRIME), method="matrix_u")

    assert isinstance(vineyard, d.VineyardU)
    assert_u_matches_recomputation(vineyard)
    assert_reconstructs_boundary(vineyard, BOUNDARY)


def test_vineyard_factory_rejects_unknown_method():
    with pytest.raises(ValueError, match="unknown vineyard method"):
        d.Vineyard(BOUNDARY, field=d.Zp(PRIME), method="unknown")


def test_vineyard_v_initializes_from_filtration():
    filtration = d.Filtration([[0], [1], [2], [0, 1], [1, 2], [0, 2]])
    columns = columns_from_filtration(filtration)

    from_filtration = d.VineyardV(filtration, field=d.Zp(PRIME))
    from_columns = d.VineyardV(columns, field=d.Zp(PRIME))

    assert vineyard_state(from_filtration) == vineyard_state(from_columns)
    assert order(from_filtration) == list(range(len(filtration)))


def test_vineyard_u_factory_initializes_from_filtration():
    filtration = d.Filtration([[0], [1], [2], [0, 1], [1, 2], [0, 2]])
    columns = columns_from_filtration(filtration)

    from_filtration = d.Vineyard(filtration, field=d.Zp(PRIME), method="matrix_u")
    from_columns = d.Vineyard(columns, field=d.Zp(PRIME), method="matrix_u")

    assert isinstance(from_filtration, d.VineyardU)
    assert vineyard_u_state(from_filtration) == vineyard_u_state(from_columns)
    assert order(from_filtration) == list(range(len(filtration)))


def test_vineyard_filtration_constructor_rejects_invalid_order():
    filtration = d.Filtration([[0, 1], [0], [1]])

    with pytest.raises(ValueError, match="boundary face does not precede"):
        d.Vineyard(filtration, field=d.Zp(PRIME))
