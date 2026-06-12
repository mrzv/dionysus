import dionysus as d
import pytest


PRIME = 5
BOUNDARY = [
    [],
    [],
    [(1, 0), (1, 1)],
    [(1, 0), (2, 1)],
]


def normalize(column):
    return [
        (coefficient % PRIME, index)
        for coefficient, index in sorted(column, key=lambda entry: entry[1])
        if coefficient % PRIME
    ]


def order(vineyard):
    return [vineyard.cell_at(i) for i in range(len(vineyard))]


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


def vineyard_state(vineyard):
    return (
        [vineyard.reduced_column(i) for i in range(len(vineyard))],
        [vineyard.chain(i) for i in range(len(vineyard))],
        [vineyard.pair(i) for i in range(len(vineyard))],
    )


def assert_matches_recomputation(vineyard):
    expected = recompute_matrix_v(BOUNDARY, order(vineyard))

    assert vineyard_state(vineyard) == expected


def test_vineyard_initializes_from_matrix_v_reduction():
    vineyard = d.Vineyard(d.Zp(PRIME), BOUNDARY)

    assert_matches_recomputation(vineyard)


def test_vineyard_transpose_repairs_duplicate_low():
    vineyard = d.Vineyard(d.Zp(PRIME), BOUNDARY)

    assert vineyard.transpose_position(0) == (0, 1)

    assert order(vineyard) == [1, 0, 2, 3]
    assert_matches_recomputation(vineyard)


def test_vineyard_handles_multiple_adjacent_transpositions():
    vineyard = d.Vineyard(d.Zp(PRIME), BOUNDARY)

    for position in [0, 0, 2, 2]:
        vineyard.transpose_position(position)
        assert_matches_recomputation(vineyard)


def test_vineyard_rejects_invalid_filtration_transposition():
    vineyard = d.Vineyard(d.Zp(PRIME), BOUNDARY)

    with pytest.raises(ValueError):
        vineyard.transpose_position(1)
