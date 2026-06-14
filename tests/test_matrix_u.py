import pickle

import dionysus as d


def chain_entries(chain):
    return [(entry.element, entry.index) for entry in chain]


def matrix_columns(matrix):
    return [chain_entries(matrix[i]) for i in range(len(matrix))]


def matrix_pairs(matrix):
    return [matrix.pair(i) for i in range(len(matrix))]


def normalize_column(column, prime):
    return [(element, index) for index, element in sorted(column.items()) if element % prime]


def multiply_reduced_by_trails(reduced, trails, prime):
    columns = [{} for _ in range(len(reduced))]
    for row, trail in enumerate(trails):
        for trail_entry in trail:
            column = columns[trail_entry.index]
            for reduced_entry in reduced[row]:
                index = reduced_entry.index
                column[index] = (column.get(index, 0) + trail_entry.element * reduced_entry.element) % prime

    return [normalize_column(column, prime) for column in columns]


def identity_matrix(size):
    return [[(1, i)] for i in range(size)]


def matrix_filtration(prime=5):
    matrix = d.ReducedMatrix(d.Zp(prime), 4)
    matrix[0] = []
    matrix[1] = []
    matrix[2] = [(1, 0), (1, 1)]
    matrix[3] = [(1, 0), (2, 1)]
    return d.MatrixFiltration(matrix, [0, 0, 1, 1], [0, 1, 2, 3])


def test_matrix_u_reconstructs_boundary_matrix():
    prime = 5
    expected_boundary = [[], [], [(1, 0), (1, 1)], [(1, 0), (2, 1)]]

    reduced, trails = d.homology_persistence(matrix_filtration(prime), prime=prime, method="matrix_u")

    assert multiply_reduced_by_trails(reduced, trails, prime) == expected_boundary


def test_matrix_u_is_inverse_of_matrix_v():
    prime = 5

    _, chains = d.homology_persistence(matrix_filtration(prime), prime=prime, method="matrix_v")
    _, trails = d.homology_persistence(matrix_filtration(prime), prime=prime, method="matrix_u")

    assert multiply_reduced_by_trails(chains, trails, prime) == identity_matrix(len(chains))


def test_matrix_u_matches_column_reduction():
    prime = 5

    reduced, _ = d.homology_persistence(matrix_filtration(prime), prime=prime, method="matrix_u")
    column_reduced = d.homology_persistence(matrix_filtration(prime), prime=prime, method="column")

    assert matrix_columns(reduced) == matrix_columns(column_reduced)
    assert matrix_pairs(reduced) == matrix_pairs(column_reduced)


def test_matrix_u_no_negative_matches_column_no_negative():
    prime = 5

    reduced, trails = d.homology_persistence(matrix_filtration(prime), prime=prime, method="matrix_u_no_negative")
    column_reduced = d.homology_persistence(matrix_filtration(prime), prime=prime, method="column_no_negative")

    assert len(trails) == len(reduced)
    assert matrix_columns(reduced) == matrix_columns(column_reduced)
    assert matrix_pairs(reduced) == matrix_pairs(column_reduced)


def test_matrix_u_skipped_columns_match_matrix_v():
    filtration = d.Filtration([[0], [1], [0, 1]])
    relative = d.Filtration([[0]])

    _, chains = d.homology_persistence(filtration, relative, method="matrix_v")
    _, trails = d.homology_persistence(filtration, relative, method="matrix_u")

    assert chain_entries(chains[0]) == []
    assert chain_entries(trails[0]) == []
    assert chain_entries(chains[1]) == [(1, 1)]
    assert chain_entries(trails[1]) == [(1, 1)]


def test_reduced_matrix_pickle_roundtrip():
    prime = 5
    reduced, _ = d.homology_persistence(matrix_filtration(prime), prime=prime, method="matrix_u")

    restored = pickle.loads(pickle.dumps(reduced))

    assert restored.field().prime() == prime
    assert matrix_columns(restored) == matrix_columns(reduced)
    assert matrix_pairs(restored) == matrix_pairs(reduced)
