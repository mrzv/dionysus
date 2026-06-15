import dionysus as d


def chain_entries(chain):
    return [(entry.element, repr(entry.index)) for entry in chain]


def matrix_columns(matrix_filtration):
    return [chain_entries(matrix_filtration[i].boundary()) for i in range(len(matrix_filtration))]


def assert_boundary_matrix(matrix_filtration):
    assert matrix_filtration.dimensions() == [0, 0, 0, 1, 1, 1, 2]
    assert matrix_filtration.values() == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    assert matrix_columns(matrix_filtration) == [
        [],
        [],
        [],
        [(-1, "Cell 0"), (1, "Cell 1")],
        [(-1, "Cell 0"), (1, "Cell 2")],
        [(-1, "Cell 1"), (1, "Cell 2")],
        [(1, "Cell 3"), (-1, "Cell 4"), (1, "Cell 5")],
    ]


def assert_coboundary_matrix(matrix_filtration):
    assert matrix_filtration.dimensions() == [2, 1, 1, 1, 0, 0, 0]
    assert matrix_filtration.values() == [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0]
    assert matrix_columns(matrix_filtration) == [
        [],
        [(1, "Cell 0")],
        [(-1, "Cell 0")],
        [(1, "Cell 0")],
        [(1, "Cell 1"), (1, "Cell 2")],
        [(-1, "Cell 1"), (1, "Cell 3")],
        [(-1, "Cell 2"), (-1, "Cell 3")],
    ]


def test_boundary_matrix_from_filtration(triangle_cells):
    assert_boundary_matrix(d.boundary(d.Filtration(triangle_cells)))


def test_coboundary_matrix_from_filtration(triangle_cells):
    assert_coboundary_matrix(d.coboundary(d.Filtration(triangle_cells)))


def test_boundary_matrix_from_multi_filtration(triangle_cells):
    assert_boundary_matrix(d.boundary(d.MultiFiltration(triangle_cells)))


def test_coboundary_matrix_from_multi_filtration(triangle_cells):
    assert_coboundary_matrix(d.coboundary(d.MultiFiltration(triangle_cells)))
