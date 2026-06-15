import math

import dionysus as d


def diagram_points(diagrams):
    return [
        [(point.birth, point.death, point.data) for point in diagram]
        for diagram in diagrams
    ]


def chain_entries(chain):
    return [(entry.element, entry.index) for entry in chain]


def edge_filtration():
    return d.Filtration([[0], [1], [0, 1]])


def test_cohomology_persistence_matches_homology_pairs_and_diagrams(prime):
    filtration = edge_filtration()
    homology = d.homology_persistence(filtration, prime=prime)
    cohomology = d.cohomology_persistence(
        filtration, prime=prime, keep_cocycles=True
    )

    assert len(cohomology) == len(homology)
    assert [cohomology.pair(i) for i in range(len(cohomology))] == [
        homology.pair(i) for i in range(len(homology))
    ]
    assert diagram_points(d.init_diagrams(cohomology, filtration)) == [
        [(0.0, math.inf, 0)],
        [],
    ]
    assert diagram_points(d.init_diagrams(cohomology, filtration)) == diagram_points(
        d.init_diagrams(homology, filtration)
    )


def test_cohomology_persistence_exposes_alive_cocycle(prime):
    cohomology = d.cohomology_persistence(
        edge_filtration(), prime=prime, keep_cocycles=True
    )
    alive = list(cohomology)

    assert len(alive) == 1
    assert alive[0].index == 0
    assert chain_entries(alive[0].cocycle) == [(1, 0), (1, 1)]
    assert chain_entries(cohomology.cocycle(0)) == [(1, 0), (1, 1)]


def test_omnifield_persistence_basic_columns_and_diagrams():
    filtration = edge_filtration()
    omnifield = d.omnifield_homology_persistence(filtration)

    assert len(omnifield) == len(filtration)
    assert omnifield.primes() == []
    assert omnifield.specials() == {}
    assert [chain_entries(omnifield.column(i, 2)) for i in range(len(omnifield))] == [
        [],
        [],
        [(1, 0), (1, 1)],
    ]
    assert [chain_entries(omnifield.column(i, 5)) for i in range(len(omnifield))] == [
        [],
        [],
        [(4, 0), (1, 1)],
    ]
    assert diagram_points(d.init_diagrams(omnifield, filtration, 2)) == [
        [(0.0, math.inf, 0)],
        [],
    ]
