import dionysus as tpl
import pytest


def test_issue72_zigzag_homology_persistence():
    f_list = [[0], [1], [0, 1], [2], [1, 2]]
    f = tpl.Filtration(f_list)
    times = [[0], [1], [1, 3, 4], [2], [2]]

    def detail(i, t, d, zz, cells):
        if d:
            action = 'added'
        else:
            action = 'removed'
        print(f'{i}) time= {t}, simplex:{action}')
        for z in zz:
            print(z, ' -> ', ' + '.join("%d * (%s)" % (x.element, f[cells[x.index]]) for x in z))

    zz, dgms, cells = tpl.zigzag_homology_persistence(f, times, callback=detail)

    assert len(zz) > 0
    assert len(dgms) > 0
    assert len(cells) == len(f)


def test_issue72_fast_zigzag_matrix_v():
    f_list = [[0], [1], [0, 1], [2], [1, 2]]
    times = [[0], [1], [1, 3, 4], [2], [2]]
    cone = tpl.fast_zigzag(f_list, times)
    r, v = tpl.homology_persistence(cone, method='matrix_v')

    assert len(cone) > 0
    assert len(r) == len(cone)
    assert len(v) == len(cone)


def test_zigzag_homology_persistence_rejects_mismatched_times_length():
    f = tpl.Filtration([[0], [1]])

    with pytest.raises(ValueError, match="times length"):
        tpl.zigzag_homology_persistence(f, [[0]])


def test_zigzag_homology_persistence_rejects_out_of_order_times():
    f = tpl.Filtration([[0]])

    with pytest.raises(ValueError, match="strictly increasing"):
        tpl.zigzag_homology_persistence(f, [[2, 1]])


def test_zigzag_homology_persistence_rejects_equal_remove_readd_time():
    f = tpl.Filtration([[0]])

    with pytest.raises(ValueError, match="strictly increasing"):
        tpl.zigzag_homology_persistence(f, [[0, 1, 1, 2]])


def test_zigzag_homology_persistence_requires_active_boundary():
    f = tpl.Filtration([[0], [1], [0, 1]])

    with pytest.raises(ValueError, match="boundary must be active"):
        tpl.zigzag_homology_persistence(f, [[0, 1], [0], [2]])


def test_zigzag_homology_persistence_rejects_removing_active_face():
    f = tpl.Filtration([[0], [1], [0, 1]])

    with pytest.raises(ValueError, match="coface is active"):
        tpl.zigzag_homology_persistence(f, [[0, 1], [0], [1]])
