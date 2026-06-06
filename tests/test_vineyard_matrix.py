import dionysus as d
import pytest


def vineyard_matrix():
    vm = d.VineyardMatrix(d.Zp(5), 4)
    vm.set_column(0, [])
    vm.set_column(1, [(1, 0)])
    vm.set_column(2, [])
    vm.set_column(3, [(2, 2), (1, 1)])
    return vm


def test_columns_are_stored_by_stable_id():
    vm = vineyard_matrix()

    assert vm[3] == [(1, 1), (2, 2)]
    assert vm.contains(3, 1)
    assert vm.contains(3, 2)
    assert not vm.contains(3, 0)
    assert not vm.contains(3, 4)


def test_set_column_rejects_out_of_range_entries():
    vm = d.VineyardMatrix(d.Zp(5), 4)

    with pytest.raises(IndexError):
        vm.set_column(0, [(1, 4)])


def test_transpose_updates_order_without_reordering_columns():
    vm = vineyard_matrix()

    swapped = vm.transpose_position(1)

    assert swapped == (1, 2)
    assert vm.cell_at(1) == 2
    assert vm.cell_at(2) == 1
    assert vm.position(1) == 2
    assert vm.position(2) == 1
    assert vm[3] == [(1, 1), (2, 2)]


def test_transpose_does_not_eagerly_refresh_lows():
    vm = vineyard_matrix()

    assert vm.low(3) == 2
    assert vm.pivot(2) == 3
    assert vm.pivot(1) == vm.unpaired

    vm.transpose_position(1)

    assert vm.low(3) == 2
    assert vm.pivot(2) == 3
    assert vm.pivot(1) == vm.unpaired


def test_refresh_low_updates_pivot_cache_for_one_column():
    vm = vineyard_matrix()

    vm.transpose_position(1)
    vm.refresh_low(3)

    assert vm.low(3) == 1
    assert vm.pivot(1) == 3
    assert vm.pivot(2) == vm.unpaired


def test_set_low_updates_pivot_cache_explicitly():
    vm = vineyard_matrix()

    vm.set_low(3, 1)

    assert vm.low(3) == 1
    assert vm.pivot(1) == 3
    assert vm.pivot(2) == vm.unpaired


def test_duplicate_lows_are_not_represented_in_global_pivot_map():
    vm = d.VineyardMatrix(d.Zp(5), 4)
    vm.set_column(2, [(1, 1)])
    vm.set_column(3, [(1, 1), (1, 2)])

    vm.transpose_position(1)
    vm.refresh_low(3)

    assert vm.low(2) == 1
    assert vm.low(3) == 1
    assert vm.pivot(1) == 3
