import math

import dionysus as d
import pytest


PRIME = 5


def order_after_events(result):
    if result.events:
        return result.events[-1].order
    return result.final_order


def test_homotopy_records_single_crossing_and_vines():
    filtration = d.Filtration([[0], [1]])

    result = d.vineyard_homotopy(
        filtration, [0.0, 1.0], [1.0, 0.0], field=d.Zp(PRIME)
    )

    assert isinstance(result.vineyard, d.VineyardV)
    assert len(result.events) == 1
    assert result.events[0].time == pytest.approx(0.5)
    assert (result.events[0].first, result.events[0].second) == (0, 1)
    assert result.final_order == [1, 0]
    assert order_after_events(result) == [1, 0]

    segments = [segment for vine in result.vines for segment in vine.segments]
    assert len(segments) == 4
    assert {segment.event1 for segment in segments if segment.t1 == pytest.approx(0.5)} == {0}
    assert {segment.event0 for segment in segments if segment.t0 == pytest.approx(0.5)} == {0}
    assert all(math.isinf(segment.death0) and math.isinf(segment.death1) for segment in segments)


def test_homotopy_processes_simultaneous_degenerate_block():
    filtration = d.Filtration([[0], [1], [2]])

    result = d.vineyard_homotopy(
        filtration, [0.0, 1.0, 2.0], [2.0, 1.0, 0.0], field=d.Zp(PRIME)
    )

    assert [event.time for event in result.events] == [pytest.approx(0.5)] * 3
    assert [(event.first, event.second) for event in result.events] == [
        (1, 2),
        (0, 2),
        (0, 1),
    ]
    assert result.final_order == [2, 1, 0]


def test_homotopy_preserves_filtration_order_with_dimension_ties():
    filtration = d.Filtration([[0], [1], [0, 1]])

    result = d.vineyard_homotopy(
        filtration, [0.0, 1.0, 1.0], [1.0, 0.0, 1.0], field=d.Zp(PRIME)
    )

    assert len(result.events) == 1
    assert result.final_order == [1, 0, 2]
    assert result.vineyard.position(0) < result.vineyard.position(2)
    assert result.vineyard.position(1) < result.vineyard.position(2)


def test_homotopy_supports_matrix_u():
    filtration = d.Filtration([[0], [1], [0, 1]])

    result = d.vineyard_homotopy(
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


def test_homotopy_rejects_non_filtration_endpoint_values():
    filtration = d.Filtration([[0], [1], [0, 1]])

    with pytest.raises(ValueError, match="not a filtration"):
        d.vineyard_homotopy(
            filtration, [0.0, 0.0, 1.0], [0.0, 1.0, 0.5], field=d.Zp(PRIME)
        )


def test_homotopy_rejects_unknown_method():
    filtration = d.Filtration([[0]])

    with pytest.raises(ValueError, match="unknown vineyard method"):
        d.vineyard_homotopy(filtration, [0.0], [1.0], method="unknown")
