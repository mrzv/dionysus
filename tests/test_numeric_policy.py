import math
import struct

import dionysus as d


def float32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def test_simplex_data_uses_float32_precision():
    value = 1.123456789

    simplex = d.Simplex([0], value)

    assert simplex.data == float32(value)
    assert simplex.data != value


def test_diagram_coordinates_use_simplex_data_precision():
    birth = 1.123456789
    death = 2.123456789
    filtration = d.Filtration([
        ([0], birth),
        ([1], death),
        ([0, 1], death),
    ])

    diagrams = d.init_diagrams(d.homology_persistence(filtration), filtration)

    assert [(pt.birth, pt.death) for pt in diagrams[0]] == [(float32(birth), math.inf)]
