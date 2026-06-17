.. _vineyards:

Vineyards
---------

Vineyards track how persistence pairings and diagram points change as
filtration order changes. Dionysus exposes two related interfaces:

* ``VineyardV`` and ``VineyardU`` maintain a reduced boundary matrix under
  explicit adjacent transpositions.
* ``vineyard_linear_homotopy`` computes the adjacent transpositions induced by
  a straight-line interpolation between two filtration functions and records
  the resulting vines.

Adjacent transpositions
~~~~~~~~~~~~~~~~~~~~~~~

The low-level vineyard state can be initialized directly from a filtration. The
stable ids are the current filtration indices. For a simple edge, the two
vertices have ids ``0`` and ``1``, and the edge has id ``2``:

.. doctest::

    >>> import dionysus as d
    >>> f = d.Filtration([[0], [1], [0, 1]])
    >>> v = d.Vineyard(f, field=d.Zp(2))
    >>> isinstance(v, d.VineyardV)
    True
    >>> [v.pair(i) == v.unpaired for i in range(len(v))]
    [True, False, False]
    >>> (v.pair(1), v.pair(2))
    (2, 1)

The default ``Vineyard`` factory uses the ``matrix_v`` method. Pass
``method='matrix_u'`` to maintain trails instead:

.. doctest::

    >>> u = d.Vineyard(f, field=d.Zp(2), method='matrix_u')
    >>> isinstance(u, d.VineyardU)
    True

The constructor does not sort the filtration for you; call
:meth:`~dionysus._dionysus.Filtration.sort` first if you want the usual
data/dimension/lexicographic order.

.. doctest::

    >>> [v.cell_at(i) for i in range(len(v))]
    [0, 1, 2]

For lower-level uses, ``Vineyard`` also accepts explicit boundary columns in
stable cell ids, e.g. ``d.Vineyard([[], [], [(1, 0), (1, 1)]], field=d.Zp(2))``.

Both vineyard states support adjacent swaps by current filtration position. The
state updates the reduced matrix, the ``V`` chains or ``U`` trails, and the
cached persistence pairing:

.. doctest::

    >>> v.transpose(0)
    (0, 1, True)
    >>> [v.cell_at(i) for i in range(len(v))]
    [1, 0, 2]

The stable ids do not change during transpositions. Methods such as ``pair()``,
``low()``, ``pivot()``, ``reduced_column()``, ``chain()``, and ``trail()`` use
those stable ids.

Linear homotopy
~~~~~~~~~~~~~~~

``vineyard_linear_homotopy`` takes a filtration and two scalar functions on its
simplices. It follows the straight-line interpolation

.. math::

   f_t(\sigma) = (1-t) f_0(\sigma) + t f_1(\sigma), \quad 0 \leq t \leq 1,

performs the adjacent transpositions where neighboring simplices exchange
order, and records both the events and the persistence vines.

Both endpoint functions must be valid filtrations: every simplex must have a
value at least as large as the values of its faces. This validation is strict.
When multiple simplices have the same value, Dionysus breaks ties by dimension
and then lexicographically so each endpoint order remains a filtration. During
event scheduling, proposed adjacent swaps that would put a coface before one of
its faces are ignored.

.. doctest::

    >>> f = d.Filtration([[0], [1], [0, 1]])
    >>> result = d.vineyard_linear_homotopy(f,
    ...                                     [0.0, 1.0, 2.0],
    ...                                     [1.0, 0.0, 2.0],
    ...                                     field=d.Zp(2))
    >>> [(round(e.time, 1), e.first, e.second, e.pairing_switched) for e in result.events]
    [(0.5, 0, 1, True)]
    >>> result.final_order
    [1, 0, 2]

The result stores:

* ``vineyard``: the final ``VineyardV`` or ``VineyardU`` state.
* ``events``: adjacent transpositions with local pairing information before and
  after the swap, plus ``pairing_switched`` to indicate whether the persistence
  pairing changed.
* ``vines``: piecewise-linear persistence vines. Each vine contains segments
  with endpoint times, birth/death values, birth/death cell ids, and the events
  that opened or closed the segment.
* ``final_order``: stable cell ids in the final filtration order.

Use ``method='matrix_u'`` to run the same linear homotopy with a MatrixU-backed
vineyard state:

.. doctest::

    >>> result = d.vineyard_linear_homotopy(f,
    ...                                     [0.0, 1.0, 2.0],
    ...                                     [1.0, 0.0, 2.0],
    ...                                     field=d.Zp(2),
    ...                                     method='matrix_u')
    >>> isinstance(result.vineyard, d.VineyardU)
    True

The vines are stored as piecewise-linear segments. Each segment records its time
interval, birth and death values at the interval endpoints, and the stable cell
ids that define the feature:

Vine identity is continued combinatorially. If an adjacent transposition does
not change the persistence pairing, existing vine segments remain open and the
event does not split them. If the pairing switches, affected feature keys are
migrated by swapping the two transposed stable cell ids, preserving the vine id.
This also treats contact with the diagonal as continuation of the same vine;
segments are only split when the defining birth/death cells change.

.. doctest::

    >>> result = d.vineyard_linear_homotopy(f,
    ...                                     [0.0, 1.0, 2.0],
    ...                                     [1.0, 0.0, 2.0],
    ...                                     field=d.Zp(2))
    >>> unpaired = result.vineyard.unpaired
    >>> for i, vine in enumerate(result.vines):
    ...     for segment in vine.segments:
    ...         if segment.death_cell == unpaired:
    ...             death_cell = 'inf'
    ...             endpoint0 = (round(segment.t0, 1), round(segment.birth0, 1), 'inf')
    ...             endpoint1 = (round(segment.t1, 1), round(segment.birth1, 1), 'inf')
    ...         else:
    ...             death_cell = segment.death_cell
    ...             endpoint0 = (round(segment.t0, 1), round(segment.birth0, 1), round(segment.death0, 1))
    ...             endpoint1 = (round(segment.t1, 1), round(segment.birth1, 1), round(segment.death1, 1))
    ...         print(i, segment.birth_cell, death_cell, endpoint0, endpoint1)
    0 0 inf (0.0, 0.0, 'inf') (0.5, 0.5, 'inf')
    0 1 inf (0.5, 0.5, 'inf') (1.0, 0.0, 'inf')
    1 1 2 (0.0, 1.0, 2.0) (0.5, 0.5, 2.0)
    1 0 2 (0.5, 0.5, 2.0) (1.0, 1.0, 2.0)
