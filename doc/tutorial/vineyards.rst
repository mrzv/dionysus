Vineyards
---------

Vineyards track how persistence pairings and persistence points change while a
filtration order changes. Dionysus exposes two related interfaces:

* ``VineyardV`` and ``VineyardU`` maintain a reduced boundary matrix under
  explicit adjacent transpositions.
* ``vineyard_linear_homotopy`` computes the adjacent transpositions induced by
  a straight-line interpolation between two filtration functions and records
  the resulting vines.

Adjacent transpositions
~~~~~~~~~~~~~~~~~~~~~~~

The low-level vineyard state is initialized from either a filtration or boundary
columns in stable cell ids. For explicit columns, the stable ids are the column
positions. For a simple edge, the two vertices have ids ``0`` and ``1``, and the
edge has id ``2``:

.. doctest::

    >>> import dionysus as d
    >>> columns = [[], [], [(1, 0), (1, 1)]]
    >>> v = d.Vineyard(columns, field=d.Zp(2))
    >>> isinstance(v, d.VineyardV)
    True
    >>> [v.pair(i) == v.unpaired for i in range(len(v))]
    [True, False, False]
    >>> (v.pair(1), v.pair(2))
    (2, 1)

The default ``Vineyard`` factory uses the ``matrix_v`` method. Pass
``method='matrix_u'`` to maintain trails instead:

.. doctest::

    >>> u = d.Vineyard(columns, field=d.Zp(2), method='matrix_u')
    >>> isinstance(u, d.VineyardU)
    True

When initialized from a :class:`~dionysus._dionysus.Filtration`, the stable ids
are the current filtration indices. Dionysus does not sort the filtration for
you; call :meth:`~dionysus._dionysus.Filtration.sort` first if you want the
usual data/dimension/lexicographic order.

.. doctest::

    >>> f = d.Filtration([[0], [1], [0, 1]])
    >>> vf = d.Vineyard(f, field=d.Zp(2))
    >>> [vf.cell_at(i) for i in range(len(vf))]
    [0, 1, 2]

Both vineyard states support adjacent swaps by current filtration position. The
state updates the reduced matrix, the ``V`` chains or ``U`` trails, and the
cached persistence pairing:

.. doctest::

    >>> v.transpose_position(0)
    (0, 1)
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
value at least as large as the values of its faces. When multiple simplices have
the same value, Dionysus breaks ties by dimension and then lexicographically so
each intermediate order remains a filtration.

.. doctest::

    >>> f = d.Filtration([[0], [1], [0, 1]])
    >>> result = d.vineyard_linear_homotopy(f,
    ...                                     [0.0, 1.0, 1.0],
    ...                                     [1.0, 0.0, 1.0],
    ...                                     field=d.Zp(2))
    >>> [(round(e.time, 1), e.first, e.second) for e in result.events]
    [(0.5, 0, 1)]
    >>> result.final_order
    [1, 0, 2]

The result stores:

* ``vineyard``: the final ``VineyardV`` or ``VineyardU`` state.
* ``events``: adjacent transpositions with local pairing information before and
  after the swap.
* ``vines``: piecewise-linear persistence vines. Each vine contains segments
  with endpoint times, birth/death values, birth/death cell ids, and the events
  that opened or closed the segment.
* ``final_order``: stable cell ids in the final filtration order.

Use ``method='matrix_u'`` to run the same linear homotopy with a MatrixU-backed
vineyard state:

.. doctest::

    >>> result = d.vineyard_linear_homotopy(f,
    ...                                     [0.0, 1.0, 1.0],
    ...                                     [1.0, 0.0, 1.0],
    ...                                     field=d.Zp(2),
    ...                                     method='matrix_u')
    >>> isinstance(result.vineyard, d.VineyardU)
    True
