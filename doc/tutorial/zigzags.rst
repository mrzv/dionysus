Zigzag Persistence
------------------

.. testsetup::

    from __future__ import print_function   # if you are using Python 2
    import dionysus as d

Carlsson and de Silva introduced `zigzag persistence
<https://arxiv.org/abs/0812.0197>`_, a generalization of ordinary persistent
homology that allows linear maps connecting homology groups to point in either direction in the sequence, e.g.:
:math:`H_*(K_1) \to H_*(K_2) \leftarrow H_*(K_3) \to H_*(K_4) \leftarrow \ldots`
To express such a zigzag filtration, we consider the maximal simplicial
complex, :math:`\cup K_i`, and encode it as
a :class:`~dionysus._dionysus.Filtration`:

.. doctest::

    >>> f = d.Filtration([[0], [1], [0,1], [2], [0,2], [1,2]])

For each simplex in the complex, we specify a list of times when it enters and
leaves the filtration. This information is provided as a list of lists,
``times``. For the `i`-th simplex in the filtration, ``times[i]`` is a list of
times, where values in even positions (counting from 0) specify when the
simplex is added to the complex and odd positions when it is removed:

.. doctest::

    >>> times = [[.4, .6, .7], [.1], [.9], [.9], [.9], [.9]]

Given the two inputs, we can compute zigzag persistent homology
of the corresponding sequence of simplicial complexes, using
:func:`~dionysus._dionysus.zigzag_homology_persistence`:

.. doctest::

    >>> zz, dgms, cells = d.zigzag_homology_persistence(f, times)

The function returns a triple: an internal representation of
:class:`~dionysus._dionysus.ZigzagPersistence`, which stores cycles still alive
in the right-most homology group in the sequence, the persistence diagrams that
represent the decomposition of the sequence, auxiliary map to translate from
internal indices used in the cycles into the indices of the simplices in the
:class:`~dionysus._dionysus.Filtration`:

.. doctest::

    >>> print(zz)
    Zigzag persistence with 2 alive cycles

    >>> for i,dgm in enumerate(dgms):
    ...     print("Dimension:", i)
    ...     for p in dgm:
    ...         print(p)
    Dimension: 0
    (0.4,0.6)
    (0.7,0.9)
    (0.1,inf)
    Dimension: 1
    (0.9,inf)

    >>> for z in zz:
    ...     print(z)
    1*0
    1*4 + 1*5 + 1*6

    >>> for x in sorted(cells):
    ...     print(x)
    (0, 1)
    (2, 0)
    (3, 3)
    (4, 2)
    (5, 4)
    (6, 5)


Representative cycles
~~~~~~~~~~~~~~~~~~~~~

The first and the third element of the triple, combined, can be used to extract
representative cycles. The third element is the map from the cycle's internal
representation to the filtration indices. The following snippet outputs the
cycles in terms of the simplices.

.. doctest::

    >>> for z in zz:
    ...     print(' + '.join("%d * (%s)" % (x.element, f[cells[x.index]]) for x in z))
    1 * (<1> 0)
    1 * (<0,1> 0) + 1 * (<0,2> 0) + 1 * (<1,2> 0)

Intermediate steps
~~~~~~~~~~~~~~~~~~

:func:`~dionysus._dionysus.zigzag_homology_persistence` takes an optional `callback` argument,
which gets called back after every step of the zigzag. The function receives four arguments, `(i,t,d,zz)`.
`i` is the index of the simplex being added or removed. `t` is the current
time. `d` is the direction: ``True`` if the simplex is being added, ``False``,
if removed. `zz` is the current state of :class:`~dionysus._dionysus.ZigzagPersistence`.

.. doctest::

    >>> def detail(i,t,d,zz,cells):
    ...     print(i,t,d)
    ...     for z in zz:
    ...         print(z, ' -> ', ' + '.join("%d * (%s)" % (x.element, f[cells[x.index]]) for x in z))

    >>> zz, dgms, cells = d.zigzag_homology_persistence(f, times, callback = detail)
    1 0.10000000149011612 True
    1*0  ->  1 * (<1> 0)
    0 0.4000000059604645 True
    1*0  ->  1 * (<1> 0)
    1*1  ->  1 * (<0> 0)
    0 0.6000000238418579 False
    1*0  ->  1 * (<1> 0)
    0 0.699999988079071 True
    1*0  ->  1 * (<1> 0)
    1*2  ->  1 * (<0> 0)
    3 0.8999999761581421 True
    1*0  ->  1 * (<1> 0)
    1*2  ->  1 * (<0> 0)
    1*3  ->  1 * (<2> 0)
    2 0.8999999761581421 True
    1*0  ->  1 * (<1> 0)
    1*3  ->  1 * (<2> 0)
    4 0.8999999761581421 True
    1*0  ->  1 * (<1> 0)
    5 0.8999999761581421 True
    1*0  ->  1 * (<1> 0)
    1*4 + 1*5 + 1*6  ->  1 * (<0,1> 0) + 1 * (<0,2> 0) + 1 * (<1,2> 0)
