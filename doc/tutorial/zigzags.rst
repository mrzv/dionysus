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

    >>> zz, dgms = d.zigzag_homology_persistence(f, times)

The function returns a pair: an internal representation of
:class:`~dionysus._dionysus.ZigzagPersistence`, which stores cycles still alive
in the right-most homology group in the sequence, and the persistence diagrams
that represent the decomposition of the sequence.

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
    1*4 + 1*5 + 1*6
    1*0
