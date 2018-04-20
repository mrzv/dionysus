.. _omni-field:

Omni-field Persistence
----------------------

.. testsetup::

    from __future__ import print_function   # if you are using Python 2
    import dionysus as d

It's possible to compute persistence over all fields :math:`\mathbb{Z}_p` at once, by doing
some extra bookkeeping. :class:`~dionysus._dionysus.OmniFieldPersistence`
performs the necessary operations. We start with a triangulation of the
Klein bottle:

.. doctest::

    >>> klein_bottle = [[0], [1], [2], [3], [4], [5], [6], [7], [8],
    ...       [0,1], [1,2], [2,0], [0,3], [3,4], [4,0],
    ...       [1,5], [5,6], [6,2], [2,7], [7,8], [8,1],
    ...       [3,5], [5,7], [7,3], [4,6], [6,8], [8,4],
    ...       [0,5], [1,7], [2,3], [3,6], [5,8], [7,4],
    ...       [4,2], [6,1], [8,0],
    ...       [0,3,5], [0,1,5], [1,5,7], [1,2,7], [2,7,3], [2,3,0],
    ...       [3,4,6], [3,5,6], [5,6,8], [5,7,8], [7,8,4], [7,3,4],
    ...       [4,0,2], [4,6,2], [6,2,1], [6,8,1], [8,1,0], [8,4,0]]
    >>> f = d.Filtration(klein_bottle)
    >>> print(f)
    Filtration with 54 simplices

We can compute persistence over all fields :math:`\mathbb{Z}_p` at once:

.. doctest::

    >>> ofp = d.omnifield_homology_persistence(f)

We can examine which primes produce special results. This does not necessarily
mean that persistence diagram over the respective finite fields looks
different, but it does indicate that something special happened during the
reduction. In contrast, for any prime not in this list, the persistence
diagrams look the same.

.. doctest::

    >>> print(ofp.primes())
    [2]


To access the results, we can construct persistence diagrams over a particular
:math:`\mathbb{Z}_p`. E.g., over :math:`\mathbb{Z}_2`, we get:

.. doctest::

    >>> dgms = d.init_diagrams(ofp, f, 2)
    >>> for i, dgm in enumerate(dgms):
    ...     print("Dimension:", i)
    ...     for pt in dgm:
    ...         print(pt)
    Dimension: 0
    (0,inf)
    Dimension: 1
    (0,inf)
    (0,inf)
    Dimension: 2
    (0,inf)

Indicating that :math:`H_0(K, \mathbb{Z}_2) = \mathbb{Z}_2, H_1(K, \mathbb{Z}_2) = \mathbb{Z}_2 \oplus \mathbb{Z}_2, H_2(K, \mathbb{Z}_2) = \mathbb{Z}_2`.

But over :math:`\mathbb{Z}_3` (as well as over any other :math:`\mathbb{Z}_p`), we get:

.. doctest::

    >>> dgms = d.init_diagrams(ofp, f, 3)
    >>> for i, dgm in enumerate(dgms):
    ...     print("Dimension:", i)
    ...     for pt in dgm:
    ...         print(pt)
    Dimension: 0
    (0,inf)
    Dimension: 1
    (0,inf)
    Dimension: 2

Indicating that :math:`H_0(K, \mathbb{Z}_3) = \mathbb{Z}_3, H_1(K, \mathbb{Z}_3) = \mathbb{Z}_3, H_2(K, \mathbb{Z}_3) = 0`.

It's possible to examine the columns of the reduced matrices over any given prime:

.. doctest::
   :options: +NORMALIZE_WHITESPACE

    >>> for i in range(len(ofp)):
    ...     if ofp.special(i, 2):
    ...         print("Column %d mod %d:" % (i,2), ofp.column(i, 2))
    ...         print("Column %d mod %d:" % (i,3), ofp.column(i, 3))
    Column 53 mod 2:
    Column 53 mod 3: 1*9 + 1*10 + 2*11
