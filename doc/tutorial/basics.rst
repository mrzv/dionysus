First, we import everything from Dionysus:

.. doctest::

    >>> from __future__ import print_function   # if you are using Python 2
    >>> from dionysus import *


.. todo::

    Explain where to get the sample data from.

Simplices
---------

.. doctest::

    >>> s = Simplex([0,1,2])
    >>> print("Dimension:", s.dimension())
    Dimension: 2

We can iterate over the vertices of the simplex.

.. doctest::

    >>> for v in s:
    ...     print(v)
    0
    1
    2

Or over the boundary:

.. doctest::

    >>> for sb in s.boundary():
    ...     print(sb)
    <1,2> 0
    <0,2> 0
    <0,1> 0

Simplices can store optional data, and the 0 reported after each boundary edge is the default value of the data.

We can use :func:`~dionysus.closure` to generate all faces of a set of
simplices. For example, an 8-sphere is the 8-dimensional skeleton of the
closure of the 9-simplex.

.. doctest::

    >>> simplex9 = Simplex([0,1,2,3,4,5,6,7,8,9])
    >>> sphere8  = closure([simplex9], 8)
    >>> print(len(sphere8))
    1022

Filtration
----------

A filtration is a nested sequence of simplicial complexes,
:math:`K_1 \subseteq K_2 \subseteq \ldots \subseteq K_n`.
Without loss of generality, we can assume that
two consecutive filtrations differ by a single simplex, so we can think of
a filtration as a sequence of simplices.

.. image:: figures/filtration.png

In Dionysus, a filtration is represented by a special class
:class:`~dionysus._dionysus.Filtration` that supports both iterating over the
simplices and looking up an index given a simplex. A filtration can be
:meth:`~dionysus._dionysus.Filtration.sort`\ ed. By default this orders
simplices by their data, breaking ties by dimension, and then
lexicographically.

.. doctest::

    >>> simplices = [([2], 4), ([1,2], 5), ([0,2], 6),
    ...              ([0], 1),   ([1], 2), ([0,1], 3)]
    >>> f = Filtration()
    >>> for vertices, time in simplices:
    ...     f.append(Simplex(vertices, time))
    >>> f.sort()
    >>> for s in f:
    ...    print(s)
    <0> 1
    <1> 2
    <0,1> 3
    <2> 4
    <1,2> 5
    <0,2> 6

We can lookup the index of a given simplex. (Indexing starts from 0.)

.. doctest::

    >>> print(f.index(Simplex([1,2])))
    4

Persistent Homology
-------------------

Applying homology functor to the filtration, we get a sequence of homology groups, connected by linear maps:
:math:`H_*(K_1) \to H_*(K_2) \to \ldots \to H_*(K_n)`. To compute decomposition of this sequence, i.e., persistence barcode,
we use :func:`~dionysus._dionysus.homology_persistence`.

.. doctest::
   :options: +NORMALIZE_WHITESPACE

    >>> m = homology_persistence(f)
    >>> for i,c in enumerate(m):
    ...     print(i, c)
    0
    1
    2 1*0 + 1*1
    3
    4 1*1 + 1*3
    5

.. image:: figures/barcode.png

.. doctest::

    >>> for i in range(len(m)):
    ...     if m.pair(i) < i: continue      # skip negative simplices
    ...     dim = f[i].dimension()
    ...     if m.pair(i) != m.unpaired:
    ...         print(dim, i, m.pair(i))
    ...     else:
    ...         print(dim, i)
    0 0
    0 1 2
    0 3 4
    1 5

Alternatively:

.. doctest::

    >>> m = homology_persistence(f, method = 'column')


**Homology.**
Dionysus doesn’t compute homology directly, but we can get it as a by-product
of persistent homology.

.. doctest::

    >>> f = Filtration(sphere8)
    >>> f.sort()
    >>> m = homology_persistence(f, prime=2)
    >>> dgms = init_diagrams(m, f)
    >>> for i, dgm in enumerate(dgms):
    ...     print("Dimension:", i)
    ...     for p in dgm:
    ...         print(p)
    Dimension: 0
    (0,inf)
    Dimension: 1
    Dimension: 2
    Dimension: 3
    Dimension: 4
    Dimension: 5
    Dimension: 6
    Dimension: 7
    Dimension: 8
    (0,inf)


Diagram Distances
-----------------

:func:`~dionysus._dionysus.wasserstein_distance` can compute `q`-th Wasserstein distance between a pair of persistence diagrams.

.. testsetup::

    import numpy as np
    np.random.seed(0)

.. doctest::

    >>> f1 = fill_rips(np.random.random((20, 2)), 2, 1)
    >>> m1 = homology_persistence(f1)
    >>> dgms1 = init_diagrams(m1, f1)
    >>> f2 = fill_rips(np.random.random((20, 2)), 2, 1)
    >>> m2 = homology_persistence(f2)
    >>> dgms2 = init_diagrams(m2, f2)
    >>> dist = wasserstein_distance(dgms1[1], dgms2[1], q=2)
    >>> print("Distance between 1-dimensional persistence diagrams:", dist)
    Distance between 1-dimensional persistence diagrams: 0.037361082536
