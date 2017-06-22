Cohomology Persistence
----------------------

.. testsetup::

    from __future__ import print_function   # if you are using Python 2
    from dionysus import *

.. doctest::

    >>> simplices = [([2], 4), ([1,2], 5), ([0,2], 6),
    ...              ([0], 1),   ([1], 2), ([0,1], 3)]
    >>> f = Filtration()
    >>> for vertices, time in simplices:
    ...     f.append(Simplex(vertices, time))
    >>> f.sort()

Applying cohomology functor to the filtration, we get a sequence of cohomology groups, connected by linear maps:
:math:`H^*(K_1) \to H^*(K_2) \to \ldots \to H^*(K_n)`. To compute decomposition of this sequence, i.e., persistence barcode,
we use :func:`~dionysus._dionysus.cohomology_persistence`.

.. doctest::

    >>> p = cohomology_persistence(f, prime=2)

The returned object stores the persistence pairs as well as the cocycles still
alive at the end of the filtration (i.e., a basis for :math:`H^*(K_n)`). To
extract persistence diagrams, we use, as before,
:func:`~dionysus._dionysus.init_diagrams`:

.. doctest::

    >>> dgms = init_diagrams(p, f)
    >>> for i,dgm in enumerate(dgms):
    ...     print(i)
    ...     for pt in dgm:
    ...         print(pt)
    0
    (1,inf)
    (2,3)
    (4,5)
    1
    (6,inf)

To access the alive cocycles, we iterate over the returned object. For each
element, `index` stores the index in the filtration when the cocycle was born,
while `cocycle` stores the cocycle itself.

.. doctest::

    >>> for c in p:
    ...     print(c.index, c.cocycle)
    0 1*0 + 1*1 + 1*3
    5 1*5


