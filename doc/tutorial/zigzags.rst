Zigzag Persistence
------------------

.. testsetup::

    from __future__ import print_function   # if you are using Python 2
    from dionysus import *

.. doctest::

    >>> f = Filtration([[0], [1], [0,1], [2], [0,2], [1,2]])
    >>> times = [[.4, .6, .7], [.1], [.9], [.9], [.9], [.9]]
    >>> zz, dgms = zigzag_homology_persistence(f, times)

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
