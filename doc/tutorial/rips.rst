Vietoris--Rips Complexes
------------------------

.. testsetup::

    from __future__ import print_function   # if you are using Python 2
    from dionysus import *

.. image:: figures/vietoris-rips.png
   :scale: 50 %
   :align: center

.. testsetup::

   import numpy as np
   np.random.seed(0)

.. doctest::
   :options: +ELLIPSIS

   >>> import numpy as np
   >>> points = np.random.random((100,2))
   >>> f = fill_rips(points, 2, .3)
   >>> print(f)
   Filtration with 5974 simplices
   >>> for s in f:
   ...     print(s)
   <0> 0
   ...
   <9,61,92> 0.0899135
   <9,72,92> 0.0899135
   <9,82,92> 0.0899135

:func:`~dionysus._dionysus.fill_rips` also accepts `condensed distance matrices <https://docs.scipy.org/doc/scipy-0.18.1/reference/spatial.distance.html>`_
(linearized lower triangular part of a symmetric matrix):

.. doctest::

   >>> from scipy.spatial.distance import pdist
   >>> dists = pdist(points)
   >>> f = fill_rips(dists, 2, .3)
   >>> print(f)
   Filtration with 5974 simplices

SciPy provides a helper function `squareform <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html>`_
to convert between redundant square matrices (:math:`n \times n`) and condensed
matrices (vectors with :math:`{n \choose 2}` elements).

.. doctest::

  >>> from scipy.spatial.distance import squareform
  >>> sq_dist = squareform(dists)
  >>> print(sq_dist.shape)
  (100, 100)
  >>> print(squareform(sq_dist).shape)
  (4950,)


