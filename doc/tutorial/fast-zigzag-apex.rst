.. _fast-zigzag-apex:

Fast Zigzag and Apex Representatives
-------------------------------------

.. testsetup::

    from __future__ import print_function   # if you are using Python 2
    import dionysus as d

Dionysus can `compute zigzag persistence via extended persistence <https://arxiv.org/abs/2204.11080>`_ and in the process recover `apex representatives <https://arxiv.org/abs/2502.17704>`_.
The filtration setup is the same as for :ref:`zigzag persistence <zigzag-persistence>`. One has to specify simplices and the times at which they appear and disappear:

.. doctest::

    >>> simplices = [[0], [1], [0,1], [2], [0,2], [1,2]]
    >>> times = [[.4,.6,.7,1.], [.1,.2,.3,1.], [.8,.95], [.5], [.8,1.], [.9,1.]]

Then the cone for extended persistence can be set up via :func:`~dionysus.fast_zigzag`, and its homology compute using :func:`~dionysus._dionysus.homology_persistence`.
If one wants to recover apex representatives, it is essential to compute matrix :math:`V` in :math:`R=DV` decomposition by passing ``method = 'matrix_v'``.

.. doctest::

    >>> cone = d.fast_zigzag(simplices, times)
    >>> r,v = d.homology_persistence(cone, method = 'matrix_v')

Then the persistence diagrams can be recovered using :func:`~dionysus.init_zigzag_diagrams`:

.. doctest::

    >>> dgms = d.init_zigzag_diagrams(r,cone,diagonal=True)

The diagrams are split by dimension, where each dimension stores a dictionary, with up to four entries, one for each extended persistence diagram type (ordinary, relative, and two types of extended), each of which corresponds to a different type of interval in zigzag persistence (``oo`` = open-open, ``co`` = closed-open, ``oc`` = open-closed, ``cc`` = closed-closed).
The points in these diagrams can be used to recover their corresponding apex representatives:

.. doctest::

    >>> max_t = max(max(t) for t in times)
    >>> for dim,type_dgm in enumerate(dgms):
    ...     print("Dimension:", dim)
    ...     for t,dgm in type_dgm.items():
    ...         print("Type:", t)
    ...         for pt in dgm:
    ...             print(pt)
    ...             apex_rep = d.apex(pt,r,v,cone)
    Dimension: 0
    Type: cc
    (0.1,0.2)
    (0.3,inf)
    (0.4,0.6)
    Type: co
    (0.5,0.8)
    (0.7,0.8)
    Type: oc
    (1,1)
    (1,1)
    Dimension: 1
    Type: cc
    (0.9,0.95)

The apex representatvie, an instance of :class:`~dionysus.zigzag.ApexRepresentative`, may be of independent interest:

.. doctest::

    >>> print("apex representative:", " + ".join(f"{cone[x]} × {time} ⋅ {c}" for (time, (x,c)) in apex_rep))
    apex representative: <0,1> 0.8 × (0.8999999761581421, 0.949999988079071) ⋅ 1 + <0,2> 0.8 × (0.8999999761581421, 0.949999988079071) ⋅ 1 + <1,2> 0.9 × (0.8999999761581421, 0.949999988079071) ⋅ 1

Given an apex representative, one can recover a (compatible) barcode representative at a given time:

.. doctest::

    >>> max_t = max(max(t) for t in times)

    >>> left = pt.birth
    >>> right = pt.death
    >>> if pt.death != float('inf'):
    ...     middle = (pt.birth + pt.death)/2
    ... else:
    ...     middle = max_t + 1
    ...     right = middle

    >>> left_representative = d.point_representative(apex_rep, left)
    >>> middle_representative = d.point_representative(apex_rep, middle)
    >>> right_representative = d.point_representative(apex_rep, right)

    >>> print(f"left ({left}) representative:", ' + '.join(f"{coeff} ⋅ {cone[idx]}" for (idx,coeff) in left_representative))
    left (0.8999999761581421) representative: 1 ⋅ <0,1> 0.8 + 1 ⋅ <0,2> 0.8 + 1 ⋅ <1,2> 0.9

    >>> print(f"midpoint ({middle}) representative:", ' + '.join(f"{coeff} ⋅ {cone[idx]}" for (idx,coeff) in middle_representative))
    midpoint (0.9249999821186066) representative: 1 ⋅ <0,1> 0.8 + 1 ⋅ <0,2> 0.8 + 1 ⋅ <1,2> 0.9

    >>> print(f"right ({right}) representative:", ' + '.join(f"{coeff} ⋅ {cone[idx]}" for (idx,coeff) in right_representative))
    right (0.949999988079071) representative: 1 ⋅ <0,1> 0.8 + 1 ⋅ <0,2> 0.8 + 1 ⋅ <1,2> 0.9
