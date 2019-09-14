.. _diagnostics:

Diagnostics
-----------

If you run into problems, the first thing to check is whether the filtration
you are giving to Dionysus is a simplicial complex, using
:func:`~dionysus.is_simplicial`:

.. testsetup::

    from __future__ import print_function   # if you are using Python 2
    import dionysus as d

.. doctest::

   >>> f = d.Filtration([([0], 1.), ([0,1], 2.)])
   >>> d.is_simplicial(f)
   False

You can also ask this function to report problems it finds:

.. doctest::

   >>> d.is_simplicial(f, report = True)
   <1> 0 in boundary of <0,1> 2 not found in the filtration
   False
