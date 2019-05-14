API
===

The following classes and functions live in the ``dionysus`` module.

Filtration
----------

.. autoclass:: dionysus._dionysus.Simplex
    :members:
    :special-members: __iter__, __getitem__, __len__, __contains__

.. autoclass:: dionysus._dionysus.Filtration
    :members:
    :special-members: __iter__, __getitem__, __len__, __contains__

.. autofunction:: dionysus._dionysus.fill_rips

.. autofunction:: dionysus.closure

.. autofunction:: dionysus._dionysus.fill_freudenthal


Persistence
-----------

.. autofunction:: dionysus._dionysus.homology_persistence

.. autofunction:: dionysus._dionysus.cohomology_persistence

.. autofunction:: dionysus._dionysus.omnifield_homology_persistence

.. autofunction:: dionysus._dionysus.zigzag_homology_persistence


Diagrams
--------

.. autofunction:: dionysus._dionysus.init_diagrams

.. autoclass:: dionysus._dionysus.Diagram
    :members:
    :special-members: __iter__, __len__

.. autoclass:: dionysus._dionysus.DiagramPoint
    :members:

.. autofunction:: dionysus._dionysus.wasserstein_distance

.. autofunction:: dionysus._dionysus.bottleneck_distance

Matrices
--------

.. autoclass:: dionysus._dionysus.Chain
    :members:
    :special-members: __iter__, __getitem__, __len__, __eq__, __ne__

.. autoclass:: dionysus._dionysus.ChainEntry
    :members:

.. autoclass:: dionysus._dionysus.ReducedMatrix
    :members:
    :special-members: __iter__, __getitem__, __len__

.. autoclass:: dionysus._dionysus.OmniFieldPersistence
    :members:
    :special-members: __len__

.. autoclass:: dionysus._dionysus.CoChain
    :members:
    :special-members: __iter__, __getitem__, __len__, __eq__, __ne__

.. autoclass:: dionysus._dionysus.CoChainEntry
    :members:

.. autoclass:: dionysus._dionysus.CohomologyPersistenceColumnHead
    :members:

    Wrapper around columns in `CohomologyPersistence`.

.. autoclass:: dionysus._dionysus.CohomologyPersistence
    :members:
    :special-members: __len__, __iter__

.. autoclass:: dionysus._dionysus.ZigzagPersistence
    :members:
    :special-members: __len__, __iter__

.. autoclass:: dionysus._dionysus.Zp


Plotting
--------

The following functions live in ``dionysus.plot`` module.

.. autofunction:: dionysus.plot.plot_diagram

.. autofunction:: dionysus.plot.plot_bars

.. autofunction:: dionysus.plot.plot_diagram_density


Diagnostics
-----------

.. autofunction:: dionysus.is_simplicial


.. Auto
   ----

   .. automodule:: dionysus._dionysus
       :members:

   .. automodule:: dionysus
       :members:
