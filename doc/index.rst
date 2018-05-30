.. Dionysus documentation master file, created by
   sphinx-quickstart on Tue Mar  7 10:16:22 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Dionysus 2 documentation!
====================================

.. sidebar:: Contents

    .. toctree::
       :maxdepth: 3

       tutorial/index
       API

Dionysus 2 is the second incarnation of the library for computing
persistent homology. As before, it's written in C++, with Python bindings.
The second version is re-written from scratch, which helps it accomplish a few goals:

  * `Modified BSD license <https://github.com/mrzv/dionysus/blob/master/LICENSE.txt>`_ (because GPL causes too many problems in academic software).
  * No dependency on Boost.Python; Dionysus 2 uses (and includes) `PyBind11 <https://github.com/pybind/pybind11>`_ instead, which greatly simplifies installation.
  * Cleaner, more consistent internal design (for example, all algorithms support arbitrary fields).
  * Some new algorithms, e.g., :ref:`omni-field` and Wasserstein and bottleneck :ref:`distance computation <diagram-distances>` from `Hera <https://bitbucket.org/grey_narn/hera>`_.
  * A few :ref:`plotting` routines, based on `Matplotlib <https://matplotlib.org/>`_.
  * Better integration with `NumPy <http://www.numpy.org/>`_.

Features that haven't (yet) made it over from `Dionysus 1 <http://mrzv.org/software/dionysus>`_ include vineyards.
Alpha shape filtrations are available via `DioDe <https://github.com/mrzv/diode>`_.

**Dependencies:**
  * `Boost <http://www.boost.org/>`_, although Dionysus 2 doesn't link any of its libraries, so it's considerably easier to build the project.
  * (Optional) `SciPy <https://www.scipy.org/>`_ for the LSQR routine used in :ref:`circular`.
  * (Optional) `Maplotlib <https://matplotlib.org/>`_ for :ref:`plotting`.

**Requirements:**
  * Boost needs to be at least version 1.55.
  * If you are using GCC, the oldest supported version is 5.4.

**Contact:**
  * please use the `dionysus mailing list <https://groups.io/g/dionysus/>`_
    for all questions and discussion related to the library;
  * GitHub's `issue tracker <https://github.com/mrzv/dionysus/issues>`_
    is a central location for bug reports and feature requests.

Get, Build, Install
-------------------

The simplest way to install Dionysus, as a Python package, is from `PyPI <https://pypi.org/project/dionysus/>`_:

.. parsed-literal::

    pip install --verbose dionysus

Pass ``--upgrade`` to ``pip``, if you have already installed some version of Dionysus.

Alternatively, you can install it directly from the development repository (this gives you the latest version):

.. parsed-literal::

    pip install --verbose `git+https://github.com/mrzv/dionysus.git <https://github.com/mrzv/dionysus.git>`_

Alternatively, you can clone and build everything by hand.
To get Dionysus 2, either clone its `repository <https://github.com/mrzv/dionysus>`_:

.. parsed-literal::

    git clone `<https://github.com/mrzv/dionysus.git>`_

or download it as a `Zip archive <https://github.com/mrzv/dionysus/archive/master.zip>`_.

To build the project::

    mkdir build
    cd build
    cmake ..
    make

To use the Python bindings, either launch Python from ``.../build/bindings/python`` or add this directory to your ``PYTHONPATH`` variable, by adding::

    export PYTHONPATH=.../build/bindings/python:$PYTHONPATH

to your ``~/.bashrc`` or ``~/.zshrc``.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
