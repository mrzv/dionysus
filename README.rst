Dionysus 2
==========

Dionysus is a computational topology package focused on
persistent homology. It is written in C++, with Python bindings.
The second version (`previous version <http://mrzv.org/software/dionysus/>`_)
is re-written from scratch, which helps it accomplish a few goals:

  * `Modified BSD license <https://github.com/mrzv/dionysus/blob/master/LICENSE.txt>`_ (because GPL causes too many problems in academic software).
  * No dependency on Boost.Python; Dionysus 2 uses (and includes) `PyBind11 <https://github.com/pybind/pybind11>`_ instead, which greatly simplifies installation.
  * Cleaner, more consistent internal design (for example, all algorithms support arbitrary fields).
  * Some new algorithms, e.g., `omni-field persistence <http://mrzv.org/software/dionysus2/tutorial/omni-field.html#omni-field>`_ and Wasserstein and bottleneck `distance computation <http://mrzv.org/software/dionysus2/tutorial/basics.html#diagram-distances>`_ from `Hera <https://bitbucket.org/grey_narn/hera>`_.
  * A few `plotting <http://mrzv.org/software/dionysus2/tutorial/plotting.html#plotting>`_ routines, based on `Matplotlib <https://matplotlib.org/>`_.
  * Better integration with `NumPy <http://www.numpy.org/>`_.

Features that haven't (yet) made it over from `Dionysus 1 <http://mrzv.org/software/dionysus>`_ include vineyards.
Alpha shape filtrations are available via `DioDe <https://github.com/mrzv/diode>`_.

**Dependencies:**
  * `Boost <http://www.boost.org/>`_, although Dionysus 2 doesn't link any of its libraries, so it's considerably easier to build the project.
  * (Optional) `SciPy <https://www.scipy.org/>`_ for the LSQR routine used in `circular coordinates <http://mrzv.org/software/dionysus2/tutorial/cohomology.html#circular>`_.
  * (Optional) `Matplotlib <https://matplotlib.org/>`_ for plotting.

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

Documentation
-------------

Documentation for Dionysus can be found `here <http://mrzv.org/software/dionysus2/>`_.
