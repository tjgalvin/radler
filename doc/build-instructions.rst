.. _buildinstructions:

Build instructions
==================

Radler can be installed as a stand-alone package, but is also installed as a part of `WSClean <https://wsclean.readthedocs.io>`_. 
If you only want to install WSClean, it is not necessary to build Radler yourself.

Installing from PyPI
~~~~~~~~~~~~~~~~~~~~
Radler can be installed from PyPI:

::

    pip install radler

This is the easiest way to install radler, but might not always contain the latest features. To be sure to have the latest features, build it from source.

Building from source
~~~~~~~~~~~~~~~~~~~~
Radler needs a number of dependencies in order to successfully compile. They can be installed with:

::

    apt install git make cmake libpython3-dev g++ casacore-dev \
    libboost-date-time-dev libcfitsio-dev libfftw3-dev libgsl-dev \
    libhdf5-dev pybind11-dev

Note that you need Ubuntu 25.04 or later. For older versions of Ubuntu, you need to build ``casacore`` from source first (see the ``docker`` directory for an example).

In order to be able to build the documentation with ``make doc`` and ``sphinx``, a few documentation tools need to be installed:

::

    apt -y install doxygen python3-pip
    pip3 install sphinx sphinx_rtd_theme breathe


Quick compilation guide
~~~~~~~~~~~~~~~~~~~~~~~

::

    git clone https://git.astron.nl/RD/Radler.git
    cd Radler
    mkdir build && cd build
    cmake -DBUILD_PYTHON_BINDINGS=On ..
    make
    make install


Installation options
~~~~~~~~~~~~~~~~~~~~

(Use :code:`ccmake` or :code:`cmake-gui` to configure all options interactively;
or use the :code:`-D` option to set the options on the command line.)

* :code:`BUILD_PYTHON_BINDINGS`: build Python module 'radler' to use Radler from Python
* :code:`BUILD_TESTING`: compile tests (requires Boost Unit Test Framework)

All other build options serve development purposes only, and can be left at the default values by a regular user.

All libraries are installed in :code:`<installpath>/lib`. The header files in
:code:`<installpath>/include`. The Python module in
:code:`<installpath>/lib/python{VERSION_MAJOR}.{VERSION_MINOR}/site-packages`. Depending on your configuration, it might be necessary to set
:code:`LD_LIBRARY_PATH` and :code:`PYTHONPATH` appropiately.
