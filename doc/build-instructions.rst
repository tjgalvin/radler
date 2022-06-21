.. _buildinstructions:

Build instructions
==================

Radler can be installed as a stand-alone package, but is also installed as a part of `WSClean <https://wsclean.readthedocs.io>`_. 
If you only want to install WSClean, it is not necessary to build Radler yourself.

Dependencies
~~~~~~~~~~~~
Radler needs a number of dependencies in order to successfully compile. On a clean (ubuntu 22.04) system,
the dependencies can be installed with (see also the ``docker`` directory):

General packages:

::

    apt-get -y install wget git make cmake g++ doxygen \
    libboost-all-dev libhdf5-dev libfftw3-dev \
    libblas-dev liblapack-dev libgsl-dev libxml2-dev \
    libgtkmm-3.0-dev libpython3-dev python3-distutils

Astronomy-specific packages:

::

    apt-get -y install casacore-dev libcfitsio-dev

In order to be able to build the documentation with ``make doc``, ``sphinx`` and some other documentation tools need to be installed:

::

    pip3 install sphinx doxygen sphinx_rtd_theme breathe myst-parser




Quick installation guide
~~~~~~~~~~~~~~~~~~~~~~~~

::

    git clone --recursive https://git.astron.nl/RD/Radler.git
    cd Radler
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=<radler_install_path> ..
    make install


Installation options
~~~~~~~~~~~~~~~~~~~~

(Use :code:`ccmake` or :code:`cmake -i` to configure all options.)

* :code:`BUILD_WITH_PYTHON`: build Python module 'radler' to use Radler from Python
* :code:`BUILD_TESTING`: compile tests

All other build options serve development purposes only, and can/should be left at the default values by a regular user.

All libraries are installed in :code:`<installpath>/lib`. The header files in
:code:`<installpath>/include`. The Python module in
:code:`<installpath>/lib/python{VERSION_MAJOR}.{VERSION_MINOR}/site-packages`. Make sure that your
:code:`LD_LIBRARY_PATH` and :code:`PYTHONPATH` are set as appropiate.
