Python interface
================

Radler provides an interface via the :code:`radler` module. In order to build
this module, make sure the cmake configuration option :code:`BUILD_PYTHON_BINDINGS` is turned :code:`ON`.

After successfully compiling and installing the python bindings, you should update your :code:`LD_LIBRARY_PATH`
and your :code:`PYTHONPATH` as follows:

::

    export LD_LIBRARY_PATH=<installpath>/lib/:$LD_LIBRARY_PATH
    export PYTHONPATH=<installpath>/lib/python<VERSION_MAJOR>.<VERSION_MINOR>/site-packages:$PYTHONPATH

where :code:`VERSION_MAJOR` and :code:`VERSION_MINOR` depend upon the specific version of python on your system.
The :code:`radler` module can now be imported in python with:

::

    import radler

.. toctree::
   :maxdepth: 1

   python/settings
