.. _installation:

Installation instructions
=========================

*libcommute* is a header-only C++ library, which means it is usable without
installation. You can obtain the latest version of the library by making a local
clone of its `GitHub repository <https://github.com/krivenko/libcommute>`_.
Assuming that `git <https://git-scm.com/>`_ package is installed and visible via
the ``PATH`` environment variable, the following shell command will create a
directory ``libcommute.src`` containing the latest *libcommute* source files.

.. code-block:: shell

    $ git clone https://github.com/krivenko/libcommute.git libcommute.src

You can then use *libcommute* in your own code by passing
'``-I/path/to/libcommute.src/include``' to the compiler command line.

If your project uses `CMake <https://cmake.org/download/>`_ as its build system,
you might want to install *libcommute* header files along with CMake
configuration scripts. This will allow you to import *libcommute* targets,
as well as to run unit tests and build example programs.
The minimum required version of CMake is **3.10.0**.

The following sequence of shell commands will build unit tests and examples.

.. code-block:: shell

    $ mkdir libcommute.build && cd libcommute.build
    $ cmake ../libcommute.src                       \
      -DCMAKE_INSTALL_PREFIX=<libcommute_prefix>    \
      -DEXAMPLES=ON                                 \
      -DTESTS=ON
    $ make

Compilation of the tests can be disabled with CMake flag ``-DTESTS=OFF``
*(not recommended)*. Examples are compiled by default and can be disabled
with ``-DEXAMPLES=OFF``.

It is generally recommended to run the unit tests before installing the library
to detect potential problems early on.

.. code-block:: shell

    $ make test

Failed tests (having anything different from ``Passed`` in the test result
column) signal host-specific compilation/linking problems or bugs in the library
itself. The following command completes installation of library's files into
``<libcommute_prefix>``.

.. code-block:: shell

    $ make install

Documentation of *libcommute* can optionally be built and installed using the
``DOCUMENTATION`` CMake flag (requires
`Sphinx 2.1.0 or newer <https://www.sphinx-doc.org>`_,
`Read the Docs Sphinx theme <http://sphinx-rtd-theme.readthedocs.io/en/stable>`_
, `MathJax <https://www.mathjax.org/>`_ and
`MyST Parser <https://myst-parser.readthedocs.io/en/latest/>`_).

The table below gives a complete lists of supported CMake options with their
meaning.

+------------------------------+-----------------------------------------------+
| Option name                  | Description                                   |
+==============================+===============================================+
| ``CMAKE_INSTALL_PREFIX``     | Path to the directory *libcommute* will be    |
|                              | installed into.                               |
+------------------------------+-----------------------------------------------+
| ``CMAKE_BUILD_TYPE``         | CMake build type (``Release``, ``Debug`` or   |
|                              | ``RelWithDebInfo``) used to compile unit      |
|                              | tests and examples.                           |
+------------------------------+-----------------------------------------------+
| ``TESTS=[ON|OFF]``           | Enable/disable compilation of unit tests      |
|                              | (default: ``ON``).                            |
+------------------------------+-----------------------------------------------+
| ``EXAMPLES=[ON|OFF]``        | Enable/disable compilation of examples        |
|                              | (default: ``ON``).                            |
+------------------------------+-----------------------------------------------+
| ``USE_BOOST=[ON|OFF]``       | Try to build unit tests that use Boost        |
|                              | (default: ``ON``).                            |
+------------------------------+-----------------------------------------------+
| ``USE_GMPXX=[ON|OFF]``       | Try to build unit tests that use GMP C++      |
|                              | bindings (default: ``ON``).                   |
+------------------------------+-----------------------------------------------+
| ``USE_EIGEN3=[ON|OFF]``      | Try to build unit tests and examples that use |
|                              | Eigen 3 (default: ``ON``).                    |
+------------------------------+-----------------------------------------------+
| ``DOCUMENTATION=[ON|OFF]``   | Enable/disable generation of *libcommute*'s   |
|                              | Sphinx documentation (default: ``OFF``).      |
+------------------------------+-----------------------------------------------+
| ``Sphinx_ROOT``              | Path to Sphinx installation.                  |
+------------------------------+-----------------------------------------------+
| ``MathJax_ROOT``             | Path to MathJax installation (directory       |
|                              | containing ``MathJax.js``).                   |
+------------------------------+-----------------------------------------------+
| ``Eigen3_ROOT``              | Path to Eigen 3 installation (required to     |
|                              | build some optional unit tests).              |
+------------------------------+-----------------------------------------------+
| ``Boost_ROOT``               | Path to Boost installation (required to build |
|                              | some optional unit tests).                    |
+------------------------------+-----------------------------------------------+
| ``STATIC_ANALYSIS=[ON|OFF]`` | [DEV] Run static analysis tools on C++        |
|                              | sources of unit tests and examples.           |
+------------------------------+-----------------------------------------------+
| ``CLANG_TIDY_EXECUTABLE``    | [DEV] Path to the ``clang-tidy`` executable.  |
+------------------------------+-----------------------------------------------+
| ``CPPCHECK_EXECUTABLE``      | [DEV] Path to the ``cppcheck`` executable.    |
+------------------------------+-----------------------------------------------+
| ``CPPCHECK_EXTRA_FLAGS``     | [DEV] Semicolon-separated list of additional  |
|                              | flags to be passed to ``cppcheck``.           |
+------------------------------+-----------------------------------------------+
