.. _usage:

How to use *libcommute* in your project
=======================================

There are two ways to use *libcommute*, depending on what build system your C++
project is based on. Either way, you would need a compiler that supports at
least C++11 and -- optionally -- C++17, if you are interested in using the
:ref:`dynamic index sequence feature <dyn_indices>`.

It is recommended to include ``<libcommute/libcommute.hpp>``, which pulls in
all *libcommute*'s headers. The entire library is under 5000 lines of code and
compilation times should normally not be a concern.

Makefiles/no build system
-------------------------

Assume that ``<libcommute>`` is either a directory with unpacked
*libcommute* sources or the :ref:`installation <installation>`
directory. Just one ``-I`` flag on compiler's command line is sufficient to
make *libcommute*'s header files visible to your code.

.. code:: shell

  g++ -O3 -I<libcommute>/include -o myprog myprog.cpp

(similar for ``clang++`` and other compilers).

pkg-config
----------

*libcommute* installs a pkg-config configuration file ``libcommute.pc`` under
``<libcommute_prefix>/share/pkgconfig``. One can run

.. code:: shell

  pkg-config --cflags libcommute

get the correct ``-I`` flag needed for compilation. This method requires
``<libcommute_prefix>/share/pkgconfig`` being part of pkg-config's lookup path
(modification of the ``PKG_CONFIG_PATH`` environment variable may be needed).

CMake
-----

Here is a minimal example of a root ``CMakeLists.txt`` script for your
project.

.. code:: cmake

  cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)

  project(myproject LANGUAGES CXX)

  # Change the C++ standard to '17' if you plan to use
  # the dynamic index sequence feature
  set(CMAKE_CXX_STANDARD 11)

  # libcommute_ROOT is the installation prefix of libcommute
  set(libcommute_DIR ${libcommute_ROOT}/lib/cmake)

  # Import libcommute target
  find_package(libcommute 0.7.2 CONFIG REQUIRED)

  # Build an executable called 'myprog'
  add_executable(myprog myprog.cpp)
  target_link_libraries(myprog PRIVATE libcommute::libcommute)

.. note::

  Starting from version 0.7.2, *libcommute* exports the namespaced CMake target
  ``libcommute::libcommute`` instead of just ``libcommute``.
