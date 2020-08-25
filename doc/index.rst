libcommute
==========

*libcommute* is a C++11/14/17 template library that includes two major parts.

* A Domain-Specific Language (DSL) designed to easily construct and manipulate
  polynomial expressions with quantum-mechanical operators,
  especially those used in the quantum many-body theory. The most commonly
  used instances of such expressions are many-body Hamiltonians and operators
  of physical observables.

* A fast intermediate representation of the quantum-mechanical operators
  that enables their action on state vectors in finite-dimensional Hilbert
  spaces. This feature provides a basis for writing highly performant Exact
  Diagonalization (ED) codes without loss of flexibility.

*libcommute* offers a bunch of features and is extensible in a number of ways.

* The :cpp:class:`libcommute::expression` class is an arbitrary polynomial
  expression that can mix bosonic/fermionic
  (`CCR and CAR algebras <https://en.wikipedia.org/wiki/CCR_and_CAR_algebras>`_)
  and spin operators. It supports all standard arithmetic operations,
  compound-assignement operators, Hermitian conjugation and gives access to
  individual monomials via an iteration interface. No matrix representation is
  internally used to store the expressions, so there is virtually no limit on
  their size.

* Algebra generators such as fermionic creation/annihilation operators
  :math:`\hat c^\dagger_{i,j,k,\ldots}`
  (:math:`\hat c_{i,j,k,\ldots}`) are templated on the number and types of the
  indices :math:`i,j,k,\ldots`. It is, therefore, easy to implement notation
  that is most natural for the problem at hand.

* With a C++17-capable compiler, it is also possible to use dynamically-typed
  indices on the generators. More specifically, every single generator in an
  expression can carry its own sequence of indices :math:`i,j,k,\ldots`.
  Neither lengths of the sequences nor types of their elements have to agree
  between different generators.

* Coefficients of the polynomial expressions can be real, complex
  (``std::complex<T>``) or of any user-defined type ``S`` (one needs to
  specialize the structure ``scalar_traits<S>`` for that type).

* User-defined commutation/anticommutation operator algebras can be added via
  inheritance from the abstract base class :cpp:class:`libcommute::generator`.

* The common shorthand notation :math:`\pm \text{H.c.}` is supported in
  the expressions (credits to
  `Kristofer Björnson <https://github.com/dafer45>`_ for his original
  implementation of this idea in the
  `TBTK library <https://github.com/dafer45/TBTK>`_).

* The :cpp:class:`libcommute::qoperator` class translates an expression into
  a form, where it can be quickly acted with on a state vector.

* A state vector is a container type modeling a special
  :cpp:concept:`StateVector` concept. By default,
  :cpp:class:`libcommute::qoperator` will
  work with any container that gives read/write access to its elements
  via ``operator[]`` and exposes its size via a ``size()`` method
  (``std::vector<T>`` is one example).
  If there is a quicker way to perform a specific operation on a user-defined
  container type (for instance, setting all vector's elements to 0),
  it can be implemented by overloading the corresponding free function
  (part of the concept).

* `libcommute::qoperator` mirrors `libcommute::expression` in the manner
  it can be extended to support more commutation/anticommutation algebras.

* As a bonus, there is an utility class
  :cpp:class:`libcommute::space_partition`,
  which implements the Hilbert space partitioning algorithm
  (`Computer Physics Communications 200, March 2016, 274-284
  <http://dx.doi.org/10.1016/j.cpc.2015.10.023>`_, Sec. 4.2).

Follow the links in the contents table below to learn more about
*libcommute*'s capabilities and see some usage examples.

Being a header-only template library, *libcommute* does not require
installation. However, you may still choose to follow the :ref:`installation`
to build and run unit tests, build examples, and install the library so that it
is discoverable by other CMake-based projects. The source code repository and
the issue tracker are hosted on
`GitHub <https://github.com/krivenko/libcommute>`_.

*libcommute* is distributed under the terms of the *Mozilla Public License,
v. 2.0*. You can obtain a copy of the MPL at http://mozilla.org/MPL/2.0/.

Contents
========

.. toctree::
    :name: mastertoc
    :maxdepth: 3

    installation
    usage
    expression/index
    qoperator/index
    examples/index
    genindex
    search

