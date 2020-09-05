.. _expr_intro:

Introduction
============

Dealing with involved Hamiltonians and other operators from the quantum
many-body theory in C++ can be tough. In some cases, quantum-mechanical
operators can be uniformly represented by finite-dimensional matrices.
This representation, however, is often impractical. Storing the matrices can
quickly prove infeasible as the amount of required memory grows exponentially
with the number of degrees of freedom. For this reason, many computational
programs in the field use sparse matrices or hard-coded procedures that describe
how said operators act on quantum states. These implementations usually accept
a few Hamiltonian parameters as input, but switching to a more general
form/adding more terms to the Hamiltonian requires a considerable code rewrite.

The goal of *libcommute*'s Domain-Specific Language (DSL) is to streamline this
coding task. It introduces an abstraction of the polynomial quantum-mechanical
operator expression, which can be manipulated as easily as an equation written
on a piece of paper.

As a primer, let us consider the following simple program that constructs
Hamiltonian of an electronic tight-binding model on a square :math:`10\times 10`
lattice with only nearest-neighbour hopping allowed,

.. math::

  \hat H_\text{e} = -t \sum_\sigma \sum_{\langle i,j\rangle}
                    c^\dagger_{i,\sigma} c_{j,\sigma}.

.. literalinclude:: ../../examples/holstein.cpp
  :language: cpp
  :lines: 18-68

Now, let us add a harmonic oscillator at each lattice site (a localized phonon),

.. math::

  \hat H_\text{ph} = \omega_0 \sum_i a^\dagger_i a_i.

.. literalinclude:: ../../examples/holstein.cpp
  :language: cpp
  :lines: 71-91

.. note::

  We had to assign an empty spin label "" to the bosons, because all operators
  in ``H_ph`` have to carry exactly three indices with the last one being a
  string. It is possible to overcome this limitation and put just two integer
  indices on :math:`a^\dagger`/:math:`a` by switching to the
  :ref:`dynamically-typed indices <dyn_indices>`.
  Be aware, however, that the dynamic indices require C++17 and may result
  in less type-safe code.

Finally, we are going to couple electrons with phonons and arrive at
Hamiltonian of the Holstein model :math:`\hat H_H`,

.. math::

  \hat H_\text{e-ph} = g \sum_\sigma \sum_i n_{i,\sigma}(a^\dagger_i + a_i),

.. math::

  \hat H_H = \hat H_\text{e} + \hat H_\text{ph} + \hat H_\text{e-ph}.

.. literalinclude:: ../../examples/holstein.cpp
  :language: cpp
  :lines: 93-

Now that we have the complete Hamiltonian object, we could proceed along one
of the following routes.

- :ref:`Make a <ed_intro>` ``loperator`` object out of :math:`\hat H_H` and
  use it to act on state vectors in a finite-dimensional Hilbert space.
  This is a common step in implementing Exact Diagonalization algorithms.
- Use the :ref:`iteration interface <expr_iteration>` to analyze the structure
  of the Hamiltonian term by term.
- :ref:`Transform <expr_iteration>` the Hamiltonian by applying a function
  to each term.
