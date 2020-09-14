.. _ed_intro:

Introduction
============

Sometimes, algebraic manipulations with quantum-mechanical operators can be an
important part of a numerical algorithm in their own right. It is, however, much
more common to use the resulting algebraic expressions as operators acting in a
space of quantum states (Hilbert space). This scenario is of particular
importance when implementing algorithms for partial or full diagonalization
of quantum Hamiltonians -- so called exact diagonalization (ED) methods.

*libcommute* gives user a handful of tools to minimize overhead when
writing ED codes. The following few steps describe a typical workflow of a
*libcommute*-based ED code.

* Define Hamiltonian of a quantum system as a
  :ref:`polynomial expression <expression>`.

* Make a finite-dimensional :ref:`Hilbert space <hilbert_space>`. This can
  be done either automatically for a given Hamiltonian, or -- in more tricky
  cases -- by explicitly constructing a product of
  :ref:`elementary spaces <elementary_spaces>`.

* Use the polynomial expression and the Hilbert space to make a
  :ref:`linear operator <loperator>` object representing the Hamiltonian.
  This object can be either a :ref:`simple linear operator <simple_loperator>`
  or a :ref:`parametric linear operator <param_loperator>`.

* Optionally :ref:`find invariant subspaces <space_partition>` of the
  Hamiltonian to reduce the numerical cost of the diagonalization procedure.
  The :cpp:class:`mapped_basis_view` state vector adapter can be acted upon by
  linear operators defined in the full Hilbert space, while saving memory by
  storing only the part of a state vector belonging to a certain invariant
  subspace.

* Act with the constructed linear operator on
  :ref:`state vectors <state_vector>` to obtain a matrix representation of the
  Hamiltonian. Instead of computing all matrix elements at once and resorting
  to a dense diagonalization method such as QR algorithm, one could opt to
  power iteration or Lanczos algorithm for sparse eigenproblems. In the latter
  case, one must act with the operator on a state at each step of the iterative
  diagonalization scheme.
