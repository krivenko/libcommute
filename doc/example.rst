.. _examples:

Neat examples
=============

Spin-1/2 Heisenberg chain and its integrals of motion
-----------------------------------------------------

The spin-1/2 Heisenberg chain is a textbook example of an integrable quantum
system. Its Hamiltonian

.. math::

  \hat H = g \sum_i \mathbf{S}_i \cdot \mathbf{S}_{i+1}

conserves three projections of the total spin

.. math::

  \mathbf{S} = \sum_i \mathbf{S}_i

as well as a series of higher order charges :math:`Q_n`. Existence of these
charges can be derived from the transfer matrix theory. Explicit expressions
for :math:`Q_n` were obtained in [Grab1994]_. The following program constructs
Hamiltonian of the Heisenberg chain with periodic boundary conditions and
checks that :math:`[\hat H, \mathbf{S}] = 0`, :math:`[\hat H, Q_n] = 0` and
:math:`[Q_n, Q_m] = 0` for :math:`m,n = 3,4,5`.

.. literalinclude:: ../example/heisenberg_chain.cpp
  :language: cpp
  :lines: 25-
  :linenos:


Application of ``parametric_qoperator``
---------------------------------------

Virasoro algebra
----------------

TODO: polynomials of parameter as coefficients

.. [Grab1994] "Quantum Integrals of Motion for the Heisenberg Spin Chain",
   M. P. Grabowski and P. Mathieu,
   Mod. Phys. Lett. A, Vol. 09, No. 24, pp. 2197-2206 (1994),
   https://doi.org/10.1142/S0217732394002057
