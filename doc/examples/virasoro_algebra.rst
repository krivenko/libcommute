.. _virasoro_algebra:

Virasoro algebra
================

Virasoro algebra is the central extension of a Lie algebra with wide applications
in the two-dimensional conformal field theory. Its generators
:math:`L_n, n \in \mathbb{Z}` satisfy commutation relations

.. math::

  [L_m, L_n] = (m - n) L_{m+n} + c(m^3 - m) \delta_{m, -n},

where :math:`c` is the central charge commuting with all generators. It was
shown in [Fair1988]_ that the Virasoro algebra can be constructed out of just
two generators :math:`L_3` and :math:`L_{-2}` using the following recurrence
relations,

.. math::

  \begin{align}
    L_1 &= \frac{1}{5}[L_3, L_{-2}],\\
    L_{-1} &= \frac{1}{3}[L_1, L_{-2}],\\
    L_2 &= \frac{1}{4}[L_3, L_{-1}],\\
    L_0 &= \frac{1}{2}[L_1, L_{-1}],\\
    L_{n+1} &= \frac{1}{n-1}[L_n, L_1] \text{ for } n>2,\\
    L_{-n-1} &= \frac{1}{1-n}[L_{-n}, L_{-1}] \text{ for } n>1.
  \end{align}

In the example below, we show how to implement the Virasoro algebra in
*libcommute*'s framework and use it to verify the recurrence relations stated
above.

.. literalinclude:: ../../examples/virasoro_algebra.cpp
  :language: cpp
  :lines: 25-
  :linenos:

.. [Fair1988] "A presentation for the Virasoro and super-Virasoro algebras",
   D. B. Fairlie, J. Nuyts and C. K. Zachos ,
   Commun. Math. Phys. **117**, pp. 595–614 (1988),
   https://doi.org/10.1007/BF01218387
