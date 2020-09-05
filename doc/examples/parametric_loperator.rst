.. _parametric_loperator:

Custom scalar types and ``parametric_loperator``
================================================

The goal of this example is twofold. First, it shows how to use user-defined
types as ``ScalarType``'s, i.e. types of coefficient in expressions.
Second, it demonstrates functionality of the ``parametric_loperator`` class.

As for the custom scalar type, we will define a class implementing the
algebra of polynomials

.. math::

  P_n(x) = C_0 + C_1 x + C_2 x^2 + \ldots + C_n x^n.

:math:`x` can be understood as a parameter (for instance, time or an external
field) a quantum-mechanical operator depends on. One could envision other
practically relevant choices of ``ScalarType`` such as interpolators based
on the Chebyshev polynomial expansion or Padé approximants. The expression
constructed in the program is the Hamiltonian of a displaced quantum harmonic
oscillator with frequency :math:`\omega_0 = 2` and the displacement depending on
an external parameter :math:`\lambda`, :math:`g(\lambda) = 3\lambda`.

.. math::

  \hat H = \frac{\omega_0}{2} +
           \omega_0 \left(a^\dagger - g(\lambda)\right)
           \left(a - g(\lambda)\right).

An instance of the ``parametric_loperator`` class is constructed out of this
expression and applied to the vacuum state for a few values of :math:`\lambda`.

.. math::

  \hat H|0\rangle = \omega_0\left(g^2(\lambda) + \frac{1}{2}\right)|0\rangle -
                    \omega_0 g(\lambda) |1\rangle.

The crucial point here is that the slow task of constructing
``parametric_loperator`` is performed only once. Substitution of :math:`\lambda`
values is delayed until the moment the operator acts on the state.

.. literalinclude:: ../../examples/parametric_loperator.cpp
  :language: cpp
  :lines: 22-
  :linenos:
