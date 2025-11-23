.. _generator_gamma:

Advanced: A user-defined algebra
================================

.. default-domain:: cpp

Introducing a new :ref:`algebra <generator>` is as easy as deriving a class from
the abstract base :class:`libcommute::generator`. *libcommute*'s DSL can work
with algebras, whose generators :math:`g_\alpha` obey commutation relations

.. math::

  g_\alpha g_\beta - c g_\beta g_\alpha = F_{\alpha\beta} +
      \sum_\gamma f_{\alpha\beta}^\gamma g_\gamma

with real constants :math:`c`, :math:`F_{\alpha\beta}` and
:math:`f_{\alpha\beta}^\gamma`. Such algebraic structures include the Lie and
Clifford algebras.

In the following example we will define the algebra of Dirac gamma matrices
:math:`{\gamma^0,\gamma^1,\gamma^2,\gamma^3}`.
The generators will carry one integer index (0, 1, 2 or 3), which will also fix
their canonical ordering. The identities we will be using in the code below are

- Anti-commutation relation
  :math:`\{\gamma^\mu, \gamma^\nu\} = 2\eta^{\mu\nu} \hat I_4`, where
  :math:`\eta^{\mu\nu}` is the Minkowski metric tensor with signature
  :math:`(+,-,-,-)` and :math:`\hat I_4` is the four-dimensional identity
  matrix.
- Squares of gamma matrices,

  .. math::

    (\gamma^0)^2 &= \hat I_4,\\
    (\gamma^k)^2 &= -\hat I_4 \text{ for } k=1,2,3,\\

- Hermitian conjugates

  .. math::

    (\gamma^0)^\dagger &= \gamma^0,\\
    (\gamma^k)^\dagger &= -\gamma^k.

.. literalinclude:: ../../examples/gamma.hpp
  :language: cpp
  :lines: 18-116

.. note::

  If the constants :math:`c`, :math:`F_{\alpha\beta}`,
  :math:`f_{\alpha\beta}^\gamma` belong to a certain category of numbers
  (integers, rationals or general real values), then the generator is compatible
  with all :var:`ScalarType <libcommute::expression::ScalarType>`'s belonging to
  the same or wider category.

  :struct:`libcommute::linear_function` internally stores its coefficients as
  instances of the variant type :struct:`libcommute::var_number`. In the
  particular example presented here, all coefficients are integers, which are
  automatically converted to :struct:`libcommute::var_number` in the expressions
  similar to ``f.set(2 * eta)``. In the case of the rational coefficients,
  it is recommended to explicitly construct the values by calling
  ``var_number(numerator, denominator)`` instead of just passing floating-point
  numbers. This way, the defined generators can be used with a rational
  :var:`ScalarType <libcommute::expression::ScalarType>`.

It is usually worth defining a factory function that creates an expression
containing one generator with a unity prefactor.

.. literalinclude:: ../../examples/gamma.hpp
  :language: cpp
  :lines: 118-132

Now we can check that generators of our algebra actually fulfil the canonical
anti-commutation relations.

.. literalinclude:: ../../examples/gamma.expression.cpp
  :language: cpp
  :lines: 26-45

Let us also define the fifth gamma matrix

.. math::

  \gamma^5 = i \gamma^0 \gamma^1 \gamma^2 \gamma^3

and check that

.. math::

  (\gamma^5)^\dagger = \gamma^5,\\
  \{\gamma^5, \gamma^\mu\} = 0.

.. literalinclude:: ../../examples/gamma.expression.cpp
  :language: cpp
  :lines: 47-61
