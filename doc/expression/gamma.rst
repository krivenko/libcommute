.. _generator_gamma:

Example: A user-defined algebra
===============================

.. default-domain:: cpp

Introducing a new :ref:`algebra <generator>` is as easy as deriving a class from
the abstract base :class:`libcommute::generator`. In this example we will define
the algebra of Dirac gamma matrices
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

.. literalinclude:: ../../examples/new_algebra.cpp
  :language: cpp
  :lines: 18-120

It is usually worth defining a factory function that creates an expression
containing one generator with a unity prefactor.

.. literalinclude:: ../../examples/new_algebra.cpp
  :language: cpp
  :lines: 121-132

Now we can check that generators of our algebra actually fulfil the canonical
anti-commutation relations.

.. literalinclude:: ../../examples/new_algebra.cpp
  :language: cpp
  :lines: 134-155

Let us also define the fifth gamma matrix

.. math::

  \gamma^5 = i \gamma^0 \gamma^1 \gamma^2 \gamma^3

and check that

.. math::

  (\gamma^5)^\dagger = \gamma^5,\\
  \{\gamma^5, \gamma^\mu\} = 0.

.. literalinclude:: ../../examples/new_algebra.cpp
  :language: cpp
  :lines: 157-
