.. _loperator:

Linear operators
================

.. default-domain:: cpp
.. namespace:: libcommute

*libcommute*'s linear operators implement action of
:ref:`polynomial expressions <expression>` on
:ref:`state vectors <state_vector>` in
:ref:`finite-dimensional Hilbert spaces <hilbert_space>`. Linear operator
classes are core devices used to build exact diagonalization routines based on
*libcommute*.

The linear operators come in two flavors, simple (:class:`loperator`)
and parametric (:class:`parametric_loperator`).

.. _simple_loperator:

Simple linear operator
----------------------

Simple linear operators are instances of class template :class:`loperator`.

Let us assume that :expr:`expr` is a real expression involving only the
predefined algebras (fermions, bosons and spins), and :expr:`hs` is a Hilbert
space object compatible with :expr:`expr`. Then creating a simple linear
operator is as easy as

.. code-block:: cpp

  auto Q = libcommute::make_loperator(expr, hs);

Futhermore, let us declare a couple of objects that model the
:ref:`state vectors <state_vector>` concept (in this particular case, they are
standard vectors)

.. code-block:: cpp

  std::vector<double> psi(hs.dim());
  std::vector<double> phi(hs.dim());

The following code fragment shows three (nearly) equivalent ways to act with our
operator on a state, :math:`|\phi\rangle = \hat Q |\psi\rangle`.

.. code-block:: cpp

  // Function call syntax
  phi = Q(psi);

  // Multiplication syntax
  phi = Q * psi;

  // 'In-place' action syntax
  Q(psi, phi);

Although all three forms are semantically equivalent, the last one is faster as
it eliminates the need for a temporary object to store
:math:`\hat Q|\psi\rangle`.

.. note::

  An important note must be made about the template parameters of
  :class:`loperator`. Expressions can dynamically accommodate new algebras via
  inheritance from the polymorphic base :class:`generator`. Unlike them, linear
  operators represent action of a certain set of algebras that is fixed at
  compile time. This design decision allows to remove the virtual function call
  overhead from the hottest parts of codes, where linear operators are
  repeatedly applied to state vectors. Template parameter pack
  :expr:`loperator::AlgebraIDs` (sorted list of integers) determines what
  algebras a given linear operator object supports. :func:`make_loperator`
  returns operators covering all three predefined algebras. For the sake of
  optimization, it may make sense to manually create :class:`loperator` objects
  with a more restricted set of IDs if some of the predefined algebras will
  never be used in expressions.

.. class:: template<typename ScalarType, int... AlgebraIDs> loperator

.. function:: template<typename ScalarType, typename... IndexTypes> \
              loperator<ScalarType, fermion, boson, spin> \
              make_loperator(expression<ScalarType, IndexTypes...> const& expr,\
              hilbert_space<IndexTypes...> const& hs)


.. _param_loperator:

Parametric linear operator
--------------------------

.. class:: template<typename ScalarType, int... AlgebraIDs> \
           parametric_loperator

.. function:: template<typename ScalarType, typename... IndexTypes> \
              parametric_loperator<ScalarType, fermion, boson, spin> \
              make_param_loperator( \
                expression<ScalarType, IndexTypes...> const& expr, \
                hilbert_space<IndexTypes...> const& hs)
