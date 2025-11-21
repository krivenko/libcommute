.. _monomial_action:

Advanced: Linear operator representation of a user-defined algebra
==================================================================

.. default-domain:: cpp

.. namespace:: libcommute

Having :ref:`defined a new algebra <generator_gamma>` one can take a step
further and implement the :ref:`linear operator representation <loperator>`
of it. Let us build on the :ref:`example <generator_gamma>` where we introduced
the algebra of 4-dimensional :math:`\gamma`-matrices. This page will outline
the extra steps needed to plug the new algebra into the
:ref:`linear operator <loperator>` framework.

* Declare a new :ref:`elementary space <elementary_spaces>` class.
  The elementary space must share algebra ID with ``generator_gamma``.
  Its dimension is 4, so the :ref:`number of occupied bits <hilbert_space>` is
  :math:`b = 2`.

  .. literalinclude:: ../../examples/gamma.loperator.cpp
    :language: cpp
    :lines: 18-63

* Specialize class template :class:`monomial_action` for the new algebra. This
  specialization will describe how a product
  :math:`\gamma^\mu \gamma^\nu \gamma^\kappa \ldots` acts on a single basis
  state :math:`|n\rangle, n = \{0,1,2,3\}`.

  .. literalinclude:: ../../examples/gamma.loperator.cpp
    :language: cpp
    :lines: 65-144

  In general, :class:`monomial_action` is defined as follows.

  .. class:: template<int AlgebraID> monomial_action

    *Defined in <libcommute/loperator/monomial_action.hpp>*

    Action of a product of generators belonging to the same algebra
    :var:`AlgebraID` on a basis vector.

    .. function:: template<typename... IndexTypes> \
                  monomial_action( \
                    monomial<IndexTypes...>::range_type const& m_range, \
                    hilbert_space<IndexTypes...> const& hs)

      :var:`m_range` is a pair of :class:`monomial::const_iterator`'s.
      This range represents the product of generators one wants to act with.
      :var:`hs` is the Hilbert space to act in.

    .. function:: template<typename ScalarType> \
                  bool act(sv_index_type & index, ScalarType & coeff) const

      Act on a basis state. It is assumed the action of a generator (and,
      therefore, of a monomial) maps a single basis state to a single basis
      state multiplied by a constant.

      :var:`index` is the index of the input and output basis states.

      :var:`coeff` must be multiplied by the overall constant factor acquired
      as a result of monomial action.

      This method should return ``false`` if the action result is the identical
      zero and ``true`` otherwise.

* Make a :ref:`linear operator object <loperator>` with the new algebra ID and
  apply it to 4-dimensional state vectors as usual.

  .. literalinclude:: ../../examples/gamma.loperator.cpp
    :language: cpp
    :lines: 150-193
