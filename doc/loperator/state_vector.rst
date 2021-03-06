.. _state_vectors:

State vectors
=============

.. _state_vector:

``StateVector`` concept
-----------------------

.. default-domain:: cpp

.. namespace:: libcommute

Let us say we want to make type :type:`SV` a *libcommute*-compatible state
vector type so that :ref:`linear operators <loperator>` can act on instances
of :type:`SV`. For this we have to make :type:`SV` model the ``StateVector``
concept.
In a nutshell, a ``StateVector`` type is a one-dimensional array of numbers
(quantum amplitudes) allowing integer indexing and implementing a certain
interface. Elements of the array do not have to be stored contiguously.
Acceptable index values must be at least 64-bit wide unsigned integers since
*libcommute* uses the following type for basis state indexing.

.. type:: sv_index_type = std::uint64_t

  *Defined in <libcommute/loperator/state_vector.hpp>*

  Type of basis state indices.

The table below shows the interface (a set of free functions and a
metafunction) that needs be implemented for an object :expr:`sv` of type
:type:`SV`.

*libcommute* provides an implementation of the ``StateVector`` concept for
:type:`std::vector` (see *<libcommute/loperator/state_vector.hpp>*).

.. list-table::
  :header-rows: 1

  * - Function/metafunction
    - Description
    - Implementation for :expr:`std::vector<T>`
  * - :expr:`element_type<SV>::type`
    - The type of the elements.
    - :expr:`T`

  * - :expr:`get_element(sv, n)`
    - Return the :expr:`n`-th element of :expr:`sv`.
    - :expr:`sv[n]`

  * - :expr:`update_add_element(sv, n, value)`
    - Add a value of some type :type:`T` to the :expr:`n`-th element of
      :expr:`sv`.
    - :expr:`sv[n] += value` or :expr:`sv[n] = sv[n] + value`

      The compound-assignment from type :expr:`T` will be used
      whenever :expr:`sv`'s elements support it. Otherwise, the implementation
      will fall back to the simple addition.

  * - :expr:`set_zeros(sv)`
    - Fill :expr:`sv` with zeros.
    - :expr:`std::fill(sv.begin(), sv.end(), zero)`.

      The zero value is created by
      :expr:`make_const(0)` as described in ":ref:`custom_scalar_type`".

  * - :expr:`zeros_like(sv)`
    - Return an object of the same type and size as :expr:`sv` but filled with
      zeros.
    - Creates a new object as :expr:`std::vector<T>(sv.size(), zero)`.

  * - :expr:`foreach(sv, f)`
    - Apply a function-like object :expr:`f` to all basis state index/non-zero
      element pairs :expr:`(n, a)` in :expr:`sv`.
    - In a for-loop, calls :expr:`f(n, a)` for all non-zero elements :expr:`a`
      as detected by :expr:`is_zero()` (see ":ref:`custom_scalar_type`").

Inclusion of *<libcommute/loperator/state_vector_eigen3.hpp>* makes some
`Eigen 3 <https://eigen.tuxfamily.org/>`_ types (`column vectors`_,
`vector blocks`_,
`column-like matrix blocks`_ and one-dimensional `Eigen::Map views`_)
compatible with the ``StateVector`` concept as well.

.. _column vectors:
  https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
  #TutorialMatrixVectors

.. _vector blocks:
  https://eigen.tuxfamily.org/dox/classEigen_1_1VectorBlock.html

.. _column-like matrix blocks:
  https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
  #TutorialBlockOperationsSyntaxColumnRows

.. _Eigen::Map views:
  https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html


.. _sparse_state_vector:

Sparse state vector
-------------------

:class:`sparse_state_vector` is a state vector that saves memory by storing only
the non-zero elements. It is essentially a wrapper around
:class:`std::unordered_map` modelling the ``StateVector`` concept. Here, we show
only the part of its interface not covered by ``StateVector``.

.. class:: template<typename ScalarType> sparse_state_vector

  State vector with a sparse storage of elements (quantum amplitudes).
  :expr:`ScalarType` is the type of the elements.

  .. function::   sparse_state_vector() = delete
                  sparse_state_vector(sv_index_type size)

    Construct a zero (empty) sparse vector with a given :expr:`size` --
    dimension of the corresponding Hilbert space.

  .. function:: sv_index_type size() const

    Size (dimension) of the vector.

  .. function:: ScalarType & operator[](sv_index_type n)

    Access the :expr:`n`-th element. If it is zero (missing from the storage),
    then a new value-initialized element will be inserted and a reference to
    it will be returned.

    .. warning::

      Improper use of this method may result in zero elements being stored in
      the unordered map. Only the non-zero values should be assigned to the
      references returned by it.

  .. function:: sv_index_type n_nonzeros() const

    Get the number of non-zero (stored) elements.

.. _mapped_basis_view:

Mapped basis view
-----------------

:class:`mapped_basis_view` is another utility type modelling the ``StateVector``
concept. It is a view of a state vector, which translates basis state
index arguments of :expr:`get_element()` and :expr:`update_add_element()`
according to a predefined map :expr:`sv_index_type` -> :expr:`sv_index_type`.
The element access functions throw :expr:`std::out_of_range` if their index
argument is missing from the map.

:class:`mapped_basis_view` can be used in situations
where a :ref:`linear operator <loperator>` acts in a small subspace of
a full Hilbert space, and it is desirable to store vector components only within
that subspace. Such a situation naturally emerges when working with
:ref:`invariant subspaces of operators <space_partition>`.

.. class:: template<typename StateVector, bool Ref = true> mapped_basis_view

  View of a :type:`StateVector` object that translates basis state indices
  according to a certain mapping.

  :type:`StateVector` - type of the underlying state vector object. Defining a
  read-only view (such that prohibits :expr:`update_add_element()` operations)
  requires using a ``const``-qualified type here. For example, one can use
  ``StateVector = std::vector<double>`` for a read-write view, and
  ``StateVector = const std::vector<double>`` for a read-only view.

  .. _mapped_basis_view_Ref:

  :type:`Ref` - by default, :type:`mapped_basis_view`
  stores a reference to the underlying state vector. Setting this option to
  ``false`` will result in a copy being created and stored instead. This feature
  can be useful when the underlying type is already a view-like object similar
  to ``Eigen::Map``.

The mapped basis views should always be constructed by means of a special
factory class :class:`basis_mapper` and its methods
:func:`make_view()`/:func:`make_const_view()`/:func:`make_view_no_ref()`/
:func:`make_const_view_no_ref`.

.. class:: basis_mapper

  Factory class for :class:`mapped_basis_view`.

  .. rubric:: Constructors

  .. function:: basis_mapper(std::vector<sv_index_type> const& \
                             basis_state_indices)

    Build a mapping from a list of basis states :expr:`basis_state_indices`
    to their positions within the list.

    .. code-block:: cpp

      std::vector<sv_index_type> basis_indices{3, 5, 6};
      basis_mapper mapper(basis_indices);

      // Views created by 'mapper' will translate basis state indices
      // according to
      // 0 -> std::out_of_range
      // 1 -> std::out_of_range
      // 2 -> std::out_of_range
      // 3 -> 0
      // 4 -> std::out_of_range
      // 5 -> 1
      // 6 -> 2
      // 7 -> std::out_of_range
      // ...

  .. function:: template<typename HSType, \
                         typename LOpScalarType, \
                         int... LOpAlgebraIDs> \
                basis_mapper(loperator<LOpScalarType,LOpAlgebraIDs...>const& O,\
                             HSType const& hs)

    Build a mapping from a set of all basis states contributing to
    :math:`\hat O|0\rangle`.

    Operator :expr:`O` acts in the Hilbert space :expr:`hs`.
    :math:`|0\rangle` is the basis state with index 0 ('vacuum' state in
    the case of fermions and bosons).
    Mapped values are assigned continuously starting from 0 without any specific
    order.

  .. function:: template<typename HSType, \
                         typename LOpScalarType, \
                         int... LOpAlgebraIDs> \
                basis_mapper( \
                std::vector<loperator<LOpScalarType, LOpAlgebraIDs...>> \
                  const& O_list, \
                HSType const& hs, int N)

    Given a list of operators
    :math:`\{\hat O_1, \hat O_2, \hat O_3, \ldots, \hat O_M\}`, build a mapping
    from all basis states contributing to all states
    :math:`\hat O_1^{n_1} \hat O_2^{n_2} \ldots \hat O_M^{n_M} |0\rangle`,
    where :math:`n_m \geq 0` and :math:`\sum_{m=1}^M n_M = N`.

    Operators in :expr:`O_list` act in the Hilbert space :expr:`hs`.
    :math:`|0\rangle` is the basis state with index 0 ('vacuum' state in
    the case of fermions and bosons).
    Mapped values are assigned continuously starting from 0 without any specific
    order.

    This constructor is useful to create a mapping from a fixed-particle-number
    subspace of a fermionic/bosonic Hilbert space.

  .. rubric:: :class:`mapped_basis_view` factory functions

  .. function:: template<typename StateVector> \
                mapped_basis_view<StateVector> \
                make_view(StateVector & sv) const
                template<typename StateVector> \
                mapped_basis_view<const StateVector> \
                make_const_view(StateVector const& sv) const

    Make a read/write or constant view of :expr:`sv`.
    Constant views will not be accepted by :expr:`update_add_element()`.

  .. function:: template<typename StateVector> \
                mapped_basis_view<StateVector, false> \
                make_view_no_ref(StateVector sv) const
                template<typename StateVector> \
                mapped_basis_view<const StateVector, false> \
                make_const_view_no_ref(StateVector sv) const

    Make a read/write or constant view
    :ref:`holding a copy <mapped_basis_view_Ref>` of :expr:`sv`. Constant views
    will not be accepted by :expr:`update_add_element()`.

    .. warning::

      To reduce memory footprint, :class:`mapped_basis_view` objects store
      a reference to the basis index map owned by their parent
      :class:`basis_mapper` object. For this reason, the views should never
      outlive the mapper.

  .. rubric:: Other methods

  .. function:: sv_index_type size() const

    Number of elements in the index map.

  .. function:: std::unordered_map<sv_index_type, sv_index_type> \
                const& map() const

    Direct access to the underlying index map.

  .. function:: std::unordered_map<sv_index_type, sv_index_type> \
                inverse_map() const

    Build and return an inverse index map. Depending on map's size, building
    the inverse can be an expensive operation. Calling this method on a
    non-invertible map is undefined behavior.
