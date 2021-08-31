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
    - Type of the elements.
    - :expr:`T`

  * - :expr:`get_element(sv, n)`
    - Return the :expr:`n`-th element of :expr:`sv`.
    - :expr:`sv[n]`

  * - :expr:`update_add_element(sv, n, value)`
    - Add a value of some type :type:`U` to the :expr:`n`-th element of
      :expr:`sv`.
    - :expr:`sv[n] += value` or :expr:`sv[n] = sv[n] + value`

      The compound-assignment from type :type:`U` will be used
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
index arguments of :func:`get_element()` and :func:`update_add_element()`
according to a predefined map :type:`sv_index_type` -> :type:`sv_index_type`.
The element access functions throw :type:`std::out_of_range` if their index
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
:func:`basis_mapper:: make_view()`/:func:`basis_mapper::make_const_view()`.

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
                HSType const& hs, unsigned int N)

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
                make_view(StateVector && sv) const
                template<typename StateVector> \
                mapped_basis_view<StateVector const> \
                make_const_view(StateVector && sv) const

    Make a read/write or constant view of :expr:`sv`.
    Constant views will not be accepted by :func:`update_add_element()`.
    If :expr:`sv` is not an lvalue reference, the resulting view will
    :ref:`hold a copy <mapped_basis_view_Ref>` of :expr:`sv`.

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

.. _n_fermion_sector_view:

N-fermion sector views
----------------------

There are two more specialized flavours of the basis mapping views called
:math:`N`-fermion sector views and :math:`N`-fermion multisector views. They
can come in handy when working with particle-number preserving models of
fermions. If the model is large, then generating and storing a basis state index
map for :type:`mapped_basis_view` may become too expensive.

.. class:: template<typename StateVector, bool Ref = true> n_fermion_sector_view

  *Defined in <libcommute/loperator/n_fermion_sector_view.hpp>*

  View of a :type:`StateVector` object that translates basis state indices from
  a full :ref:`Hilbert space <hilbert_space>` to its subspace (sector) with a
  fixed total occupation of fermionic degrees of freedom :math:`N`. The full
  Hilbert space does not have to be purely fermionic.

  :type:`n_fermion_sector_view` is generally less performant than
  :type:`mapped_basis_view` in terms of the index translation speed. However,
  its required storage space scales only as :math:`O(M \min(N, M - N))`, where
  :math:`M` is the total number of the fermionic degrees of freedom. This
  scaling law is much milder that the exponential growth of the sector size.

  :type:`StateVector` - type of the underlying state vector object. Defining a
  read-only view (such that prohibits :expr:`update_add_element()` operations)
  requires using a ``const``-qualified type here. For example, one can use
  ``StateVector = std::vector<double>`` for a read-write view, and
  ``StateVector = const std::vector<double>`` for a read-only view.

  .. _n_fermion_sector_view_Ref:

  :type:`Ref` - by default, :type:`n_fermion_sector_view`
  stores a reference to the underlying state vector. Setting this option to
  ``false`` will result in a copy being created and stored instead. This feature
  can be useful when the underlying type is already a view-like object similar
  to ``Eigen::Map``.

  .. function:: template <typename SV, typename HSType> \
                         n_fermion_sector_view(SV&& sv, \
                         HSType const& hs, unsigned int N)

    Construct a view of the state vector :expr:`sv`, defined in the
    :expr:`N`-fermion sector of the full Hilbert space :expr:`hs`.

  .. function:: sv_index_type map_index(sv_index_type index) const

    Translate a basis state :expr:`index` from the full Hilbert space to the
    sector.

.. struct:: template <typename HSType> sector_descriptor

  Description of an :math:`N`-fermion sector defined over a subset of fermionic
  degrees of freedom.

  :type:`HSType` - type of the full Hilbert space this sector belongs to.

  .. member:: std::set<typename HSType::index_types> indices

    Set of indices corresponding to the relevant fermionic degrees of freedom.

  .. member:: unsigned int N

    Total occupation of the sector.

.. class:: template<typename StateVector, bool Ref = true> \
           n_fermion_multisector_view

  *Defined in <libcommute/loperator/n_fermion_sector_view.hpp>*

  View of a :type:`StateVector` object that translates basis state indices from
  a full :ref:`Hilbert space <hilbert_space>` to an :math:`N`-fermion
  multisector. A multisector is a set of all basis states, which have
  :math:`N_1` particles within a subset of fermionic modes :math:`\{S_1\}`,
  :math:`N_2` particles within another subset :math:`\{S_2\}` and so on. There
  can be any number of individual pairs :math:`(\{S_i\}, N_i)` (sectors
  contributing to the multisector) as long as all subsets :math:`\{S_i\}` are
  disjoint. The full Hilbert space does not have to be purely fermionic.

  :type:`n_fermion_multisector_view` is generally less performant than
  :type:`mapped_basis_view` in terms of the index translation speed. However,
  its required storage space scales only as
  :math:`O(\sum_i M_i \min(N_i, M_i - N_i))`, where
  :math:`M_i = |\{S_i\}|`. This scaling law is much milder that the exponential
  growth of the multisector size.

  It is advised to use :type:`n_fermion_sector_view` instead, if there is only
  one contributing sector that also spans all fermionic degrees of freedom.

  :type:`StateVector` - type of the underlying state vector object. Defining a
  read-only view (such that prohibits :expr:`update_add_element()` operations)
  requires using a ``const``-qualified type here. For example, one can use
  ``StateVector = std::vector<double>`` for a read-write view, and
  ``StateVector = const std::vector<double>`` for a read-only view.

  .. _n_fermion_multisector_view_Ref:

  :type:`Ref` - by default, :type:`n_fermion_multisector_view`
  stores a reference to the underlying state vector. Setting this option to
  ``false`` will result in a copy being created and stored instead. This feature
  can be useful when the underlying type is already a view-like object similar
  to ``Eigen::Map``.

  .. function:: template <typename SV, typename HSType> \
                n_fermion_multisector_view(SV&& sv, HSType const& hs, \
                std::vector<sector_descriptor<HSType>> const& sectors)

    Construct a view of the state vector :expr:`sv`, defined in the
    :math:`N`-fermion multisector of the full Hilbert space :expr:`hs`.
    The multisector is defined via a list of contributing :expr:`sectors`
    (list of :math:`(\{S_i\}, N_i)` pairs).

  .. function:: sv_index_type map_index(sv_index_type index) const

    Translate a basis state :expr:`index` from the full Hilbert space to the
    multisector.

Besides :class:`n_fermion_sector_view` and :class:`n_fermion_multisector_view`,
*<libcommute/loperator/n_fermion_sector_view.hpp>* defines a few supplemental
utility functions that help working with (multi)sectors.

.. function:: template <typename StateVector, typename HSType> \
              auto make_nfs_view(StateVector&& sv, HSType const& hs, \
              unsigned int N)
              template <typename StateVector, typename HSType> \
              auto make_const_nfs_view(StateVector&& sv, HSType const& hs, \
              unsigned int N)

  Make and return a read/write or constant :expr:`N`-fermion sector view of
  :expr:`sv` within the full Hilbert space :expr:`hs`. If :expr:`sv` is not an
  lvalue reference, the resulting view will
  :ref:`hold a copy <n_fermion_sector_view_Ref>` of :expr:`sv`.

.. function:: template <typename StateVector, typename HSType> \
              auto make_nfms_view(StateVector&& sv, HSType const& hs, \
              std::vector<sector_descriptor<HSType>> const& sectors)
              template <typename StateVector, typename HSType> \
              auto make_const_nfms_view(StateVector&& sv, HSType const& hs, \
              std::vector<sector_descriptor<HSType>> const& sectors)

  Make and return a read/write or constant :math:`N`-fermion multisector view of
  :expr:`sv` within the full Hilbert space :expr:`hs`. The multisector is
  defined via a list of contributing :expr:`sectors` (list of
  :math:`(\{S_i\}, N_i)` pairs). If :expr:`sv` is not an lvalue reference,
  the resulting view will
  :ref:`hold a copy <n_fermion_sector_view_Ref>` of :expr:`sv`.

.. function:: template <typename HSType> sv_index_type \
              n_fermion_sector_size(HSType const& hs, unsigned int N)

  Size of the :expr:`N`-fermion sector within the full Hilbert space :expr:`hs`.

.. function:: template <typename HSType> sv_index_type \
              n_fermion_multisector_size(HSType const& hs, \
              std::vector<sector_descriptor<HSType>> const& sectors)

  Size of the :math:`N`-fermion multisector within the full Hilbert space
  :expr:`hs`. The multisector is defined via a list of contributing
  :expr:`sectors` (list of :math:`(\{S_i\}, N_i)` pairs).

.. function:: template <typename HSType> std::vector<sv_index_type> \
              n_fermion_sector_basis_states(HSType const& hs, unsigned int N)

  Build and return a list of basis state indices forming the :expr:`N`-fermion
  sector within the full Hilbert space :expr:`hs`. The order of the indices in
  the list is consistent with the results of
  :func:`n_fermion_sector_view::map_index()`.

  .. code-block:: cpp

    auto basis_states = n_fermion_sector_basis_states(hs, N);
    auto view = n_fermion_sector_view(st, hs, N);

    for(sv_index_type n = 0; n < basis_states.size(); ++n) {
      view.map_index(basis_states[n]) == n; // true for all n
    }

.. function:: template <typename HSType> std::vector<sv_index_type> \
              n_fermion_multisector_basis_states(HSType const& hs, \
              std::vector<sector_descriptor<HSType>> const& sectors)

  Build and return a list of basis state indices forming an :math:`N`-fermion
  multisector within the full Hilbert space :expr:`hs`. The multisector is
  defined via a list of contributing :expr:`sectors` (list of
  :math:`(\{S_i\}, N_i)` pairs). The order of the indices in the list is
  consistent with the results of
  :func:`n_fermion_multisector_view::map_index()`.

  .. code-block:: cpp

    auto basis_states = n_fermion_multisector_basis_states(hs, sectors);
    auto view = n_fermion_multisector_view(st, hs, sectors);

    for(sv_index_type n = 0; n < basis_states.size(); ++n) {
      view.map_index(basis_states[n]) == n; // true for all n
    }
