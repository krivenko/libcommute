# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - Unreleased

- Better support for user-defined scalar types in ``expression``. In particular,
  it is now possible to use integer-like and rational-like numeric types, as
  requested by Dr. Cezary Śliwa. Possible incompatibility of a scalar type with
  a certain algebra generator (e.g. spin operators vs an integer scalar) is
  checked at runtime. This is achieved by storing the structure constants of the
  algebra in the tagged union type ``var_number``, which retains information
  about the category of the stored value (integer, rational or floating point).
  The compatibility is then checked in the trait method
  ``scalar_traits<ScalarType>::make_const(var_number const& vn)``.
- Added optional support for [boost::rational](
  https://www.boost.org/doc/libs/latest/libs/rational/rational.html) and
  [GMP](https://gmplib.org) C++ types ``mpz_class``, ``mpq_class``,
  ``mpf_class`` as ``ScalarType``.
- ``n_fermion_sector_view`` and ``n_fermion_multisector_view`` are now
  parameterized on the type of ranking algorithm used to map basis state indices
  from a full Hilbert space to a sector. Supported ranking algorithms are
  ``combination_ranking`` (selected by default), ``staggered_ranking`` and
  ``trie_ranking``. All three algorithms are described in M. Wallerberger,
  [K. Held, Phys. Rev. Research 4, 033238 (2022)](
  https://doi.org/10.1103/PhysRevResearch.4.033238).
- Reduced the maximum allowed number of bits in the binary representation of
  a basis state index (``hilbert_space::max_n_bits``) to 63. This way
  ``hilbert_space::dim()`` can return a valid value of type ``sv_index_type``
  even when all 63 bits are used up.
- New pure virtual method ``elementary_space::dim()`` and its implementation in
  derived classes.
- New method ``hilbert_space::is_sparse()`` that returns ``true`` if some of
  the constituent elementary spaces have non-power-of-two dimensions.
- New method ``hilbert_space::vec_size()`` that returns the minimal size of a
  state vector compatible with this Hilbert space.
- Semantics of the existing method ``hilbert_space::dim()`` has been changed:
  Now it returns the exact dimension of the Hilbert space, which is smaller than
  ``hilbert_space::vec_size()`` if the Hilbert space is sparse.
- New method overload ``hilbert_space::dim(elementary_space<...> const& es)``
  that returns dimension of a constituent elementary space.
- Changed ``space_partition`` to store a constant reference to the Hilbert space
  being partitioned. This change is necessary to enable support for the sparse
  ``hilbert_space`` objects. ``space_partition`` is now templated on the Hilbert
  space type, whereas its constructors are not.
- ``space_partition::merge_subspaces()`` and
  ``space_partition::find_connections()`` no longer accept the ``hs`` argument
  and use the stored reference instead.
- Two helper factory functions ``make_space_partition()`` have been added.
- Improved performance and stability of ``space_partition::merge_subspaces()``
  by switching to a non-recursive variant of the algorithm.
- Fixed a negative index bug in ``n_fermion_sector_view``.
  Credits to Dr. Cezary Śliwa for providing the patch.
- Whenever possible, use compiler intrinsics to speed up complex bit
  manipulation operations (``popcount``, ``tzcount``, ``pdep``, ``pext``).
- Change base type of ``spin_component`` to the 1-byte wide ``std::uint8_t``.
- New CMake option ``CPPCHECK_EXTRA_FLAGS``. It can be used to pass additional
  command line flags to ``cppcheck``.

## [0.7.2] - 2022-11-12

- Export namespaced CMake target ``libcommute::libcommute`` instead of
  plain ``libcommute``.
- Install CMake configuration files into
  ``${CMAKE_INSTALL_PREFIX}/lib/cmake/libcommute``, which is the recommended
  location.
- Upgraded bundled Catch2 to version 2.13.9 (this fixes issue #2 a.k.a.
  [catchorg/Catch2#2178](https://github.com/catchorg/Catch2/issues/2178)).
- Fixed compilation with clang/libc++ 15 (issue #5).
- Added project citation information.

## [0.7.1] - 2021-12-17

- New methods ``space_partition::subspace_basis()`` and
  ``space_partition::subspace_bases()``.
- New methods ``sparse_state_vector::prune()`` (two overloads).
- New example ``hubbard_holstein_1d`` and minor updates to the documentation.

## [0.7.0] - 2021-10-09

- Added an implementation of the ``StateVector`` concept for some
  [Eigen 3](https://eigen.tuxfamily.org/) types (vectors, vector segments,
  column-like matrix blocks and one-dimensional ``Eigen::Map`` views).
  The corresponding header file is ``loperator/state_vector_eigen3.hpp``.
- New classes ``n_fermion_sector_view`` and ``n_fermion_multisector_view`` that
  implement N-fermion (multi)sector views of a state vector. The classes are
  supplemented with free utility functions ``make_nf(m)s_view()``,
  ``make_const_nf(m)s_view()``, ``n_fermion_(multi)sector_size()`` and
  ``n_fermion_(multi)sector_basis_states()``.
- New member typedef ``hilbert_space::index_types``.
- New method ``hilbert_space::has_algebra()``.
- New method ``space_partition::find_connections()``.
- Renamed the ``constexpr`` integer ``LIBCOMMUTE_MIN_USER_DEFINED_ALGEBRA_ID``
  to ``min_user_defined_algebra_id`` so that it does not appear to be a macro.
- In the 3-argument constructor of ``basis_mapper``, change the type of the last
  argument ``N`` from ``int`` to ``unsigned int``.
- Removed methods ``basis_mapper::make_(const_)view_no_ref()``.
  ``basis_mapper::make_(const_)view()`` will now return a non-reference view
  when supplied with a non lvalue-reference argument.
- New CMake option ``STATIC_ANALYSIS``. When enabled, the ``clang-tidy`` and
  ``cppcheck`` static analysis tools will be run on the C++ sources of unit
  tests and examples as part of build process.
- Minor bugfixes in unit tests and improvements in coding style.
- Reformatted C++ sources using ``clang-format`` and a style based on ``LLVM``.

## [0.6.1] - 2021-03-30

- New method ``hilbert_space::index()``.
