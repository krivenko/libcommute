# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

- Added an implementation of the ``StateVector`` concept for some
  [Eigen 3](https://eigen.tuxfamily.org/) types (vectors, vector segments,
  column-like matrix blocks and one-dimensional ``Eigen::Map`` views).
  The corresponding header file is ``loperator/state_vector_eigen3.hpp``.
- New method ``hilbert_space.has_algebra()``.
- New method ``space_partition.find_connections()``.
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

- New method ``hilbert_space.index()``.
