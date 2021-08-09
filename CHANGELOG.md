# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

- Added an implementation of the ``StateVector`` concept for some
  [Eigen 3](https://eigen.tuxfamily.org/) types (vectors, vector segments,
  column-like matrix blocks and one-dimensional ``Eigen::Map`` views).
  The corresponding header file is ``loperator/state_vector_eigen3.hpp``.
- New method ``space_partition.find_connections()``.
- New CMake option ``STATIC_ANALYSIS``. When enabled, the ``clang-tidy`` and
  ``cppcheck`` static analysis tools will be run on the C++ sources of unit
  tests and examples as part of build process.
- Minor bugfixes in unit tests and improvements in coding style.
- Reformatted C++ sources using ``clang-format`` and a style based on ``LLVM``.

## [0.6.1] - 2021-03-30

- New method ``hilbert_space.index()``.
