# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

- New method ``space_partition.find_connections()``.
- Added an implementation of the ``StateVector`` concept for some
  [Eigen 3](https://eigen.tuxfamily.org/) types (vectors, vector segments,
  column-like matrix blocks and one-dimensional ``Eigen::Map`` views).
  The corresponding header file is ``loperator/state_vector_eigen3.hpp``.

## [0.6.1] - 2021-03-30

- New method ``hilbert_space.index()``.
