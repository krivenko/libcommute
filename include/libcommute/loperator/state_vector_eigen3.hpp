/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_LOPERATOR_STATE_VECTOR_EIGEN3_HPP_
#define LIBCOMMUTE_LOPERATOR_STATE_VECTOR_EIGEN3_HPP_

#include "../scalar_traits.hpp"
#include "state_vector.hpp"

#include <Eigen/Core>

#include <cassert>
#include <type_traits>

namespace libcommute {

//
// Implementation of the StateVector interface for Eigen3 types
//

//
// element_type specializations
//

// Vector types

template <typename ScalarType, int Rows, int Opts, int MRows>
struct element_type<Eigen::Matrix<ScalarType, Rows, 1, Opts, MRows, 1>> {
  using type = ScalarType;
};

template <typename ScalarType, int Rows, int Opts, int MRows>
struct element_type<Eigen::Matrix<ScalarType, Rows, 1, Opts, MRows, 1> const> {
  using type = ScalarType;
};

// Eigen::VectorBlock (vector segment)

template <typename XprType> struct element_type<Eigen::VectorBlock<XprType>> {
  using type = typename XprType::Scalar;
};

template <typename XprType>
struct element_type<Eigen::VectorBlock<XprType> const> {
  using type = typename XprType::Scalar;
};

// Column of a matrix

template <typename XprType, int BlockRows, bool InnerPanel>
struct element_type<Eigen::Block<XprType, BlockRows, 1, InnerPanel>> {
  using type = typename XprType::Scalar;
};

template <typename XprType, int BlockRows, bool InnerPanel>
struct element_type<Eigen::Block<XprType, BlockRows, 1, InnerPanel> const> {
  using type = typename XprType::Scalar;
};

// Nx1-submatrix of a matrix

template <typename XprType, int BlockRows, bool InnerPanel>
struct element_type<
    Eigen::Block<XprType, BlockRows, Eigen::Dynamic, InnerPanel>> {
  using type = typename XprType::Scalar;
};

template <typename XprType, int BlockRows, bool InnerPanel>
struct element_type<
    Eigen::Block<XprType, BlockRows, Eigen::Dynamic, InnerPanel> const> {
  using type = typename XprType::Scalar;
};

// Vector-like Eigen::Map view

template <typename ScalarType, int Rows, int Opts, int MRows>
struct element_type<
    Eigen::Map<Eigen::Matrix<ScalarType, Rows, 1, Opts, MRows, 1>>> {
  using type = ScalarType;
};

template <typename ScalarType, int Rows, int Opts, int MRows>
struct element_type<
    Eigen::Map<Eigen::Matrix<ScalarType, Rows, 1, Opts, MRows, 1>> const> {
  using type = ScalarType;
};

//
// Specializations of free functions for Eigen::DenseBase
//

template <typename Derived>
inline auto get_element(Eigen::DenseBase<Derived> const& sv, sv_index_type n) ->
    typename Eigen::DenseBase<Derived>::Scalar {
  return sv(n);
}

template <typename Derived, typename T>
inline void
update_add_element(Eigen::DenseBase<Derived>& sv, sv_index_type n, T value) {
  sv(n) += value;
}

template <typename Derived>
inline void set_zeros(Eigen::DenseBase<Derived>& sv) {
  sv.setZero();
}

template <typename Derived>
inline Eigen::
    Matrix<typename Derived::Scalar, Eigen::Dynamic, 1, 0, Eigen::Dynamic, 1>
    zeros_like(Eigen::DenseBase<Derived> const& sv) {
  using ret_t = Eigen::
      Matrix<typename Derived::Scalar, Eigen::Dynamic, 1, 0, Eigen::Dynamic, 1>;
  return ret_t::Zero(sv.size());
}

template <typename Derived, typename Functor>
inline void foreach(Eigen::DenseBase<Derived> const& sv, Functor&& f) {
  using ScalarType = typename Eigen::DenseBase<Derived>::Scalar;
  sv_index_type size = sv.size();
  for(sv_index_type n = 0; n < size; ++n) {
    auto const& a = sv(n);
    if(scalar_traits<ScalarType>::is_zero(a))
      continue;
    else
      f(n, a);
  }
}

//
// Specializations of free functions for Eigen::Block with
// a dynamic number of columns
//

template <typename XprType, int BlockRows, bool InnerPanel>
inline auto get_element(
    Eigen::Block<XprType, BlockRows, Eigen::Dynamic, InnerPanel> const& sv,
    sv_index_type n) -> typename XprType::Scalar {
  assert(sv.cols() == 1);
  return sv(n, 0);
}

template <typename XprType, int BlockRows, bool InnerPanel, typename T>
inline void update_add_element(
    Eigen::Block<XprType, BlockRows, Eigen::Dynamic, InnerPanel>& sv,
    sv_index_type n,
    T value) {
  assert(sv.cols() == 1);
  sv(n, 0) += value;
}

template <typename XprType, int BlockRows, bool InnerPanel, typename Functor>
inline void
foreach(Eigen::Block<XprType, BlockRows, Eigen::Dynamic, InnerPanel> const& sv,
        Functor&& f) {
  using ScalarType = typename XprType::Scalar;
  assert(sv.cols() == 1);
  sv_index_type size = sv.rows();
  for(sv_index_type n = 0; n < size; ++n) {
    auto const& a = sv(n, 0);
    if(scalar_traits<ScalarType>::is_zero(a))
      continue;
    else
      f(n, a);
  }
}

} // namespace libcommute

#endif
