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
#ifndef LIBCOMMUTE_LOPERATOR_SPARSE_STATE_VECTOR_HPP_
#define LIBCOMMUTE_LOPERATOR_SPARSE_STATE_VECTOR_HPP_

#include "../scalar_traits.hpp"
#include "state_vector.hpp"

#include <cassert>
#include <unordered_map>
#include <utility>

namespace libcommute {

//
// Implementation of the StateVector concept based on a sparse storage
//
template <typename ScalarType> class sparse_state_vector {

  sv_index_type size_;
  std::unordered_map<sv_index_type, ScalarType> data_;

public:
  sparse_state_vector() = delete;
  explicit sparse_state_vector(sv_index_type size) : size_(size) {}

  // Size of the vector
  inline sv_index_type size() const { return size_; }

  // Number of non-zero amplitudes
  inline sv_index_type n_nonzeros() const { return data_.size(); }

  // Element access
  inline ScalarType& operator[](sv_index_type n) { return data_[n]; }

  // Get n-th state amplitude
  inline friend ScalarType get_element(sparse_state_vector const& sv,
                                       sv_index_type n) {
    assert(n < sv.size_);
    auto it = sv.data_.find(n);
    if(it == sv.data_.end())
      return scalar_traits<ScalarType>::make_const(0);
    else
      return it->second;
  }

  // Add a constant to the n-th state amplitude
  template <typename T>
  inline friend void
  update_add_element(sparse_state_vector& sv, sv_index_type n, T&& value) {
    auto it = sv.data_.find(n);
    if(it == sv.data_.end()) {
      if(!scalar_traits<ScalarType>::is_zero(value))
        sv.data_[n] = std::forward<T>(value);
    } else {
      it->second += std::forward<T>(value);
      if(scalar_traits<ScalarType>::is_zero(it->second)) sv.data_.erase(it);
    }
  }

  // Set all amplitudes to zero
  inline friend void set_zeros(sparse_state_vector& sv) { sv.data_.clear(); }

  // Make an empty vector with the same Hilbert space dimension
  inline friend sparse_state_vector zeros_like(sparse_state_vector const& sv) {
    return sparse_state_vector(sv.size_);
  }

  // Apply functor `f` to all index/non-zero amplitude pairs
  template <typename Functor>
  inline friend void foreach(sparse_state_vector const& sv, Functor&& f) {
    for(auto const& p : sv.data_)
      std::forward<Functor>(f)(p.first, p.second);
  }

  // Force removal of all zero amplitudes from the storage
  inline void prune() {
    for(auto it = data_.begin(); it != data_.end();) {
      if(scalar_traits<ScalarType>::is_zero(it->second))
        it = data_.erase(it);
      else
        ++it;
    }
  }

  // Force removal of all amplitudes meeting a specified criterion
  template <typename UnaryPredicate> inline void prune(UnaryPredicate&& p) {
    for(auto it = data_.begin(); it != data_.end();) {
      if(std::forward<UnaryPredicate>(p)(it->second))
        it = data_.erase(it);
      else
        ++it;
    }
  }
};

// Get element type of a sparse_state_vector<ScalarType> object
template <typename ScalarType>
struct element_type<sparse_state_vector<ScalarType>> {
  using type = ScalarType;
};

} // namespace libcommute

#endif
