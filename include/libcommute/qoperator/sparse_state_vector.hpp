/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_QOPERATOR_SPARSE_STATE_VECTOR_HPP_
#define LIBCOMMUTE_QOPERATOR_SPARSE_STATE_VECTOR_HPP_

#include "state_vector.hpp"
#include "../scalar_traits.hpp"

#include <unordered_map>

namespace libcommute {

//
// Implementation of the StateVector concept based on a sparse storage
//
template<typename ScalarType> class sparse_state_vector {

  sv_index_type size_;
  std::unordered_map<sv_index_type, ScalarType> data_;

public:

  sparse_state_vector() = delete;
  sparse_state_vector(sv_index_type size) : size_(size) {}

  // Dimension of Hilbert space
  inline friend sv_index_type get_size(sparse_state_vector const& sv) {
    return sv.size_;
  }

  // Number of non-zero amplitudes
  inline sv_index_type n_amplitudes() const { return data_.size(); }

  // Element access
  inline ScalarType const& amplitude(sv_index_type n) const {
    return data_[n];
  }
  inline ScalarType & amplitude(sv_index_type n) { return data_[n]; }

  // Get n-th state amplitude
  inline friend ScalarType
  get_element(sparse_state_vector const& sv, sv_index_type n) {
    assert(n < sv.size_);
    auto it = sv.data_.find(n);
    if(it == sv.data_.end())
      return scalar_traits<ScalarType>::make_const(0);
    else
      return it->second;
  }

  // Add a constant to the n-th state amplitude
  template<typename T>
  inline friend void update_add_element(sparse_state_vector & sv,
                                        sv_index_type n,
                                        T&& value) {
    auto it = sv.data_.find(n);
    if(it == sv.data_.end()) {
      if(!scalar_traits<ScalarType>::is_zero(value))
        sv.data_[n] = value;
    } else {
      it->second += value;
      if(scalar_traits<ScalarType>::is_zero(it->second))
        sv.data_.erase(it);
    }
  }

  // Set all amplitudes to zero
  inline friend void set_zeros(sparse_state_vector & sv) { sv.data_.clear(); }

  // Make an empty vector with the same Hilbert space dimension
  inline friend sparse_state_vector zeros_like(sparse_state_vector const& sv) {
    return sparse_state_vector(sv.size_);
  }

  // Apply functor `f` to all index/non-zero amplitude pairs
  template<typename Functor>
  inline friend void foreach(sparse_state_vector const& sv, Functor&& f) {
    for(auto const& p : sv.data_) f(p.first, p.second);
  }
};

// Get element type of a sparse_state_vector<ScalarType> object
template<typename ScalarType>
struct element_type<sparse_state_vector<ScalarType>> {
  using type = ScalarType;
};

} // namespace libcommute

#endif
