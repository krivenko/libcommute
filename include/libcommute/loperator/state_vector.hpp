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
#ifndef LIBCOMMUTE_LOPERATOR_STATE_VECTOR_HPP_
#define LIBCOMMUTE_LOPERATOR_STATE_VECTOR_HPP_

#include "../scalar_traits.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace libcommute {

// Type of index into a state vector
using sv_index_type = std::uint64_t;

// Bit width of the state index type
static constexpr int sv_index_width =
    std::numeric_limits<sv_index_type>::digits;

//
// Implementation of the StateVector interface for std::vector
//

// Get element type of a StateVector object
template <typename StateVector> struct element_type {};
template <typename StateVector>
using element_type_t = typename element_type<StateVector>::type;

template <typename T> struct element_type<std::vector<T>> {
  using type = T;
};
template <typename T> struct element_type<std::vector<T> const> {
  using type = T;
};

// Get n-th state amplitude stored in a standard vector
template <typename T>
inline auto get_element(std::vector<T> const& sv, sv_index_type n) -> T const& {
  return sv[n];
}

// Add a constant to the n-th state amplitude stored in a standard vector
template <typename T1, typename T2>
inline void
update_add_element(std::vector<T1>& sv, sv_index_type n, T2 const& value) {
  add_assign(sv[n], value);
}

// Set all amplitudes stored in a standard vector to zero
template <typename T> inline void set_zeros(std::vector<T>& sv) {
  std::fill(sv.begin(), sv.end(), scalar_traits<T>::make_const(0));
}

// Create a standard vector of the same size as `sv`,
// with all amplitudes set to zero
template <typename T>
inline std::vector<T> zeros_like(std::vector<T> const& sv) {
  return std::vector<T>(sv.size(), scalar_traits<T>::make_const(0));
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in a standard vector
template <typename T, typename Functor>
inline void foreach(std::vector<T> const& sv, Functor&& f) {
  for(sv_index_type n = 0; n < sv.size(); ++n) {
    auto const& a = sv[n];
    if(scalar_traits<T>::is_zero(a))
      continue;
    else
      std::forward<Functor>(f)(n, a);
  }
}

} // namespace libcommute

#endif
