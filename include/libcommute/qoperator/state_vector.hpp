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
#ifndef LIBCOMMUTE_QOPERATOR_STATE_VECTOR_HPP_
#define LIBCOMMUTE_QOPERATOR_STATE_VECTOR_HPP_

#include "../metafunctions.hpp"
#include "../scalar_traits.hpp"

#include <cstdint>
#include <type_traits>
#include <utility>

namespace libcommute {

// Type of index into a state vector
using sv_index_type = std::uint64_t;

//
// Default (unoptimal) implementation of StateVector interface
//

// Get element type of a StateVector object
template<typename StateVector>
struct element_type {
  using type = remove_cvref_t<decltype(std::declval<StateVector>()[0])>;
};
template<typename StateVector>
using element_type_t = typename element_type<StateVector>::type;

// Get n-th state amplitude stored in a StateVector object
template<typename StateVector>
inline auto get_element(StateVector const& sv, sv_index_type n)
  -> decltype(sv[n]) {
  return sv[n];
}

// Add a constant to the n-th state amplitude stored in a StateVector object
template<typename StateVector, typename T>
inline void update_add_element(StateVector & sv, sv_index_type n, T&& value) {
  auto & sv_n = sv[n];
  sv_n = sv_n + value;
}

// Set all amplitudes stored in a StateVector object to zero
template<typename StateVector>
inline void set_zeros(StateVector & sv) {
  sv_index_type size = sv.size();
  using T = element_type_t<StateVector>;
  for(sv_index_type n = 0; n < size; ++n)
    sv[n] = scalar_traits<T>::make_const(0);
}

// Create a StateVector object of the same size as `sv`,
// with all amplitudes set to zero
template<typename StateVector>
StateVector zeros_like(StateVector const& sv) {
  StateVector res(sv.size());
  set_zeros(res);
  return res;
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in a StateVector object
template<typename StateVector, typename Functor>
inline void foreach(StateVector const& sv, Functor&& f) {
  sv_index_type size = sv.size();
  using T = element_type_t<StateVector>;
  for(sv_index_type n = 0; n < size; ++n) {
    // Emulate decltype(auto)
    decltype(get_element(sv, n)) a = get_element(sv, n);
    if(scalar_traits<T>::is_zero(a))
      continue;
    else
      f(n, a);
  }
}

} // namespace libcommute

#endif

