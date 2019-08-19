/*******************************************************************************
 *
 * This file is part of libcommute, a C++11/14/17 header-only library allowing
 * to manipulate polynomial expressions with quantum-mechanical operators.
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

// Get size of a StateVector object (Hilbert space dimension)
template<typename StateVector>
inline sv_index_type get_size(StateVector const& sv) {
  return sv.size();
}

// Get n-th state amplitude stored in a StateVector object
template<typename StateVector>
inline auto get_element(StateVector const& sv, sv_index_type n)
  -> element_type_t<StateVector> const& {
  return sv[n];
}

// Add a constant to the n-th state amplitude stored in a StateVector object
template<typename StateVector, typename T>
inline void update_add_element(StateVector & sv, sv_index_type n, T&& value) {
  sv[n] += value;
}

// Set all amplitudes stored in a StateVector object to zero
template<typename StateVector>
inline void set_zeros(StateVector & sv) {
  sv_index_type size = get_size(sv);
  using T = element_type_t<StateVector>;
  for(sv_index_type n = 0; n < size; ++n)
    sv[n] = scalar_traits<T>::make_const(0);
}

// Create a StateVector object of the same size as `sv`,
// with all amplitudes set to zero
template<typename StateVector>
StateVector zeros_like(StateVector const& sv) {
  StateVector res(get_size(sv));
  set_zeros(res);
  return res;
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in a StateVector object
template<typename StateVector, typename Functor>
inline void foreach(StateVector const& sv, Functor&& f) {
  sv_index_type size = get_size(sv);
  using T = element_type_t<StateVector>;
  for(sv_index_type n = 0; n < size; ++n) {
    auto a = get_element(sv, n);
    if(scalar_traits<T>::is_zero(a))
      continue;
    else
      f(n, a);
  }
}

} // namespace libcommute

#endif

