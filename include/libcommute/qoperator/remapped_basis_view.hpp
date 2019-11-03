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
#ifndef LIBCOMMUTE_QOPERATOR_REMAPPED_BASIS_VIEW_HPP_
#define LIBCOMMUTE_QOPERATOR_REMAPPED_BASIS_VIEW_HPP_

#include "state_vector.hpp"

#include <type_traits>
#include <unordered_map>
#include <utility>

//
// remapped_basis_view partially models the StateVector concept while adapting
// an existing state vector. Such operations as `get_element()` and
// `update_add_element()` are forwarded to the adapted vector with state
// indices being translated according to a given map. remapped_basis_view can be
// used in situations where a qoperator object acts only on a subspace of a full
// Hilbert space and it is desirable to store vector components only within this
// subspace.
//

namespace libcommute {

template<typename StateVector, bool Const = false>
struct remapped_basis_view {

  typename std::conditional<Const,
                            StateVector const&,
                            StateVector&>::type state_vector;
  using map_t = std::unordered_map<sv_index_type, sv_index_type>;
  map_t const& map;

  using scalar_type = typename element_type<StateVector>::type;

  remapped_basis_view(decltype(state_vector) sv, map_t const& map)
    : state_vector(sv), map(map)
  {}
};

// Get element type of the StateVector object adapted by a given
// remapped_basis_view object.
template<typename StateVector, bool Const>
struct element_type<remapped_basis_view<StateVector, Const>> {
  using type = typename remapped_basis_view<StateVector, Const>::scalar_type;
};

// Get state amplitude of the adapted StateVector object at index view.map[n]
template<typename StateVector, bool Const>
inline auto get_element(remapped_basis_view<StateVector, Const> const& view,
                        sv_index_type n)
  -> typename remapped_basis_view<StateVector, Const>::scalar_type {
  return get_element(view.state_vector, view.map.at(n));
}

// Add a constant to a state amplitude stored in the adapted StateVector object
// at index view.map[n].
template<typename StateVector, typename T>
inline
void update_add_element(remapped_basis_view<StateVector, false> & view,
                        sv_index_type n,
                        T&& value) {
  update_add_element(view.state_vector, view.map.at(n), std::forward<T>(value));
}

// update_add_element() is not defined for constant views
template<typename StateVector, typename T>
inline
void update_add_element(remapped_basis_view<StateVector, true> &,
                        sv_index_type,
                        T&&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "update_add_element() is not supported for constant views");
}

// zeros_like() is not defined for views
template<typename StateVector, bool Const>
StateVector zeros_like(remapped_basis_view<StateVector, Const> const&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "zeros_like() is not supported for views");
}

// Set all amplitudes stored in the adapted StateVector object to zero
template<typename StateVector>
inline void set_zeros(remapped_basis_view<StateVector, false> & view) {
  set_zeros(view.state_vector);
}

// set_zeros() is not defined for constant views
template<typename StateVector>
inline void set_zeros(remapped_basis_view<StateVector, true> &) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "set_zeros() is not supported for constant views");
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in the adapted StateVector object. This functions iterates over all values
// stored in view.map.
template<typename StateVector, bool Const, typename Functor>
inline void foreach(remapped_basis_view<StateVector, Const> const& view,
                    Functor&& f) {
  using T = element_type_t<StateVector>;
  for(auto const& p : view.map) {
    // Emulate decltype(auto)
    decltype(get_element(view.state_vector, p.second)) a =
      get_element(view.state_vector, p.second);
    if(scalar_traits<T>::is_zero(a))
      continue;
    else
      f(p.first, a);
  }
}

} // namespace libcommute

#endif
