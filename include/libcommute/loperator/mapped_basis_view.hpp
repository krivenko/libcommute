/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_LOPERATOR_MAPPED_BASIS_VIEW_HPP_
#define LIBCOMMUTE_LOPERATOR_MAPPED_BASIS_VIEW_HPP_

#include "loperator.hpp"
#include "sparse_state_vector.hpp"
#include "state_vector.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

//
// mapped_basis_view partially models the StateVector concept while adapting
// an existing state vector. Such operations as `get_element()` and
// `update_add_element()` are forwarded to the adapted vector with state
// indices being translated according to a given map. mapped_basis_view can be
// used in situations where a loperator object acts only on a subspace of a full
// Hilbert space and it is desirable to store vector components only within this
// subspace.
//

namespace libcommute {

template<typename StateVector, bool Const = false>
struct mapped_basis_view {

  typename std::conditional<Const,
                            StateVector const&,
                            StateVector&>::type state_vector;
  using map_t = std::unordered_map<sv_index_type, sv_index_type>;
  map_t const& map;

  using scalar_type = typename element_type<StateVector>::type;

  mapped_basis_view(decltype(state_vector) sv, map_t const& map)
    : state_vector(sv), map(map)
  {}
};

// Get element type of the StateVector object adapted by a given
// mapped_basis_view object.
template<typename StateVector, bool Const>
struct element_type<mapped_basis_view<StateVector, Const>> {
  using type = typename mapped_basis_view<StateVector, Const>::scalar_type;
};

// Get state amplitude of the adapted StateVector object at index view.map[n]
template<typename StateVector, bool Const>
inline auto get_element(mapped_basis_view<StateVector, Const> const& view,
                        sv_index_type n)
  -> typename mapped_basis_view<StateVector, Const>::scalar_type {
  return get_element(view.state_vector, view.map.at(n));
}

// Add a constant to a state amplitude stored in the adapted StateVector object
// at index view.map[n].
template<typename StateVector, typename T>
inline
void update_add_element(mapped_basis_view<StateVector, false> & view,
                        sv_index_type n,
                        T&& value) {
  update_add_element(view.state_vector, view.map.at(n), std::forward<T>(value));
}

// update_add_element() is not defined for constant views
template<typename StateVector, typename T>
inline
void update_add_element(mapped_basis_view<StateVector, true> &,
                        sv_index_type,
                        T&&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "update_add_element() is not supported for constant views");
}

// zeros_like() is not defined for views
template<typename StateVector, bool Const>
StateVector zeros_like(mapped_basis_view<StateVector, Const> const&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "zeros_like() is not supported for views");
}

// Set all amplitudes stored in the adapted StateVector object to zero
template<typename StateVector>
inline void set_zeros(mapped_basis_view<StateVector, false> & view) {
  set_zeros(view.state_vector);
}

// set_zeros() is not defined for constant views
template<typename StateVector>
inline void set_zeros(mapped_basis_view<StateVector, true> &) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "set_zeros() is not supported for constant views");
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in the adapted StateVector object. This functions iterates over all values
// stored in view.map.
template<typename StateVector, bool Const, typename Functor>
inline void foreach(mapped_basis_view<StateVector, Const> const& view,
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

//
// Factory class for mapped_basis_view
//
class basis_mapper {

  std::unordered_map<sv_index_type, sv_index_type> map_;

  template<typename LOpScalarType, int... LOpAlgebraIDs>
  void compositions_constructor_impl(
    std::vector<loperator<LOpScalarType, LOpAlgebraIDs...>> const& O_list,
    sparse_state_vector<
      typename loperator<LOpScalarType, LOpAlgebraIDs...>::scalar_type
    > const& st,
    int m,
    int sum_n,
    int N) {

    if(sum_n == N) {
      foreach(st, [&](
        sv_index_type out_index,
        typename loperator<LOpScalarType, LOpAlgebraIDs...>::scalar_type const&
      ) {
        map_.emplace(out_index, map_.size());
      });
    }

    if(std::size_t(m) == O_list.size()) return;

    compositions_constructor_impl(O_list, st, m + 1, sum_n, N);

    auto from_st = st;
    auto to_st = zeros_like(st);
    for(; sum_n < N; ++sum_n) {
      O_list[m](from_st, to_st);
      if(to_st.n_nonzeros() == 0) return;
      compositions_constructor_impl(O_list, to_st, m + 1, sum_n + 1, N);
      from_st = to_st;
    }
  }

public:

  // Build a mapping from a list of basis states
  // to their positions within the list
  basis_mapper(std::vector<sv_index_type> const& basis_state_indices) {
    std::transform(
      basis_state_indices.begin(),
      basis_state_indices.end(),
      std::inserter(map_, map_.end()),
      [this](sv_index_type n) { return std::make_pair(n, map_.size()); }
    );
  }

  // Build a mapping from a set of all basis states contributing to O|vac>.
  // Mapped values are assigned continuously but without any specific order.
  template<typename HSType, typename LOpScalarType, int... LOpAlgebraIDs>
  basis_mapper(loperator<LOpScalarType, LOpAlgebraIDs...> const& O,
               HSType const& hs) {
    using scalar_type =
      typename loperator<LOpScalarType, LOpAlgebraIDs...>::scalar_type;
    sv_index_type dim = get_dim(hs);
    sparse_state_vector<scalar_type> vac(dim);
    vac[0] = 1;
    auto st =  O(vac);
    foreach(st, [&](sv_index_type out_index, scalar_type const&) {
      map_.emplace(out_index, map_.size());
    });
  }

  // Given a list of operators {O_1, O_2, O_3, ... , O_M}, build a mapping
  // including all basis states contributing to all states
  // O_1^{n_1} O_2^{n_2} ... O_M^{n_M} |vac>, where n_m >= 0 and
  // \sum_{m=1}^M n_M = N.
  // Mapped values are assigned continuously but without any specific order.
  template<typename HSType, typename LOpScalarType, int... LOpAlgebraIDs>
  basis_mapper(
    std::vector<loperator<LOpScalarType, LOpAlgebraIDs...>> const& O_list,
    HSType const& hs,
    int N
  ) {
    if(N == 0 || O_list.size() == 0) {
      map_.emplace(0, 0);
      return;
    }
    using scalar_type =
      typename loperator<LOpScalarType, LOpAlgebraIDs...>::scalar_type;
    sparse_state_vector<scalar_type> vac(get_dim(hs));
    vac[0] = 1;
    compositions_constructor_impl(O_list, vac, 0, 0, N);
  }

  // Number of basis states in the mapping
  inline sv_index_type size() const { return map_.size(); }

  // Direct access to the mapping
  inline std::unordered_map<sv_index_type, sv_index_type> const& map() const {
    return map_;
  }

  // Direct access to the inverse mapping (slow!)
  inline std::unordered_map<sv_index_type, sv_index_type> inverse_map() const {
    std::unordered_map<sv_index_type, sv_index_type> inv_map;
    std::transform(
      map_.begin(),
      map_.end(),
      std::inserter(inv_map, inv_map.end()),
      [](std::pair<sv_index_type, sv_index_type> const& p) {
        return std::make_pair(p.second, p.first);
      }
    );
    // Check that map_ is actually invertible
    assert(inv_map.size() == map_.size());
    return inv_map;
  }

  // Make a non-constant basis mapping view
  template<typename StateVector>
  mapped_basis_view<StateVector, false> make_view(StateVector & sv) const {
    return mapped_basis_view<StateVector, false>(sv, map_);
  }

  // Make a constant basis mapping view
  template<typename StateVector>
  mapped_basis_view<StateVector, true>
  make_const_view(StateVector const& sv) const {
    return mapped_basis_view<StateVector, true>(sv, map_);
  }
};

} // namespace libcommute

#endif
