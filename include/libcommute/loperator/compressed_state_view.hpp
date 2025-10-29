/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2025 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_LOPERATOR_COMPRESSED_STATE_VIEW_HPP_
#define LIBCOMMUTE_LOPERATOR_COMPRESSED_STATE_VIEW_HPP_

#include "../metafunctions.hpp"
#include "bit_ops.hpp"
#include "hilbert_space.hpp"
#include "state_vector.hpp"

#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>

//
// compressed_state_view partially models the StateVector concept while adapting
// an existing state vector. Such operations as `get_element()` and
// `update_add_element()` are forwarded to the adapted vector. Basis state
// indices are translated according to a map from a Hilbert space to a
// continuous range, which means that the adapted state vector must store
// exactly hilbert_space::dim() amplitudes. For the sparse Hilbert spaces,
// this can lead to dramatic memory saving.
//

namespace libcommute {

template <typename StateVector, bool Ref = true> struct compressed_state_view {

  typename std::conditional<Ref, StateVector&, StateVector>::type state_vector;

  // Bit extraction masks to be used in map_index()
  std::vector<sv_index_type> masks;
  // Strides to be used in map_index()
  std::vector<sv_index_type> strides;

  // Implementation of foreach() algorithm
  detail::sparse_foreach_basis_state foreach_alg;

  using scalar_type = typename element_type<
      typename std::remove_const<StateVector>::type>::type;

  template <typename SV, typename HSType>
  compressed_state_view(SV&& sv, HSType const& hs)
    : state_vector(std::forward<SV>(sv)),
      foreach_alg(init(hs, masks, strides)) {}

  inline sv_index_type map_index(sv_index_type index) const {
    sv_index_type res = 0;
    for(std::size_t d = 0; d < strides.size(); ++d) {
      res += detail::extract_bits(index, masks[d]) * strides[d];
    }
    return res;
  }

private:
  template <typename HSType>
  static detail::sparse_foreach_basis_state
  init(HSType const& hs,
       std::vector<sv_index_type>& masks,
       std::vector<sv_index_type>& strides) {
    std::vector<sv_index_type> sizes;
    sizes.reserve(hs.size());

    // Fill 'sizes' and 'masks' while merging (replacing by the product)
    // dimensions of a non-power-of-two elementary space and all power-of-two
    // spaces preceding it.
    sv_index_type dim = 1;
    sv_index_type mask = 0;
    hs.foreach_elementary_space(
        [&hs, &sizes, &masks, &dim, &mask](
            typename HSType::elementary_space_t const& es) {
          int n_bits = es.n_bits();
          sv_index_type es_dim = es.dim();

          dim *= es_dim;
          mask += ((sv_index_type(1) << n_bits) - 1) << hs.bit_range(es).first;

          // Non-power-of-two elementary space: Add merged dimension to 'sizes'
          // and a new extraction mask to 'masks'
          if(es_dim != detail::pow2(n_bits)) {
            sizes.push_back(dim);
            dim = 1;
            masks.push_back(mask);
            mask = 0;
          }
        });
    if(dim != 1) {
      sizes.push_back(dim);
      masks.push_back(mask);
    }

    // Compute strides from sizes
    if(!sizes.empty()) {
      strides.resize(sizes.size());
      strides[0] = 1;
      for(std::size_t d = 1; d < strides.size(); ++d) {
        strides[d] = strides[d - 1] * sizes[d - 1];
      }
    }

    return detail::sparse_foreach_basis_state(sizes);
  }
};

// Get element type of the StateVector object adapted by a given
// compressed_state_view object.
template <typename StateVector, bool Ref>
struct element_type<compressed_state_view<StateVector, Ref>> {
  using type = typename compressed_state_view<StateVector>::scalar_type;
};

// Get state amplitude of the adapted StateVector object
template <typename StateVector, bool Ref>
inline auto get_element(compressed_state_view<StateVector, Ref> const& view,
                        sv_index_type n) ->
    typename compressed_state_view<StateVector>::scalar_type {
  return get_element(view.state_vector, view.map_index(n));
}

// Add a constant to a state amplitude stored in the adapted StateVector object
template <typename StateVector, bool Ref, typename T>
inline void update_add_element(compressed_state_view<StateVector, Ref>& view,
                               sv_index_type n,
                               T&& value) {
  update_add_element(view.state_vector,
                     view.map_index(n),
                     std::forward<T>(value));
}

// update_add_element() is not defined for constant views
template <typename StateVector, bool Ref, typename T>
inline void
update_add_element(compressed_state_view<StateVector const, Ref>&,
                   sv_index_type,
                   T&&) { // NOLINT(cppcoreguidelines-missing-std-forward)
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "update_add_element() is not supported for constant views");
}

// zeros_like() is not defined for views
template <typename StateVector, bool Ref>
inline StateVector zeros_like(compressed_state_view<StateVector, Ref> const&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "zeros_like() is not supported for views");
}

// Set all amplitudes stored in the adapted StateVector object to zero
template <typename StateVector, bool Ref>
inline void set_zeros(compressed_state_view<StateVector, Ref>& view) {
  set_zeros(view.state_vector);
}

// set_zeros() is not defined for constant views
template <typename StateVector, bool Ref>
inline void set_zeros(compressed_state_view<StateVector const, Ref>&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "set_zeros() is not supported for constant views");
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in the adapted StateVector object.
template <typename StateVector, bool Ref, typename Functor>
inline void foreach(compressed_state_view<StateVector, Ref> const& view,
                    Functor&& f) {
  using T = typename compressed_state_view<StateVector, Ref>::scalar_type;

  sv_index_type n = 0;
  view.foreach_alg([&view, &f, &n](sv_index_type st) {
    // Emulate decltype(auto)
    decltype(get_element(view.state_vector, n)) a =
        get_element(view.state_vector, n);
    ++n;

    if(scalar_traits<T>::is_zero(a))
      return;
    else
      std::forward<Functor>(f)(st, a);
  });
}

// Make a non-constant compressed state view
template <typename StateVector, typename HSType>
auto make_comp_state_view(StateVector&& sv, HSType const& hs)
    -> compressed_state_view<remove_cvref_t<StateVector>,
                             std::is_lvalue_reference<StateVector>::value> {
  return {std::forward<StateVector>(sv), hs};
}

// Make a constant compressed state view
template <typename StateVector, typename HSType>
auto make_const_comp_state_view(StateVector&& sv, HSType const& hs)
    -> compressed_state_view<remove_cvref_t<StateVector> const,
                             std::is_lvalue_reference<StateVector>::value> {
  return {std::forward<StateVector>(sv), hs};
}

} // namespace libcommute

#endif
