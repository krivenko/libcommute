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
#ifndef LIBCOMMUTE_QOPERATOR_SPACE_PARTITION_HPP_
#define LIBCOMMUTE_QOPERATOR_SPACE_PARTITION_HPP_

#include "disjoint_sets.hpp"
#include "sparse_state_vector.hpp"
#include "qoperator.hpp"
#include "../scalar_traits.hpp"

#include <map>
#include <stdexcept>
#include <utility>

//
// Implementation of the automatic partitioning algorithm
//

namespace libcommute {

// Partition of a Hilbert space into a set of disjoint subspaces invariant
// under action of a given Hermitian operator (Hamiltonian).
//
/* For a detailed description of the algorithm see
 * `Computer Physics Communications 200, March 2016, 274-284
 * <http://dx.doi.org/10.1016/j.cpc.2015.10.023>`_ (section 4.2).
*/
class space_partition {

  // Space partition
  disjoint_sets ds;
  // Map representative basis state to subspace index
  std::map<sv_index_type, sv_index_type> root_to_subspace;

public:

  template<typename ScalarType>
  using matrix_elements_type = std::map<std::pair<sv_index_type, sv_index_type>,
                                        ScalarType>;

  space_partition() = delete;

  // Partition Hilbert space `hs` using Hermitian operator `h`.
  //
  // `hs` can be of any type, for which `get_dim(hs)` returns the dimension of
  // the corresponding Hilbert space, and `foreach(hs, f)` applies functor `f`
  // to each basis state index in `hs`.
  template<typename HSType, typename... QOperatorParams>
  space_partition(qoperator<QOperatorParams...> const& h, HSType const& hs)
    : ds(get_dim(hs)) {
    using scalar_type = typename qoperator<QOperatorParams...>::scalar_type;
    sv_index_type dim = get_dim(hs);

    sparse_state_vector<scalar_type> in_state(dim);
    sparse_state_vector<scalar_type> out_state(dim);
    foreach(hs, [&](sv_index_type in_index) {
      in_state.amplitude(in_index) = scalar_traits<scalar_type>::make_const(1);
      h(in_state, out_state);
      foreach(out_state, [&](sv_index_type out_index, scalar_type const&) {
        ds.set_union(in_index, out_index);
      });
      set_zeros(in_state);
    });

    update_root_to_subspace();
  }

  // Partition Hilbert space `hs` using Hermitian operator `h`. Save
  // non-vanishing matrix elements of `h` into the sparse storage object
  // `matrix_elements`.
  //
  // `hs` can be of any type, for which `get_dim(hs)` returns the dimension of
  // the corresponding Hilbert space, and `foreach(hs, f)` applies functor `f`
  // to each basis state index in `hs`.
  template<typename HSType, typename... QOperatorParams>
  space_partition(qoperator<QOperatorParams...> const& h,
                  HSType const& hs,
                  matrix_elements_type<
                    typename qoperator<QOperatorParams...>::scalar_type
                  > & me)
    : ds(get_dim(hs)) {
    using scalar_type = typename qoperator<QOperatorParams...>::scalar_type;
    sv_index_type dim = get_dim(hs);

    sparse_state_vector<scalar_type> in_state(dim);
    sparse_state_vector<scalar_type> out_state(dim);
    foreach(hs, [&](sv_index_type in_index) {
      in_state.amplitude(in_index) = scalar_traits<scalar_type>::make_const(1);
      h(in_state, out_state);

      foreach(out_state, [&](sv_index_type out_index, scalar_type const& a) {
        ds.set_union(in_index, out_index);
        me[std::make_pair(out_index, in_index)] = a;
      });

      set_zeros(in_state);
    });

    update_root_to_subspace();
  }

  // Hilbert space dimension
  sv_index_type dim() const { return ds.size(); }

  // Number of subspaces
  sv_index_type n_subspaces() const { return ds.n_sets(); }

  // Find what invariant subspace a given basis state belongs to
  sv_index_type operator[](sv_index_type index) const {
    auto it = root_to_subspace.find(ds.find_root(index));
    if(it == root_to_subspace.end())
      throw std::runtime_error("Unexpected basis state index " +
                               std::to_string(index));
    else
      return it->second;
  }

 // Apply a functor `f` to all basis states in a given space partition.
 // The functor must take two arguments, index of the basis state,
 // and index of the subspace this basis state belongs to.
 template<typename F>
 friend void foreach(space_partition const& sp, F&& f) {
   for(sv_index_type n = 0; n < sp.dim(); ++n) { f(n, sp[n]); }
 }

private:

  void update_root_to_subspace() {
    ds.compress_sets();
    ds.normalize_sets();

    root_to_subspace.clear();
    for(sv_index_type n = 0; n < dim(); ++n) {
      root_to_subspace.emplace(ds.find_root(n), root_to_subspace.size());
    }
  }

};

} // namespace libcommute

#endif
