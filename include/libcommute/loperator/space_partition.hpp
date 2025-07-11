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
#ifndef LIBCOMMUTE_LOPERATOR_SPACE_PARTITION_HPP_
#define LIBCOMMUTE_LOPERATOR_SPACE_PARTITION_HPP_

#include "../scalar_traits.hpp"
#include "disjoint_sets.hpp"
#include "loperator.hpp"
#include "sparse_state_vector.hpp"

#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

//
// Implementation of the automatic partitioning algorithm
//

namespace libcommute {

// Sparse storage for matrix elements of a quantum operator
template <typename ScalarType>
using matrix_elements_map =
    std::map<std::pair<sv_index_type, sv_index_type>, ScalarType>;

// Connections between subspaces
using connections_map = std::set<std::pair<sv_index_type, sv_index_type>>;

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

  template <typename LOpScalarType, int... LOpAlgebraIDs>
  using loperator_melem_t = matrix_elements_map<
      typename loperator<LOpScalarType, LOpAlgebraIDs...>::scalar_type>;

public:
  space_partition() = delete;

  // Partition Hilbert space `hs` using Hermitian operator `h`.
  //
  // `hs` can be of any type, for which `get_dim(hs)` returns the dimension of
  // the corresponding Hilbert space, and `foreach(hs, f)` applies functor `f`
  // to each basis state index in `hs`.
  template <typename HSType, typename LOpScalarType, int... LOpAlgebraIDs>
  space_partition(loperator<LOpScalarType, LOpAlgebraIDs...> const& h,
                  HSType const& hs)
    : ds(get_dim(hs)) {
    using scalar_type =
        typename loperator<LOpScalarType, LOpAlgebraIDs...>::scalar_type;
    sv_index_type d = get_dim(hs);

    sparse_state_vector<scalar_type> in_state(d);
    sparse_state_vector<scalar_type> out_state(d);
    foreach(hs, [&](sv_index_type in_index) {
      in_state[in_index] = scalar_traits<scalar_type>::make_const(1);
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
  template <typename HSType, typename LOpScalarType, int... LOpAlgebraIDs>
  space_partition(loperator<LOpScalarType, LOpAlgebraIDs...> const& h,
                  HSType const& hs,
                  loperator_melem_t<LOpScalarType, LOpAlgebraIDs...>& me)
    : ds(get_dim(hs)) {
    using scalar_type =
        typename loperator<LOpScalarType, LOpAlgebraIDs...>::scalar_type;
    sv_index_type d = get_dim(hs);

    sparse_state_vector<scalar_type> in_state(d);
    sparse_state_vector<scalar_type> out_state(d);
    foreach(hs, [&](sv_index_type in_index) {
      in_state[in_index] = scalar_traits<scalar_type>::make_const(1);
      h(in_state, out_state);

      foreach(out_state, [&](sv_index_type out_index, scalar_type const& a) {
        ds.set_union(in_index, out_index);
        me[std::make_pair(out_index, in_index)] = a;
      });

      set_zeros(in_state);
    });

    update_root_to_subspace();
  }

  // Perform Phase II of the automatic partition algorithm
  //
  // Merge some of the invariant subspaces together, to ensure that a given
  // operator `Cd` and its Hermitian conjugate `C` generate only one-to-one
  // connections between the subspaces.
  template <typename HSType, typename LOpScalarType, int... LOpAlgebraIDs>
  auto merge_subspaces(loperator<LOpScalarType, LOpAlgebraIDs...> const& Cd,
                       loperator<LOpScalarType, LOpAlgebraIDs...> const& C,
                       HSType const& hs,
                       bool store_matrix_elements = true)
      -> std::pair<loperator_melem_t<LOpScalarType, LOpAlgebraIDs...>,
                   loperator_melem_t<LOpScalarType, LOpAlgebraIDs...>>

  {
    using loperator_t = loperator<LOpScalarType, LOpAlgebraIDs...>;
    using scalar_type = typename loperator_t::scalar_type;
    matrix_elements_map<scalar_type> Cd_elements, C_elements;

    using conn_t = std::multimap<sv_index_type, sv_index_type>;
    conn_t Cd_conn, C_conn;

    sv_index_type d = get_dim(hs);

    // Fill connection multimaps
    sparse_state_vector<scalar_type> in_state(d);
    sparse_state_vector<scalar_type> out_state(d);
    foreach(hs, [&](sv_index_type in_index) {
      in_state[in_index] = scalar_traits<scalar_type>::make_const(1);
      sv_index_type in_subspace = ds.find_root(in_index);

      auto fill_conn = [&, this](loperator_t const& lop,
                                 conn_t& conn,
                                 matrix_elements_map<scalar_type>& elem) {
        lop(in_state, out_state);
        // Iterate over non-zero final amplitudes
        foreach(out_state, [&](sv_index_type out_index, scalar_type const& a) {
          sv_index_type out_subspace = ds.find_root(out_index);
          conn.insert({in_subspace, out_subspace});
          if(store_matrix_elements) elem[{out_index, in_index}] = a;
        });
      };

      fill_conn(Cd, Cd_conn, Cd_elements);
      fill_conn(C, C_conn, C_elements);

      set_zeros(in_state);
    });

    // Merge all 'out' subspaces corresponding to the same 'in' subspace
    // in 'conn'.
    auto merge_conn_targets = [this](conn_t const& conn) -> bool {
      if(conn.empty()) return false;

      bool subspaces_linked = false;

      auto conn_it = conn.cbegin();
      sv_index_type in_subspace = ds.find_root(conn_it->first);
      sv_index_type out_subspace = ds.find_root(conn_it->second);
      ++conn_it;
      for(; conn_it != conn.cend(); ++conn_it) {
        if(ds.find_root(conn_it->first) == ds.find_root(in_subspace)) {
          sv_index_type out_subspace1 = ds.find_root(out_subspace);
          sv_index_type out_subspace2 = ds.find_root(conn_it->second);
          if(out_subspace1 != out_subspace2) {
            ds.root_union(out_subspace1, out_subspace2);
            subspaces_linked = true;
          }
        } else {
          in_subspace = ds.find_root(conn_it->first);
          out_subspace = ds.find_root(conn_it->second);
        }
      }

      return subspaces_linked;
    };

    // Repeatedly call merge_conn_targets() until no further subspaces can
    // be merged.
    // NOLINTNEXTLINE(cppcoreguidelines-init-variables)
    bool sl_Cd, sl_C;
    do {
      sl_Cd = merge_conn_targets(Cd_conn);
      sl_C = merge_conn_targets(C_conn);
    } while(sl_Cd || sl_C);

    update_root_to_subspace();

    return std::make_pair(Cd_elements, C_elements);
  }

  // Hilbert space dimension
  sv_index_type dim() const { return ds.size(); }

  // Number of subspaces
  sv_index_type n_subspaces() const { return ds.n_sets(); }

  // Find what invariant subspace a given basis state belongs to
  sv_index_type operator[](sv_index_type index) const {
    auto it = root_to_subspace.find(ds.find_root(index));
    if(it == root_to_subspace.end())
      throw std::out_of_range("Unexpected basis state index " +
                              std::to_string(index));
    else
      return it->second;
  }

  // Find all subspace-to-subspace connections generated by a given operator
  template <typename HSType, typename LOpScalarType, int... LOpAlgebraIDs>
  connections_map
  find_connections(loperator<LOpScalarType, LOpAlgebraIDs...> const& op,
                   HSType const& hs) const {
    using scalar_type =
        typename loperator<LOpScalarType, LOpAlgebraIDs...>::scalar_type;

    connections_map connections;

    sv_index_type d = get_dim(hs);

    sparse_state_vector<scalar_type> in_state(d);
    sparse_state_vector<scalar_type> out_state(d);
    foreach(hs, [&](sv_index_type in_index) {
      sv_index_type in_subspace =
          root_to_subspace.find(ds.find_root(in_index))->second;

      in_state[in_index] = scalar_traits<scalar_type>::make_const(1);
      op(in_state, out_state);

      foreach(out_state, [&](sv_index_type out_index, scalar_type const&) {
        sv_index_type out_subspace =
            root_to_subspace.find(ds.find_root(out_index))->second;
        connections.emplace(in_subspace, out_subspace);
      });

      set_zeros(in_state);
    });

    return connections;
  }

  // Build a list of all basis states spanning a given subspace 'index'.
  std::vector<sv_index_type> subspace_basis(sv_index_type index) const {
    if(index >= n_subspaces())
      throw std::runtime_error("Wrong subspace index " + std::to_string(index));
    std::vector<sv_index_type> basis;
    foreach(*this, [&](sv_index_type i, sv_index_type subspace) {
      if(subspace == index) basis.emplace_back(i);
    });
    return basis;
  }

  // Build lists of basis states spanning subspaces in this partition.
  std::vector<std::vector<sv_index_type>> subspace_bases() const {
    std::vector<std::vector<sv_index_type>> bases(n_subspaces(),
                                                  std::vector<sv_index_type>{});
    foreach(*this, [&](sv_index_type i, sv_index_type subspace) {
      bases[subspace].emplace_back(i);
    });
    return bases;
  }

  // Apply a functor `f` to all basis states in a given space partition.
  // The functor must take two arguments, index of the basis state,
  // and index of the subspace this basis state belongs to.
  template <typename F> friend void foreach(space_partition const& sp, F&& f) {
    for(sv_index_type n = 0; n < sp.dim(); ++n) {
      std::forward<F>(f)(n, sp[n]);
    }
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
