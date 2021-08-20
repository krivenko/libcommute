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
#ifndef LIBCOMMUTE_LOPERATOR_N_FERMION_SECTOR_VIEW_HPP_
#define LIBCOMMUTE_LOPERATOR_N_FERMION_SECTOR_VIEW_HPP_

#include "../algebra_ids.hpp"
#include "../utility.hpp"

#include "elementary_space_fermion.hpp"
#include "hilbert_space.hpp"
#include "state_vector.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace libcommute {

namespace detail {

// 2^n
inline sv_index_type pow2(unsigned int n) {
  return sv_index_type(1) << n;
}

// Binomial coefficient C(n, k)
inline sv_index_type binomial(unsigned int n, unsigned int k) {
  if(k > n) return 0;
  if(k > n / 2) k = n - k;
  sv_index_type C = 1;
  for(unsigned int i = 0; i < k; ++i)
    C = (C * (n - i)) / (i + 1);
  return C;
}

// Evaluator for sums of the following form,
//
// $$
//   \sum_{j=2}^{\lambda} {{m+1-j}\choose{n-1}} =
//   [m {{m-1}\choose{n-1}} - (m+1-\lambda){{m-\lambda}\choose{n-1}}] / n
// $$
//
// The binomial coefficients are pre-computed and stored upon construction.
struct binomial_sum_t {

  // Pre-computed shifted binomial coefficients stored in a flattened 2D array
  // coeffs[q * (M - N_counted + 1) + p] = Binomial[p + q, q],
  // where $q \in [0, N_counted - 1]$ and $p \in [0, M - N_counted]$.
  std::vector<sv_index_type> coeffs = {};

  // M - N_counted + 1
  unsigned int coeffs_stride;

  inline binomial_sum_t(unsigned int M, unsigned int N_counted)
    : coeffs_stride(M - N_counted + 1) {
    assert(N_counted <= M);
    if(N_counted > 0) {
      // Fill 'coeffs'
      coeffs.resize(N_counted * coeffs_stride);
      std::fill(coeffs.begin(), coeffs.begin() + coeffs_stride, 1);
      for(unsigned int q = 1; q < N_counted; ++q) {
        for(unsigned int p = 0; p < coeffs_stride; ++p) {
          coeffs[q * coeffs_stride + p] =
              ((p + q) * coeffs[(q - 1) * coeffs_stride + p]) / q;
        }
      }
    }
  }

  inline sv_index_type
  operator()(unsigned int n, unsigned int m, unsigned int lambda) const {
    unsigned int lambdam1 = lambda - 1;
    unsigned int coeffs_index_1 = (n - 1) * coeffs_stride + (m - n);
    unsigned int coeffs_index_2 = coeffs_index_1 - lambdam1;
    return (m * coeffs[coeffs_index_1] -
            (m - lambdam1) * coeffs[coeffs_index_2]) /
           n;
  }
};

// Input parameters of an N-fermion sector
struct n_fermion_sector_params_t {
  // Total number of fermionic modes
  unsigned int M;

  // Count occupied or unoccupied fermionic modes
  bool count_occupied;

  // Number of counted modes (either occupied or unoccupied)
  unsigned int N_counted;

  template <typename HSType>
  n_fermion_sector_params_t(HSType const& hs, unsigned int N)
    : M(init_m(hs)),
      count_occupied(N <= M / 2),
      N_counted(count_occupied ? N : M - N) {
    if(N > M)
      throw std::runtime_error("Sector with " + std::to_string(N) +
                               " fermions does not exist "
                               "(there are " +
                               std::to_string(M) +
                               " fermionic modes in total)");
  }

private:
  template <typename HSType> static unsigned int init_m(HSType const& hs) {
    if(hs.total_n_bits() == 0) throw std::runtime_error("Empty Hilbert space");
    if(hs.has_algebra(fermion)) {
      auto fermion_bit_range = hs.algebra_bit_range(fermion);
      // We rely on the following assumption in map_index().
      // The only way it could be violated is by introducing a new elementary
      // space type with an algebra ID < min_user_defined_algebra_id, which is
      // an API misuse anyway.
      assert(fermion_bit_range.first == 0);
      return fermion_bit_range.second - fermion_bit_range.first + 1;
    } else
      return 0;
  }
};

// This generator object uses a non-recursive procedure to enumerate all basis
// state indices in a (pure) N-fermion sector. It does heavy lifting for
// foreach() and n_fermion_sector_basis_states().
struct n_fermion_sector_basis_generator {

  // Total number of fermionic modes
  n_fermion_sector_params_t const& params;

  // Restricted partition of M
  std::vector<unsigned int> lambdas;

  // Upper bounds for elements of 'lambdas'
  std::vector<unsigned int> lambdas_max;

  // XOR-mask used when counting unoccupied states
  const sv_index_type index_mask;

  inline explicit n_fermion_sector_basis_generator(
      n_fermion_sector_params_t const& params)
    : params(params),
      lambdas(params.N_counted, 0),
      lambdas_max(params.N_counted, params.M - (params.N_counted - 1)),
      index_mask(params.count_occupied ? sv_index_type(0) :
                                         (pow2(params.M) - 1)) {}

  // Apply 'f' to each sector basis state
  template <typename F> void operator()(F&& f) {
    if(params.N_counted == 0) {
      f(sv_index_type(0) ^ index_mask);
      return;
    }

    std::fill(lambdas.begin(), lambdas.end(), 0);
    std::fill(lambdas_max.begin(),
              lambdas_max.end(),
              params.M - (params.N_counted - 1));
    for(int i = 0; i >= 0;) {
      ++lambdas[i];
      if(lambdas[i] > lambdas_max[i]) {
        --i;
        continue;
      }

      if(static_cast<unsigned int>(i) + 1 == lambdas.size()) {
        sv_index_type index = 0;
        std::for_each(lambdas.rbegin(),
                      lambdas.rend(),
                      [&index](unsigned int lambda) {
                        index <<= 1;
                        index += 1;
                        index <<= lambda - 1;
                      });
        f(index ^ index_mask);
      } else {
        ++i;
        lambdas[i] = 0;
        lambdas_max[i] = lambdas_max[i - 1] + 1 - lambdas[i - 1];
      }
    }
  }
};

} // namespace detail

//
// N-fermion sector
//

// Size of the fermionic sector with N particles
template <typename HSType>
inline sv_index_type n_fermion_sector_size(HSType const& hs, unsigned int N) {
  auto total_n_bits = hs.total_n_bits();
  if(total_n_bits == 0) return 0;
  unsigned int M = 0;
  if(hs.has_algebra(fermion)) {
    auto fermion_bit_range = hs.algebra_bit_range(fermion);
    M = fermion_bit_range.second - fermion_bit_range.first + 1;
  }
  return detail::binomial(M, N) * detail::pow2(total_n_bits - M);
}

template <typename StateVector, bool Ref = true>
struct n_fermion_sector_view : public detail::n_fermion_sector_params_t {

  // The underlying state vector
  typename std::conditional<Ref, StateVector&, StateVector>::type state_vector;

  // Sum of binomial coefficients
  detail::binomial_sum_t binomial_sum;

  // Number of bits corresponding to the non-fermionic modes
  unsigned int M_nonfermion;

  // Although this object is not used by this class' methods, storing it here
  // allows to eliminate memory allocations in foreach().
  mutable detail::n_fermion_sector_basis_generator basis_generator;

  using scalar_type = typename element_type<
      typename std::remove_const<StateVector>::type>::type;

  template <typename SV, typename HSType>
  n_fermion_sector_view(SV&& sv, HSType const& hs, unsigned int N)
    : detail::n_fermion_sector_params_t(hs, N),
      state_vector(std::forward<SV>(sv)),
      binomial_sum(M, N_counted),
      M_nonfermion(hs.total_n_bits() - M),
      basis_generator(*this) {}

  sv_index_type map_index(sv_index_type index) const {
    unsigned int m = M;
    unsigned int n = N_counted;
    sv_index_type index_f = 0;
    unsigned int lambda = 1;
    // Translate the fermionic part of 'index' into the sector index
    while(n > 0) {
      if((index & sv_index_type(1)) == count_occupied) {
        index_f += binomial_sum(n, m, lambda);
        m -= lambda;
        --n;
        lambda = 1;
      } else {
        ++lambda;
      }
      index >>= 1;
    }

    // Remove the remaining fermionic bits from 'index'
    index >>= m;

    // Here, 'index' contains only non-fermionic bits
    return (index_f << M_nonfermion) + index;
  }
};

// Get element type of the StateVector object adapted by a given
// n_fermion_sector_view object.
template <typename StateVector, bool Ref>
struct element_type<n_fermion_sector_view<StateVector, Ref>> {
  using type = typename n_fermion_sector_view<StateVector>::scalar_type;
};

// Get state amplitude of the adapted StateVector object
// at index view.map_index(n)
template <typename StateVector, bool Ref>
inline auto get_element(n_fermion_sector_view<StateVector, Ref> const& view,
                        sv_index_type n) ->
    typename n_fermion_sector_view<StateVector>::scalar_type {
  return get_element(view.state_vector, view.map_index(n));
}

// Add a constant to a state amplitude stored in the adapted StateVector object
// at index view.map_index(n)
template <typename StateVector, bool Ref, typename T>
inline void update_add_element(n_fermion_sector_view<StateVector, Ref>& view,
                               sv_index_type n,
                               T&& value) {
  update_add_element(view.state_vector,
                     view.map_index(n),
                     std::forward<T>(value));
}

// update_add_element() is not defined for constant views
template <typename StateVector, bool Ref, typename T>
inline void update_add_element(n_fermion_sector_view<StateVector const, Ref>&,
                               sv_index_type,
                               T&&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "update_add_element() is not supported for constant views");
}

// zeros_like() is not defined for views
template <typename StateVector, bool Ref>
inline StateVector zeros_like(n_fermion_sector_view<StateVector, Ref> const&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "zeros_like() is not supported for views");
}

// Set all amplitudes stored in the adapted StateVector object to zero
template <typename StateVector, bool Ref>
inline void set_zeros(n_fermion_sector_view<StateVector, Ref>& view) {
  set_zeros(view.state_vector);
}

// set_zeros() is not defined for constant views
template <typename StateVector, bool Ref>
inline void set_zeros(n_fermion_sector_view<StateVector const, Ref>&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "set_zeros() is not supported for constant views");
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in the adapted StateVector object
template <typename StateVector, bool Ref, typename Functor>
inline void foreach(n_fermion_sector_view<StateVector, Ref> const& view,
                    Functor&& f) {
  auto dim_nonfermion = detail::pow2(view.M_nonfermion);
  sv_index_type sector_index = 0;
  view.basis_generator([&](sv_index_type index_f) {
    for(sv_index_type index_nf = 0; index_nf < dim_nonfermion; ++index_nf) {

      // Emulate decltype(auto)
      decltype(get_element(view.state_vector, sector_index)) a =
          get_element(view.state_vector, sector_index);

      ++sector_index;

      using T = typename n_fermion_sector_view<StateVector, Ref>::scalar_type;
      if(scalar_traits<T>::is_zero(a))
        continue;
      else
        f(index_f + (index_nf << view.M), a);
    }
  });
}

// Make a list of basis state indices spanning the N-fermion sector
// The order of indices is consistent with that used by n_fermion_sector_view
template <typename HSType>
inline std::vector<sv_index_type>
n_fermion_sector_basis_states(HSType const& hs, unsigned int N) {
  detail::n_fermion_sector_params_t params(hs, N);
  detail::n_fermion_sector_basis_generator gen(params);

  unsigned int M_nonfermion = hs.total_n_bits() - params.M;
  sv_index_type dim_nonfermion = detail::pow2(M_nonfermion);

  std::vector<sv_index_type> basis_states;
  basis_states.reserve(detail::binomial(params.M, N) * dim_nonfermion);
  gen([&](sv_index_type index_f) {
    for(sv_index_type index_nf = 0; index_nf < dim_nonfermion; ++index_nf) {
      basis_states.push_back(index_f + (index_nf << params.M));
    }
  });

  return basis_states;
}

template <typename StateVector>
using make_nfs_view_ret_t =
    n_fermion_sector_view<remove_cvref_t<StateVector>,
                          std::is_lvalue_reference<StateVector>::value>;

// Make a non-constant N-fermion sector view
template <typename StateVector, typename HSType>
auto make_nfs_view(StateVector&& sv, HSType const& hs, unsigned int N)
    -> make_nfs_view_ret_t<StateVector> {
  return make_nfs_view_ret_t<StateVector>(std::forward<StateVector>(sv), hs, N);
}

template <typename StateVector>
using make_const_nfs_view_ret_t =
    n_fermion_sector_view<remove_cvref_t<StateVector> const,
                          std::is_lvalue_reference<StateVector>::value>;

// Make a constant N-fermion sector view
template <typename StateVector, typename HSType>
auto make_const_nfs_view(StateVector&& sv, HSType const& hs, unsigned int N)
    -> make_const_nfs_view_ret_t<StateVector> {
  return make_const_nfs_view_ret_t<StateVector>(std::forward<StateVector>(sv),
                                                hs,
                                                N);
}

//
// N-fermion multisector
//

// Description of a N-particle sector defined over a subset of fermionic
// degrees of freedom
template <typename HSType> struct sector_descriptor {

  // Indices defining this sector
  std::set<typename HSType::index_types> indices;

  // Number of particles in the sector
  unsigned int N;
};

namespace detail {

template <typename... IndexTypes>
void validate_sectors(
    hilbert_space<IndexTypes...> const& hs,
    std::vector<sector_descriptor<hilbert_space<IndexTypes...>>> const&
        sectors) {

  unsigned int M_total = 0;
  std::set<typename hilbert_space<IndexTypes...>::index_types> all_indices;
  for(auto const& sector : sectors) {
    for(auto const& ind : sector.indices) {
      if(!hs.has(elementary_space_fermion<IndexTypes...>(ind))) {
        std::stringstream ss;
        print_tuple(ss, ind);
        throw std::runtime_error("Fermionic elementary space with indices " +
                                 ss.str() +
                                 " is not part of this Hilbert space");
      }
      all_indices.insert(ind);
    }
    M_total += sector.indices.size();
  }
  if(all_indices.size() != M_total)
    throw std::runtime_error("Some of the sectors overlap");
}

} // namespace detail

// Size of a fermionic multisector
template <typename HSType>
inline sv_index_type n_fermion_multisector_size(
    HSType const& hs,
    std::vector<sector_descriptor<HSType>> const& sectors) {

  detail::validate_sectors(hs, sectors);

  auto total_n_bits = hs.total_n_bits();
  if(total_n_bits == 0) return 0;

  unsigned int M_total = 0;
  sv_index_type multisector_size = 1;
  for(auto const& sector : sectors) {
    unsigned int M = sector.indices.size();
    if(sector.N > M) return 0;
    multisector_size *= detail::binomial(M, sector.N);
    M_total += M;
  }

  return multisector_size * detail::pow2(total_n_bits - M_total);
}

} // namespace libcommute

#endif
