/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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
#include <functional>
#include <iterator>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace libcommute {

template <typename HSType> struct sector_descriptor;

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
    coeffs.shrink_to_fit();
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
    validate_n(N);
  }

  template <typename HSType>
  n_fermion_sector_params_t(HSType const& hs,
                            sector_descriptor<HSType> const& sector)
    : M(sector.indices.size()),
      count_occupied(sector.N <= M / 2),
      N_counted(count_occupied ? sector.N : M - sector.N) {
    validate_n(sector.N);
  }

private:
  template <typename HSType> static unsigned int init_m(HSType const& hs) {
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

  void validate_n(unsigned int N) const {
    if(N > M)
      throw std::runtime_error("Sector with " + std::to_string(N) +
                               " fermions does not exist "
                               "(there are " +
                               std::to_string(M) + " fermionic modes in it)");
  }
};

// This object uses a non-recursive procedure to generate all restricted N-part
// compositions of an integer.
// It does heavy lifting for foreach() and n_fermion_sector_basis_states().
struct for_each_composition {

  // Parameters of the sector
  n_fermion_sector_params_t const& params;

  // Restricted partition of M
  std::vector<unsigned int> mutable lambdas;

  // Upper bounds for elements of 'lambdas'
  std::vector<unsigned int> mutable lambdas_max;

  inline explicit for_each_composition(n_fermion_sector_params_t const& params)
    : params(params),
      lambdas(params.N_counted, 0),
      lambdas_max(params.N_counted, params.M - params.N_counted + 1) {
    lambdas.shrink_to_fit();
    lambdas_max.shrink_to_fit();
  }

  // Apply 'f' to each composition
  template <typename F> void operator()(F&& f) const {
    if(params.N_counted == 0) {
      f(lambdas);
      return;
    }

    std::fill(lambdas.begin(), lambdas.end(), 0);
    std::fill(lambdas_max.begin(),
              lambdas_max.end(),
              params.M - params.N_counted + 1);
    for(int i = 0; i >= 0;) {
      ++lambdas[i];
      if(lambdas[i] > lambdas_max[i]) {
        --i;
        continue;
      }

      if(static_cast<unsigned int>(i + 1) == lambdas.size()) {
        f(lambdas);
      } else {
        ++i;
        lambdas[i] = 0;
        lambdas_max[i] = lambdas_max[i - 1] + 1 - lambdas[i - 1];
      }
    }
  }
};

// Combine a composition {\lambda_i} and a basis state index from
// the non-fermionic subspace into the index in the full Hilbert space
struct composition_to_full_hs_index {

  // Total number of fermionic modes
  unsigned int M;

  // XOR-mask used when counting unoccupied states
  sv_index_type const index_mask;

  explicit composition_to_full_hs_index(n_fermion_sector_params_t const& params)
    : M(params.M),
      index_mask(params.count_occupied ? sv_index_type(0) : (pow2(M) - 1)) {}

  inline sv_index_type operator()(std::vector<unsigned int> const& lambdas,
                                  sv_index_type index_nf) const {
    sv_index_type index_f = 0;
    if(lambdas.size() != 0) {
      std::for_each(lambdas.rbegin(),
                    lambdas.rend(),
                    [&index_f](unsigned int lambda) {
                      index_f <<= 1;
                      index_f += 1;
                      index_f <<= lambda - 1;
                    });
    }
    return (index_f ^ index_mask) + (index_nf << M);
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

  // Although the following two objects are not used by this class' methods,
  // storing them here allows to eliminate memory allocations in foreach().
  detail::for_each_composition for_each_comp;
  detail::composition_to_full_hs_index comp_to_index;

  using scalar_type = typename element_type<
      typename std::remove_const<StateVector>::type>::type;

  template <typename SV, typename HSType>
  n_fermion_sector_view(SV&& sv, HSType const& hs, unsigned int N)
    : detail::n_fermion_sector_params_t(hs, N),
      state_vector(std::forward<SV>(sv)),
      binomial_sum(M, N_counted),
      M_nonfermion(hs.total_n_bits() - M),
      for_each_comp(*this),
      comp_to_index(*this) {}

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
  if(view.M == 0 && view.M_nonfermion == 0) return;

  auto dim_nonfermion = detail::pow2(view.M_nonfermion);
  sv_index_type sector_index = 0;
  view.for_each_comp([&](std::vector<unsigned int> const& lambdas) {
    for(sv_index_type index_nf = 0; index_nf < dim_nonfermion; ++index_nf) {

      // Emulate decltype(auto)
      decltype(get_element(view.state_vector, sector_index)) a =
          get_element(view.state_vector, sector_index);

      ++sector_index;

      using T = typename n_fermion_sector_view<StateVector, Ref>::scalar_type;
      if(scalar_traits<T>::is_zero(a))
        continue;
      else
        f(view.comp_to_index(lambdas, index_nf), a);
    }
  });
}

// Make a list of basis state indices spanning the N-fermion sector
// The order of indices is consistent with that used by n_fermion_sector_view
template <typename HSType>
inline std::vector<sv_index_type>
n_fermion_sector_basis_states(HSType const& hs, unsigned int N) {
  detail::n_fermion_sector_params_t params(hs, N);

  unsigned int M_nonfermion = hs.total_n_bits() - params.M;
  if(params.M == 0 && M_nonfermion == 0) return {};

  detail::for_each_composition for_each_comp(params);
  detail::composition_to_full_hs_index comp_to_index(params);

  sv_index_type dim_nonfermion = detail::pow2(M_nonfermion);

  std::vector<sv_index_type> basis_states;
  basis_states.reserve(detail::binomial(params.M, N) * dim_nonfermion);
  for_each_comp([&](std::vector<unsigned int> const& lambdas) {
    for(sv_index_type index_nf = 0; index_nf < dim_nonfermion; ++index_nf) {
      basis_states.push_back(comp_to_index(lambdas, index_nf));
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

using sector_params_vec_t = std::vector<n_fermion_sector_params_t>;

template <typename HSType>
sector_params_vec_t
make_sector_params(HSType const& hs,
                   std::vector<sector_descriptor<HSType>> const& sectors) {

  validate_sectors(hs, sectors);

  sector_params_vec_t params;
  for(auto const& sector : sectors) {
    n_fermion_sector_params_t p(hs, sector);
    // Initialize only non-empty sectors
    if(p.M != 0) params.emplace_back(p);
  }
  params.shrink_to_fit();
  return params;
}

template <typename... IndexTypes>
unsigned int get_fermion_bit(
    hilbert_space<IndexTypes...> const& hs,
    typename hilbert_space<IndexTypes...>::index_types const& indices) {
  return hs.bit_range(elementary_space_fermion<IndexTypes...>(indices)).first;
}

// This object uses a non-recursive procedure to generate all restricted N-part
// compositions of multiple integers in parallel.
// It does heavy lifting for foreach() and n_fermion_multisector_basis_states().
struct for_each_composition_multi {

  // Parameters of sectors
  sector_params_vec_t const& sector_params;

  // Index of the first sector with N_computed != 0
  int s_min = -1;

  // Index of the last sector with N_computed != 0
  int s_max = -1;

  // Restricted partitions of M_s
  std::vector<std::vector<unsigned int>> mutable lambdas;

  // Upper bounds for elements of 'lambdas'
  std::vector<std::vector<unsigned int>> mutable lambdas_max;

  inline explicit for_each_composition_multi(
      sector_params_vec_t const& sector_params)
    : sector_params(sector_params) {
    for(int s = 0; s < static_cast<int>(sector_params.size()); ++s) {
      auto const& p = sector_params[s];
      lambdas.emplace_back(p.N_counted, 0);
      lambdas.back().shrink_to_fit();
      lambdas_max.emplace_back(p.N_counted, p.M - (p.N_counted - 1));
      lambdas_max.back().shrink_to_fit();
      if(p.N_counted != 0) {
        if(s_min == -1) s_min = s;
        s_max = s;
      }
    }
    lambdas.shrink_to_fit();
    lambdas_max.shrink_to_fit();
  }

  // Apply 'f' to each composition
  template <typename F> void operator()(F&& f) const {
    if(s_min == -1) { // All compositions are empty
      f(lambdas);
      return;
    }

    for(std::size_t s = 0; s < lambdas.size(); ++s) {
      std::fill(lambdas[s].begin(), lambdas[s].end(), 0);
      std::fill(lambdas_max[s].begin(),
                lambdas_max[s].end(),
                sector_params[s].M - sector_params[s].N_counted + 1);
    }

    int s = s_min; // Index of current sector
    int i = 0;     // Index within the current sector

    while(s >= s_min) {
      ++lambdas[s][i];
      if(lambdas[s][i] > lambdas_max[s][i]) {
        --i;
        if(i < 0) {
          do {
            --s;
          } while(s >= s_min && sector_params[s].N_counted == 0);
          if(s >= s_min)
            i = static_cast<int>(sector_params[s].N_counted) - 1;
        }
        continue;
      }

      if(s == s_max &&
         static_cast<unsigned int>(i + 1) == sector_params[s].N_counted) {
        f(lambdas);
      } else {
        ++i;
        if(i == static_cast<int>(sector_params[s].N_counted)) {
          do {
            ++s;
          } while(sector_params[s].N_counted == 0);
          auto const& p = sector_params[s];
          i = 0;
          lambdas_max[s][i] = p.M - p.N_counted + 1;
        } else {
          lambdas_max[s][i] = lambdas_max[s][i - 1] + 1 - lambdas[s][i - 1];
        }
        lambdas[s][i] = 0;
      }
    }
  }
};

// 'non_multisector_bit' stands for a bit corresponding to a non-multisector
// fermionic degree of freedom.
static constexpr int non_multisector_bit = -1;

// Combine a composition {\lambda_i} and a basis state index from
// the non-multisector subspace into the index in the full Hilbert space
struct compositions_to_full_hs_index {

  // Parameters of sectors
  sector_params_vec_t const& sector_params;

  // Mapping from bit position to sector index
  std::vector<int> const& bit_to_sector;

  // XOR-mask used to flip bits in all count_occupied == false sectors
  sv_index_type const index_mask;

  // Index withing each sector
  std::vector<sv_index_type> mutable sector_indices;

  compositions_to_full_hs_index(sector_params_vec_t const& sector_params,
                                std::vector<int> const& bit_to_sector)
    : sector_params(sector_params),
      bit_to_sector(bit_to_sector),
      index_mask(init_index_mask(sector_params, bit_to_sector)),
      sector_indices(sector_params.size()) {
    sector_indices.shrink_to_fit();
  }

  inline sv_index_type
  operator()(std::vector<std::vector<unsigned int>> const& lambdas,
             sv_index_type index_nf) const {

    for(unsigned int s = 0; s < lambdas.size(); ++s) {
      if(lambdas[s].size() != 0) {
        auto& index_f = sector_indices[s];
        std::for_each(lambdas[s].rbegin(),
                      lambdas[s].rend(),
                      [&index_f](unsigned int lambda) {
                        index_f <<= 1;
                        index_f += 1;
                        index_f <<= lambda - 1;
                      });
      } else
        sector_indices[s] = 0;
    }

    sv_index_type index = 0;

    // Consume bits from sector_indices and index_nf
    for(unsigned int b = 0; b < bit_to_sector.size(); ++b) {
      int s = bit_to_sector[b];
      if(s != non_multisector_bit) {
        auto& index_f = sector_indices[s];
        index += (index_f & sv_index_type(1)) << b;
        index_f >>= 1;
      } else {
        index += (index_nf & sv_index_type(1)) << b;
        index_nf >>= 1;
      }
    }

    index += index_nf << bit_to_sector.size();

    return index ^ index_mask;
  }

private:
  static sv_index_type init_index_mask(sector_params_vec_t const& sector_params,
                                       std::vector<int> const& bit_to_sector) {
    sv_index_type mask = 0;
    for(unsigned int b = 0; b < bit_to_sector.size(); ++b) {
      int s = bit_to_sector[b];
      if(s != non_multisector_bit && (!sector_params[s].count_occupied)) {
        mask += pow2(b);
      }
    }
    return mask;
  }
};

template <typename HSType>
static std::vector<int>
make_bit_to_sector(HSType const& hs,
                   std::vector<sector_descriptor<HSType>> const& sectors) {

  if(!hs.has_algebra(fermion)) return {};

  std::vector<int> result(hs.algebra_bit_range(fermion).second + 1,
                          non_multisector_bit);
  unsigned int result_size = 0;

  int sector_i = 0;
  for(auto const& sector : sectors) {
    // Process only non-empty sectors
    if(!sector.indices.empty()) {
      for(auto const& ind : sector.indices) {
        unsigned int b = detail::get_fermion_bit(hs, ind);
        result[b] = sector_i;
        result_size = std::max(result_size, b + 1);
      }
      ++sector_i;
    }
  }

  result.resize(result_size);
  result.shrink_to_fit();
  return result;
}

template <typename HSType>
unsigned int
compute_m_nonmultisector(HSType const& hs,
                         sector_params_vec_t const& sector_params) {
  return hs.total_n_bits() -
         std::accumulate(
             sector_params.begin(),
             sector_params.end(),
             0,
             [](unsigned int M, n_fermion_sector_params_t const& data) {
               return M + data.M;
             });
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

template <typename StateVector, bool Ref = true>
struct n_fermion_multisector_view {

  // The underlying state vector
  typename std::conditional<Ref, StateVector&, StateVector>::type state_vector;

  // Parameters of sectors
  detail::sector_params_vec_t sector_params;

  // Various per-sector data used by map_index()
  struct sector_map_index_data_t {

    // Sum of binomial coefficients
    detail::binomial_sum_t binomial_sum;

    // Current value of m
    unsigned int m = 0;

    // Current value of n
    unsigned int n = 0;

    // Current value of \lambda
    unsigned int lambda = 1;

    // Current index within this sector
    sv_index_type index = 0;

    explicit sector_map_index_data_t(
        detail::n_fermion_sector_params_t const& params)
      : binomial_sum(params.M, params.N_counted) {}

    // Reset all current values
    void reset(detail::n_fermion_sector_params_t const& params) {
      m = params.M;
      n = params.N_counted;
      index = 0;
      lambda = 1;
    }
  };
  std::vector<sector_map_index_data_t> mutable sector_map_index_data;

  // Mapping from bit position to sector index.
  std::vector<int> bit_to_sector;

  // Strides used to combine sector indices into a multisector index
  std::vector<sv_index_type> sector_strides;

  // The number of bits not belonging to the multisector
  unsigned int M_nonmultisector;

  // Although the following two objects are not used by this class' methods,
  // storing them here allows to eliminate memory allocations in foreach().
  detail::for_each_composition_multi for_each_comp;
  detail::compositions_to_full_hs_index comp_to_index;

  using scalar_type = typename element_type<
      typename std::remove_const<StateVector>::type>::type;

  template <typename SV, typename HSType>
  n_fermion_multisector_view(
      SV&& sv,
      HSType const& hs,
      std::vector<sector_descriptor<HSType>> const& sectors)
    : state_vector(std::forward<SV>(sv)),
      sector_params(detail::make_sector_params(hs, sectors)),
      bit_to_sector(detail::make_bit_to_sector(hs, sectors)),
      sector_strides(init_sector_strides(sector_params)),
      M_nonmultisector(detail::compute_m_nonmultisector(hs, sector_params)),
      for_each_comp(sector_params),
      comp_to_index(sector_params, bit_to_sector) {
    sector_map_index_data.reserve(sector_params.size());
    for(auto const& p : sector_params) {
      // cppcheck-suppress useStlAlgorithm
      sector_map_index_data.emplace_back(p);
    }
  }

  sv_index_type map_index(sv_index_type index) const {

    // Translate the fermionic part of 'index' into a list of sector indices
    for(std::size_t s = 0; s < sector_params.size(); ++s)
      sector_map_index_data[s].reset(sector_params[s]);

    sv_index_type index_nonmultisector = 0;
    unsigned int b_nonmultisector = 0;

    for(int s : bit_to_sector) {
      if(s !=
         detail::non_multisector_bit) { // This bit belongs to the multisector
        auto& d = sector_map_index_data[s];
        if((index & sv_index_type(1)) == sector_params[s].count_occupied) {
          d.index += d.binomial_sum(d.n, d.m, d.lambda);
          d.m -= d.lambda;
          --d.n;
          d.lambda = 1;
        } else {
          ++d.lambda;
        }
      } else { // a non-multisector bit
        index_nonmultisector += (index & sv_index_type(1)) << b_nonmultisector;
        ++b_nonmultisector;
      }
      index >>= 1;
    }

    index_nonmultisector += index << b_nonmultisector;

    auto index_multisector = std::inner_product(
        sector_map_index_data.begin(),
        sector_map_index_data.end(),
        sector_strides.begin(),
        0,
        std::plus<sv_index_type>(),
        [](sector_map_index_data_t const& data, sv_index_type stride) {
          return data.index * stride;
        });
    return (index_multisector << M_nonmultisector) + index_nonmultisector;
  }

private:
  static std::vector<sv_index_type>
  init_sector_strides(detail::sector_params_vec_t const& sector_params) {
    std::vector<sv_index_type> strides(sector_params.size());
    strides.shrink_to_fit();
    if(strides.empty()) return strides;

    strides[strides.size() - 1] = 1;
    for(std::size_t i = strides.size() - 1; i-- != 0;)
      strides[i] =
          strides[i + 1] * detail::binomial(sector_params[i + 1].M,
                                            sector_params[i + 1].N_counted);

    return strides;
  }
};

// Get element type of the StateVector object adapted by a given
// n_fermion_multisector_view object.
template <typename StateVector, bool Ref>
struct element_type<n_fermion_multisector_view<StateVector, Ref>> {
  using type = typename n_fermion_multisector_view<StateVector>::scalar_type;
};

// Get state amplitude of the adapted StateVector object
// at index view.map_index(n)
template <typename StateVector, bool Ref>
inline auto
get_element(n_fermion_multisector_view<StateVector, Ref> const& view,
            sv_index_type n) ->
    typename n_fermion_sector_view<StateVector>::scalar_type {
  return get_element(view.state_vector, view.map_index(n));
}

// Add a constant to a state amplitude stored in the adapted StateVector object
// at index view.map_index(n)
template <typename StateVector, bool Ref, typename T>
inline void
update_add_element(n_fermion_multisector_view<StateVector, Ref>& view,
                   sv_index_type n,
                   T&& value) {
  update_add_element(view.state_vector,
                     view.map_index(n),
                     std::forward<T>(value));
}

// update_add_element() is not defined for constant views
template <typename StateVector, bool Ref, typename T>
inline void
update_add_element(n_fermion_multisector_view<StateVector const, Ref>&,
                   sv_index_type,
                   T&&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "update_add_element() is not supported for constant views");
}

// zeros_like() is not defined for views
template <typename StateVector, bool Ref>
inline StateVector
zeros_like(n_fermion_multisector_view<StateVector, Ref> const&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "zeros_like() is not supported for views");
}

// Set all amplitudes stored in the adapted StateVector object to zero
template <typename StateVector, bool Ref>
inline void set_zeros(n_fermion_multisector_view<StateVector, Ref>& view) {
  set_zeros(view.state_vector);
}

// set_zeros() is not defined for constant views
template <typename StateVector, bool Ref>
inline void set_zeros(n_fermion_multisector_view<StateVector const, Ref>&) {
  static_assert(!std::is_same<StateVector, StateVector>::value,
                "set_zeros() is not supported for constant views");
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in the adapted StateVector object
template <typename StateVector, bool Ref, typename Functor>
inline void foreach(n_fermion_multisector_view<StateVector, Ref> const& view,
                    Functor&& f) {
  if(view.sector_params.size() == 0 && view.M_nonmultisector == 0) return;

  auto dim_nonmultisector = detail::pow2(view.M_nonmultisector);
  sv_index_type multisector_index = 0;
  view.for_each_comp([&](std::vector<std::vector<unsigned int>> const&
                             lambdas) {
    for(sv_index_type index_nms = 0; index_nms < dim_nonmultisector;
        ++index_nms) {

      // Emulate decltype(auto)
      decltype(get_element(view.state_vector, multisector_index)) a =
          get_element(view.state_vector, multisector_index);

      ++multisector_index;

      using T =
          typename n_fermion_multisector_view<StateVector, Ref>::scalar_type;
      if(scalar_traits<T>::is_zero(a))
        continue;
      else
        f(view.comp_to_index(lambdas, index_nms), a);
    }
  });
}

// Make a list of basis state indices spanning the N-fermion multisector
// The order of indices is consistent with that used by
// n_fermion_multisector_view
template <typename HSType>
inline std::vector<sv_index_type> n_fermion_multisector_basis_states(
    HSType const& hs,
    std::vector<sector_descriptor<HSType>> const& sectors) {
  auto sector_params = detail::make_sector_params(hs, sectors);
  if(hs.total_n_bits() == 0 && sector_params.empty()) return {};

  detail::for_each_composition_multi for_each_comp(sector_params);
  auto bit_to_sector = detail::make_bit_to_sector(hs, sectors);
  detail::compositions_to_full_hs_index comp_to_index(sector_params,
                                                      bit_to_sector);

  sv_index_type dim_multisector = std::accumulate(
      sector_params.begin(),
      sector_params.end(),
      1,
      [](sv_index_type dim, detail::n_fermion_sector_params_t const& p) {
        return dim * detail::binomial(p.M, p.N_counted);
      });

  unsigned int M_nonmultisector = compute_m_nonmultisector(hs, sector_params);
  sv_index_type dim_nonmultisector = detail::pow2(M_nonmultisector);

  std::vector<sv_index_type> basis_states;
  basis_states.reserve(dim_multisector * dim_nonmultisector);
  for_each_comp([&](std::vector<std::vector<unsigned int>> const& lambdas) {
    for(sv_index_type index_nms = 0; index_nms < dim_nonmultisector;
        ++index_nms) {
      basis_states.push_back(comp_to_index(lambdas, index_nms));
    }
  });

  return basis_states;
}

template <typename StateVector>
using make_nfms_view_ret_t =
    n_fermion_multisector_view<remove_cvref_t<StateVector>,
                               std::is_lvalue_reference<StateVector>::value>;

// Make a non-constant N-fermion sector view
template <typename StateVector, typename HSType>
auto make_nfms_view(StateVector&& sv,
                    HSType const& hs,
                    std::vector<sector_descriptor<HSType>> const& sectors)
    -> make_nfms_view_ret_t<StateVector> {
  return make_nfms_view_ret_t<StateVector>(std::forward<StateVector>(sv),
                                           hs,
                                           sectors);
}

template <typename StateVector>
using make_const_nfms_view_ret_t =
    n_fermion_multisector_view<remove_cvref_t<StateVector> const,
                               std::is_lvalue_reference<StateVector>::value>;

// Make a constant N-fermion sector view
template <typename StateVector, typename HSType>
auto make_const_nfms_view(StateVector&& sv,
                          HSType const& hs,
                          std::vector<sector_descriptor<HSType>> const& sectors)
    -> make_const_nfms_view_ret_t<StateVector> {
  return make_const_nfms_view_ret_t<StateVector>(std::forward<StateVector>(sv),
                                                 hs,
                                                 sectors);
}

} // namespace libcommute

#endif
