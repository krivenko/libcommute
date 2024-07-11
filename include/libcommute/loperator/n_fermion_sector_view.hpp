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

#if defined(__GNUC__) || defined(__clang__)
#include <immintrin.h>
#endif

/* Ranking / unranking algorithms implemented here are described in
 *
 * `"Trie-based ranking of quantum many-body states",
 * M. Wallerberger and K. Held,
 * Phys. Rev. Research 4, 033238 (2022)
 * <https://doi.org/10.1103/PhysRevResearch.4.033238>`_.
 */

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

// Parallel bits deposit
inline sv_index_type deposit_bits(sv_index_type src, sv_index_type mask) {
#if defined(__GNUC__) || defined(__clang__)
  return _pdep_u64(src, mask);
#else
  // https://www.chessprogramming.org/BMI2#Serial_Implementation
  sv_index_type res = 0;
  for(sv_index_type bb = 1; mask; bb += bb) {
    if(src & bb) res |= mask & -mask;
    mask &= mask - 1;
  }
  return res;
#endif
}

// Parallel bits extract
inline sv_index_type extract_bits(sv_index_type val, sv_index_type mask) {
#if defined(__GNUC__) || defined(__clang__)
  return _pext_u64(val, mask);
#else
  // From https://www.chessprogramming.org/BMI2#Serial_Implementation_2
  sv_index_type res = 0;
  for(sv_index_type bb = 1; mask; bb += bb) {
    if(val & mask & -mask) res |= bb;
    mask &= mask - 1;
  }
  return res;
#endif
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

// Parameters of an N-fermion sector
struct n_fermion_sector_params_t {

  // Number of fermionic modes in the sector
  unsigned int M;

  // Bit mask indicating positions of the modes
  sv_index_type mask;

  // Number of occupied modes
  unsigned int N;

  template <typename HSType>
  n_fermion_sector_params_t(HSType const& hs, unsigned int N)
    : M(init_m(hs)), mask((sv_index_type(1) << M) - 1), N(validated_n(N)) {}

  template <typename HSType>
  n_fermion_sector_params_t(HSType const& hs,
                            sector_descriptor<HSType> const& sector)
    : M(sector.indices.size()),
      mask(init_mask(hs, sector)),
      N(validated_n(sector.N)) {}

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

  template <typename... IndexTypes>
  static unsigned int
  init_mask(hilbert_space<IndexTypes...> const& hs,
            sector_descriptor<hilbert_space<IndexTypes...>> const& sector) {
    sv_index_type mask = 0;
    for(auto const& i : sector.indices) {
      unsigned int b =
          hs.bit_range(elementary_space_fermion<IndexTypes...>(i)).first;
      mask += sv_index_type(1) << b;
    }
    return mask;
  }

  inline unsigned int validated_n(unsigned int N) const {
    if(N > M)
      throw std::runtime_error("Sector with " + std::to_string(N) +
                               " fermions does not exist "
                               "(there are " +
                               std::to_string(M) + " fermionic modes in it)");
    return N;
  }
};

} // namespace detail

//
// Ranking and unranking algorithms
//

// Combination ranking algorithm
class combination_ranking {

  // Parameters of the sector
  detail::n_fermion_sector_params_t const& sector_params;

  // Count occupied or unoccupied fermionic modes
  bool count_occupied;

  // Number of counted modes (either occupied or unoccupied)
  unsigned int N_counted;

  // Sum of binomial coefficients
  detail::binomial_sum_t binomial_sum;

public:
  explicit combination_ranking(detail::n_fermion_sector_params_t const& params)
    : sector_params(params),
      count_occupied(params.N <= params.M / 2),
      N_counted(count_occupied ? params.N : params.M - params.N),
      binomial_sum(params.M, N_counted) {}

  // Rank a fermionic many-body state
  sv_index_type operator()(sv_index_type index) const {
    unsigned int m = sector_params.M;
    unsigned int n = N_counted;
    sv_index_type r = 0;
    unsigned int lambda = 1;
    while(n > 0) {
      if((index & sv_index_type(1)) == count_occupied) {
        r += binomial_sum(n, m, lambda);
        m -= lambda;
        --n;
        lambda = 1;
      } else {
        ++lambda;
      }
      index >>= 1;
    }
    return r;
  }
};

// Generator of unranked basis states
class unranking_generator {

  // Parameters of the sector
  detail::n_fermion_sector_params_t const& sector_params;

  // Number of counted modes (either occupied or unoccupied)
  unsigned int N_counted;

  // Restricted partition of params.M
  std::vector<unsigned int> mutable lambdas;

  // Upper bounds for elements of 'lambdas'
  std::vector<unsigned int> mutable lambdas_max;

  // Index of the element of 'lambdas' to be updated to generate
  // the next state
  int mutable i = 0;

  // XOR-mask used when counting unoccupied states
  sv_index_type const unoccupied_mask;

public:
  explicit unranking_generator(detail::n_fermion_sector_params_t const& params)
    : sector_params(params),
      N_counted(params.N <= params.M / 2 ? params.N : params.M - params.N),
      lambdas(N_counted, 0),
      lambdas_max(N_counted, params.M - N_counted + 1),
      unoccupied_mask(params.N <= params.M / 2 ? sv_index_type(0) :
                                                 (detail::pow2(params.M) - 1)) {
    lambdas.shrink_to_fit();
    lambdas_max.shrink_to_fit();
  }

  // Initialize iteration
  inline void init() const {
    std::fill(lambdas.begin(), lambdas.end(), 0);
    std::fill(lambdas_max.begin(),
              lambdas_max.end(),
              sector_params.M - N_counted + 1);
    i = 0;
  }

  // Return the next state
  inline sv_index_type next() const {
    if(N_counted == 0) {
      --i;
      return unoccupied_mask;
    }

    while(i >= 0) {
      ++lambdas[i];
      if(lambdas[i] > lambdas_max[i]) {
        --i;
        continue;
      }

      if(static_cast<unsigned int>(i + 1) == lambdas.size()) break;

      ++i;
      lambdas[i] = 0;
      lambdas_max[i] = lambdas_max[i - 1] + 1 - lambdas[i - 1];
    }

    sv_index_type index = 0;
    std::for_each(lambdas.rbegin(),
                  lambdas.rend(),
                  [&index](unsigned int lambda) {
                    index <<= 1;
                    index += 1;
                    index <<= lambda - 1;
                  });

    return index ^ unoccupied_mask;
  }

  // Are there still states to be returned by next()?
  inline bool done() const {
    if(N_counted == 0)
      return i != 0;
    else
      return lambdas[0] >= lambdas_max[0];
  }

  // Total number of generated unranked states
  inline sv_index_type size() const {
    return detail::binomial(sector_params.M, sector_params.N);
  }
};

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

template <typename StateVector,
          bool Ref = true,
          typename RankingAlgorithm = combination_ranking>
struct n_fermion_sector_view {

  // The underlying state vector
  typename std::conditional<Ref, StateVector&, StateVector>::type state_vector;

  // Parameters of the sector
  detail::n_fermion_sector_params_t sector_params;

  // State ranking algorithm
  RankingAlgorithm ranking;

  // Generator of unranked states
  unranking_generator unranking;

  // Mask used to select fermionic modes
  sv_index_type f_mask;

  // Number of bits corresponding to the non-fermionic modes
  unsigned int M_nonfermion;

  using scalar_type = typename element_type<
      typename std::remove_const<StateVector>::type>::type;

  template <typename SV, typename HSType>
  n_fermion_sector_view(SV&& sv, HSType const& hs, unsigned int N)
    : state_vector(std::forward<SV>(sv)),
      sector_params(hs, N),
      ranking(sector_params),
      unranking(sector_params),
      f_mask((sv_index_type(1) << sector_params.M) - 1),
      M_nonfermion(hs.total_n_bits() - sector_params.M) {}

  sv_index_type map_index(sv_index_type index) const {
    sv_index_type ranked = ranking(index & f_mask);
    // Place the non-fermionic bits of 'index' to the least significant
    // positions and put the computed rank of the fermionic part after them.
    return (ranked << M_nonfermion) + (index >> sector_params.M);
  }
};

// Get element type of the StateVector object adapted by a given
// n_fermion_sector_view object.
template <typename StateVector, bool Ref, typename RankingAlgorithm>
struct element_type<n_fermion_sector_view<StateVector, Ref, RankingAlgorithm>> {
  using type = typename n_fermion_sector_view<StateVector,
                                              Ref,
                                              RankingAlgorithm>::scalar_type;
};

// Get state amplitude of the adapted StateVector object
// at index view.map_index(n)
template <typename StateVector, bool Ref, typename RankingAlgorithm>
inline auto get_element(
    n_fermion_sector_view<StateVector, Ref, RankingAlgorithm> const& view,
    sv_index_type n) ->
    typename n_fermion_sector_view<StateVector, Ref, RankingAlgorithm>::
        scalar_type {
  return get_element(view.state_vector, view.map_index(n));
}

// Add a constant to a state amplitude stored in the adapted StateVector object
// at index view.map_index(n)
template <typename StateVector, bool Ref, typename RankingAlgorithm, typename T>
inline void update_add_element(
    n_fermion_sector_view<StateVector, Ref, RankingAlgorithm>& view,
    sv_index_type n,
    T&& value) {
  update_add_element(view.state_vector,
                     view.map_index(n),
                     std::forward<T>(value));
}

// update_add_element() is not defined for constant views
template <typename StateVector, bool Ref, typename RankingAlgorithm, typename T>
inline void update_add_element(
    n_fermion_sector_view<StateVector const, Ref, RankingAlgorithm>&,
    sv_index_type,
    T&&) { // NOLINT(cppcoreguidelines-missing-std-forward)
  static_assert(false,
                "update_add_element() is not supported for constant views");
}

// zeros_like() is not defined for views
template <typename StateVector, bool Ref, typename RankingAlgorithm>
inline StateVector
zeros_like(n_fermion_sector_view<StateVector, Ref, RankingAlgorithm> const&) {
  static_assert(false, "zeros_like() is not supported for views");
}

// Set all amplitudes stored in the adapted StateVector object to zero
template <typename StateVector, bool Ref, typename RankingAlgorithm>
inline void
set_zeros(n_fermion_sector_view<StateVector, Ref, RankingAlgorithm>& view) {
  set_zeros(view.state_vector);
}

// set_zeros() is not defined for constant views
template <typename StateVector, bool Ref, typename RankingAlgorithm>
inline void
set_zeros(n_fermion_sector_view<StateVector const, Ref, RankingAlgorithm>&) {
  static_assert(false, "set_zeros() is not supported for constant views");
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in the adapted StateVector object
template <typename StateVector,
          bool Ref,
          typename RankingAlgorithm,
          typename Functor>
inline void
foreach(n_fermion_sector_view<StateVector, Ref, RankingAlgorithm> const& view,
        Functor&& f) { // NOLINT(cppcoreguidelines-missing-std-forward)
  if(view.sector_params.M == 0 && view.M_nonfermion == 0) return;

  auto dim_nonfermion = detail::pow2(view.M_nonfermion);
  sv_index_type sector_index = 0;

  view.unranking.init();
  while(!view.unranking.done()) {
    sv_index_type unranked = view.unranking.next();
    for(sv_index_type index_nf = 0; index_nf < dim_nonfermion; ++index_nf) {

      // Emulate decltype(auto)
      decltype(get_element(view.state_vector, sector_index)) a =
          get_element(view.state_vector, sector_index);

      ++sector_index;

      using T = typename n_fermion_sector_view<StateVector,
                                               Ref,
                                               RankingAlgorithm>::scalar_type;
      if(scalar_traits<T>::is_zero(a))
        continue;
      else
        std::forward<Functor>(f)(unranked + (index_nf << view.sector_params.M),
                                 a);
    }
  }
}

// Make a list of basis state indices spanning the N-fermion sector
// The order of indices is consistent with that used by n_fermion_sector_view.
template <typename HSType, typename RankingAlgorithm = combination_ranking>
inline std::vector<sv_index_type>
n_fermion_sector_basis_states(HSType const& hs, unsigned int N) {
  detail::n_fermion_sector_params_t sector_params(hs, N);

  unsigned int M_nonfermion = hs.total_n_bits() - sector_params.M;
  if(sector_params.M == 0 && M_nonfermion == 0) return {};

  unranking_generator unranking(sector_params);

  sv_index_type dim_nonfermion = detail::pow2(M_nonfermion);

  std::vector<sv_index_type> basis_states;
  basis_states.reserve(detail::binomial(sector_params.M, N) * dim_nonfermion);

  unranking.init();
  while(!unranking.done()) {
    sv_index_type unranked = unranking.next();
    for(sv_index_type index_nf = 0; index_nf < dim_nonfermion; ++index_nf) {
      basis_states.push_back(unranked + (index_nf << sector_params.M));
    }
  }

  return basis_states;
}

template <typename StateVector, typename RankingAlgorithm>
using make_nfs_view_ret_t =
    n_fermion_sector_view<remove_cvref_t<StateVector>,
                          std::is_lvalue_reference<StateVector>::value,
                          RankingAlgorithm>;

// Make a non-constant N-fermion sector view
template <typename StateVector,
          typename HSType,
          typename RankingAlgorithm = combination_ranking>
auto make_nfs_view(StateVector&& sv, HSType const& hs, unsigned int N)
    -> make_nfs_view_ret_t<StateVector, RankingAlgorithm> {
  return make_nfs_view_ret_t<StateVector, RankingAlgorithm>(
      std::forward<StateVector>(sv),
      hs,
      N);
}

template <typename StateVector, typename RankingAlgorithm>
using make_const_nfs_view_ret_t =
    n_fermion_sector_view<remove_cvref_t<StateVector> const,
                          std::is_lvalue_reference<StateVector>::value,
                          RankingAlgorithm>;

// Make a constant N-fermion sector view
template <typename StateVector,
          typename HSType,
          typename RankingAlgorithm = combination_ranking>
auto make_const_nfs_view(StateVector&& sv, HSType const& hs, unsigned int N)
    -> make_const_nfs_view_ret_t<StateVector, RankingAlgorithm> {
  return make_const_nfs_view_ret_t<StateVector, RankingAlgorithm>(
      std::forward<StateVector>(sv),
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

template <typename HSType>
unsigned int make_nonmultisector_mask(HSType const& hs,
                                      sector_params_vec_t const& sector_params,
                                      sv_index_type& mask) {
  unsigned int M = hs.total_n_bits();
  mask = (sv_index_type(1) << hs.total_n_bits()) - 1;
  for(auto const& p : sector_params) {
    mask -= p.mask;
    M -= p.M;
  }
  return M;
}

// Generator of multi-sector unranked basis states
class multisector_unranking_generator {

  // Parameters of sectors
  detail::sector_params_vec_t const& sector_params;

  // State unranking generators for each sector
  std::vector<unranking_generator> sector_gens;

  // Total number of generated unranked states
  sv_index_type total;

  // Index of current unranked state
  sv_index_type mutable i = 0;

  // Index of the sector currently being updated
  int mutable s = 0;

  // s-th element corresponds to the unranked state with bits from the first s
  // sectors having been deposited to it
  std::vector<sv_index_type> mutable partial_unranked;

public:
  explicit multisector_unranking_generator(
      detail::sector_params_vec_t const& params)
    : sector_params(params),
      sector_gens(params.begin(), params.end()),
      total(std::accumulate(sector_gens.begin(),
                            sector_gens.end(),
                            1,
                            [](sv_index_type n, unranking_generator const& g) {
                              return n * g.size();
                            })),
      partial_unranked(sector_gens.size(), 0) {}

  inline void init() const {
    for(auto const& g : sector_gens)
      g.init();
    i = 0;
    s = 0;
  }

  inline sv_index_type next() const {
    if(sector_gens.size() == 0) {
      s = -1;
      ++i;
      return 0;
    }

    while(s >= 0) {
      if(sector_gens[s].done()) {
        sector_gens[s].init();
        s -= 1;
        continue;
      }

      sv_index_type val =
          detail::deposit_bits(sector_gens[s].next(), sector_params[s].mask);

      if(s == 0)
        partial_unranked[s] = val;
      else
        partial_unranked[s] = partial_unranked[s - 1] + val;

      if(s < sector_gens.size() - 1)
        s += 1;
      else
        break;
    }
    ++i;
    return partial_unranked.back();
  }

  inline bool done() const { return i == total; }

  inline sv_index_type size() const { return total; }
};

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

template <typename StateVector,
          bool Ref = true,
          typename RankingAlgorithm = combination_ranking>
struct n_fermion_multisector_view {

  // The underlying state vector
  typename std::conditional<Ref, StateVector&, StateVector>::type state_vector;

  // Parameters of sectors
  detail::sector_params_vec_t sector_params;

  // State ranking algorithm for each sector
  std::vector<RankingAlgorithm> ranking;

  // Strides used to combine sector indices into a multisector index
  std::vector<sv_index_type> sector_strides;

  // Generator of unranked states
  detail::multisector_unranking_generator unranking;

  // The number of bits not belonging to the multisector
  unsigned int M_nonmultisector;

  // Mask used to select non-multisector bits
  sv_index_type nonmultisector_bits_mask = 0;

  using scalar_type = typename element_type<
      typename std::remove_const<StateVector>::type>::type;

  template <typename SV, typename HSType>
  n_fermion_multisector_view(
      SV&& sv,
      HSType const& hs,
      std::vector<sector_descriptor<HSType>> const& sectors)
    : state_vector(std::forward<SV>(sv)),
      sector_params(detail::make_sector_params(hs, sectors)),
      ranking(sector_params.begin(), sector_params.end()),
      sector_strides(init_sector_strides(sector_params)),
      unranking(sector_params) {
    M_nonmultisector =
        make_nonmultisector_mask(hs, sector_params, nonmultisector_bits_mask);
  }

  sv_index_type map_index(sv_index_type index) const {
    sv_index_type index_multisector = 0;
    for(std::size_t s = 0; s < sector_params.size(); ++s) {
      sv_index_type ranked =
          ranking[s](detail::extract_bits(index, sector_params[s].mask));
      index_multisector += ranked * sector_strides[s];
    }
    sv_index_type index_nonmultisector =
        detail::extract_bits(index, nonmultisector_bits_mask);

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
      strides[i] = strides[i + 1] * detail::binomial(sector_params[i + 1].M,
                                                     sector_params[i + 1].N);

    return strides;
  }
};

// Get element type of the StateVector object adapted by a given
// n_fermion_multisector_view object.
template <typename StateVector, bool Ref, typename RankingAlgorithm>
struct element_type<
    n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm>> {
  using type =
      typename n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm>::
          scalar_type;
};

// Get state amplitude of the adapted StateVector object
// at index view.map_index(n)
template <typename StateVector, bool Ref, typename RankingAlgorithm>
inline auto get_element(
    n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm> const& view,
    sv_index_type n) ->
    typename n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm>::
        scalar_type {
  return get_element(view.state_vector, view.map_index(n));
}

// Add a constant to a state amplitude stored in the adapted StateVector object
// at index view.map_index(n)
template <typename StateVector, bool Ref, typename RankingAlgorithm, typename T>
inline void update_add_element(
    n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm>& view,
    sv_index_type n,
    T&& value) {
  update_add_element(view.state_vector,
                     view.map_index(n),
                     std::forward<T>(value));
}

// update_add_element() is not defined for constant views
template <typename StateVector, bool Ref, typename RankingAlgorithm, typename T>
inline void update_add_element(
    n_fermion_multisector_view<StateVector const, Ref, RankingAlgorithm>&,
    sv_index_type,
    T&&) { // NOLINT(cppcoreguidelines-missing-std-forward)
  static_assert(false,
                "update_add_element() is not supported for constant views");
}

// zeros_like() is not defined for views
template <typename StateVector, bool Ref, typename RankingAlgorithm>
inline StateVector zeros_like(
    n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm> const&) {
  static_assert(false, "zeros_like() is not supported for views");
}

// Set all amplitudes stored in the adapted StateVector object to zero
template <typename StateVector, bool Ref, typename RankingAlgorithm>
inline void set_zeros(
    n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm>& view) {
  set_zeros(view.state_vector);
}

// set_zeros() is not defined for constant views
template <typename StateVector, bool Ref, typename RankingAlgorithm>
inline void set_zeros(
    n_fermion_multisector_view<StateVector const, Ref, RankingAlgorithm>&) {
  static_assert(false, "set_zeros() is not supported for constant views");
}

// Apply functor `f` to all index/non-zero amplitude pairs
// in the adapted StateVector object
template <typename StateVector,
          bool Ref,
          typename RankingAlgorithm,
          typename Functor>
inline void foreach(
    n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm> const& view,
    Functor&& f) {
  if(view.sector_params.size() == 0 && view.M_nonmultisector == 0) return;

  auto dim_nonmultisector = detail::pow2(view.M_nonmultisector);
  sv_index_type multisector_index = 0;

  view.unranking.init();
  while(!view.unranking.done()) {
    sv_index_type unranked = view.unranking.next();
    for(sv_index_type index_nms = 0; index_nms < dim_nonmultisector;
        ++index_nms) {

      // Emulate decltype(auto)
      decltype(get_element(view.state_vector, multisector_index)) a =
          get_element(view.state_vector, multisector_index);

      ++multisector_index;

      using T =
          typename n_fermion_multisector_view<StateVector,
                                              Ref,
                                              RankingAlgorithm>::scalar_type;
      if(scalar_traits<T>::is_zero(a))
        continue;
      else {
        sv_index_type index =
            unranked +
            detail::deposit_bits(index_nms, view.nonmultisector_bits_mask);
        std::forward<Functor>(f)(index, a);
      }
    }
  }
}

// Make a list of basis state indices spanning the N-fermion multisector
// The order of indices is consistent with that used by
// n_fermion_multisector_view
template <typename HSType, typename RankingAlgorithm = combination_ranking>
inline std::vector<sv_index_type> n_fermion_multisector_basis_states(
    HSType const& hs,
    std::vector<sector_descriptor<HSType>> const& sectors) {

  auto sector_params = detail::make_sector_params(hs, sectors);
  if(hs.total_n_bits() == 0 && sector_params.empty()) return {};

  detail::multisector_unranking_generator unranking(sector_params);

  sv_index_type nonmultisector_bits_mask = 0;
  unsigned int M_nonmultisector =
      make_nonmultisector_mask(hs, sector_params, nonmultisector_bits_mask);

  sv_index_type dim_nonmultisector = detail::pow2(M_nonmultisector);

  std::vector<sv_index_type> basis_states;
  sv_index_type dim_multisector_fermion = unranking.size();
  basis_states.reserve(dim_multisector_fermion * dim_nonmultisector);

  unranking.init();
  while(!unranking.done()) {
    sv_index_type unranked = unranking.next();
    for(sv_index_type index_nms = 0; index_nms < dim_nonmultisector;
        ++index_nms) {

      sv_index_type index =
          unranked + detail::deposit_bits(index_nms, nonmultisector_bits_mask);
      basis_states.push_back(index);
    }
  }

  return basis_states;
}

template <typename StateVector, typename RankingAlgorithm>
using make_nfms_view_ret_t =
    n_fermion_multisector_view<remove_cvref_t<StateVector>,
                               std::is_lvalue_reference<StateVector>::value,
                               RankingAlgorithm>;

// Make a non-constant N-fermion sector view
template <typename StateVector,
          typename HSType,
          typename RankingAlgorithm = combination_ranking>
auto make_nfms_view(StateVector&& sv,
                    HSType const& hs,
                    std::vector<sector_descriptor<HSType>> const& sectors)
    -> make_nfms_view_ret_t<StateVector, RankingAlgorithm> {
  return make_nfms_view_ret_t<StateVector, RankingAlgorithm>(
      std::forward<StateVector>(sv),
      hs,
      sectors);
}

template <typename StateVector, typename RankingAlgorithm>
using make_const_nfms_view_ret_t =
    n_fermion_multisector_view<remove_cvref_t<StateVector> const,
                               std::is_lvalue_reference<StateVector>::value,
                               RankingAlgorithm>;

// Make a constant N-fermion sector view
template <typename StateVector,
          typename HSType,
          typename RankingAlgorithm = combination_ranking>
auto make_const_nfms_view(StateVector&& sv,
                          HSType const& hs,
                          std::vector<sector_descriptor<HSType>> const& sectors)
    -> make_const_nfms_view_ret_t<StateVector, RankingAlgorithm> {
  return make_const_nfms_view_ret_t<StateVector, RankingAlgorithm>(
      std::forward<StateVector>(sv),
      hs,
      sectors);
}

} // namespace libcommute

#endif
