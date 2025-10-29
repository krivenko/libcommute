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
#ifndef LIBCOMMUTE_LOPERATOR_N_FERMION_SECTOR_VIEW_HPP_
#define LIBCOMMUTE_LOPERATOR_N_FERMION_SECTOR_VIEW_HPP_

#include "../algebra_ids.hpp"
#include "../utility.hpp"

#include "bit_ops.hpp"
#include "elementary_space_fermion.hpp"
#include "hilbert_space.hpp"
#include "state_vector.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
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

// Binomial coefficient C(n, k)
inline sv_index_type binomial(unsigned int n, unsigned int k) {
  if(k > n) return 0;
  if(k > n / 2) k = n - k;
  sv_index_type C = 1;
  for(unsigned int i = 0; i < k; ++i)
    C = (C * (n - i)) / (i + 1);
  return C;
}

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
    : M(init_m(hs)), mask(detail::pow2(M) - 1), N(validated_n(N)) {}

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
    sv_index_type m = 0;
    for(auto const& i : sector.indices) {
      unsigned int b =
          hs.bit_range(elementary_space_fermion<IndexTypes...>(i)).first;
      m += detail::pow2(b);
    }
    return m;
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

  // Number of counted modes (either occupied or unoccupied)
  unsigned int N_counted;

  // Precomputed binomial coefficients
  std::vector<sv_index_type> binomial;

  // XOR-mask to be applied to input states
  sv_index_type mask;

public:
  explicit combination_ranking(detail::n_fermion_sector_params_t const& params)
    : N_counted((params.N <= params.M / 2) ? params.N : params.M - params.N),
      binomial(params.M * N_counted),
      mask((params.N == N_counted) ? sv_index_type(0) :
                                     (detail::pow2(params.M) - 1)) {
    for(unsigned int m = 0; m < params.M; ++m) {
      for(unsigned int n = 1; n <= N_counted; ++n) {
        binomial[m * N_counted + n - 1] = detail::binomial(m, n);
      }
    }
  }

  // Rank a fermionic many-body state
  inline sv_index_type operator()(sv_index_type index) const {
    index = index ^ mask;
    sv_index_type i = 0;
    unsigned int k = 0;
    while(index != 0) {
      unsigned int const c = detail::count_trailing_zeros(index);
      index = index & (index - 1);
      i += binomial[c * N_counted + k];
      ++k;
    }
    return i;
  }
};

// Staggered combination ranking algorithm
template <unsigned int R> class staggered_ranking {

  static_assert((R > 0) && (R < 64), "Chunk size R must belong to [1; 63]");

  // 2^R
  constexpr static sv_index_type rdim = detail::pow2(R);

  // Number of chunks + 1
  unsigned int mpmax;

  // Positions of heads (N' = 0, \alpha = 0) within lookup_table corresponding
  // to each M'
  std::vector<sv_index_type> mp_block_heads;

  // Table of precomputed values of rank(\alpha, M', N')
  std::vector<sv_index_type> lookup_table;

  // Population count for each of the 2^R length-R states
  std::vector<unsigned int> popcount_table;

  // XOR-mask to be applied to input states
  sv_index_type mask;

  // rank(\alpha, M', N')
  inline sv_index_type
  rank(sv_index_type r, unsigned int Mp, unsigned int Np) const {
    sv_index_type i = 0;
    unsigned int k = 1;
    while(r != 0) {
      unsigned int c = detail::count_trailing_zeros(r);
      r = r & (r - 1);
      i += detail::binomial(Mp + c, Np + k);
      ++k;
    }
    return i;
  }

public:
  explicit staggered_ranking(detail::n_fermion_sector_params_t const& params)
    : mpmax(params.M / R),
      mask((params.N <= params.M / 2) ? sv_index_type(0) :
                                        (detail::pow2(params.M) - 1)) {

    // Fill mp_block_heads
    mp_block_heads.reserve(mpmax);
    for(unsigned int mp = 0; mp <= mpmax; ++mp) {
      mp_block_heads.emplace_back((mp + R * (mp * mp - mp) / 2) * rdim);
    }

    // Fill lookup_table
    lookup_table.reserve((mpmax + 1) * (R * mpmax + 2) * rdim / 2);
    for(unsigned int mp = 0; mp <= mpmax; ++mp) {
      unsigned int const Mp = mp * R;
      for(unsigned int Np = 0; Np <= Mp; ++Np) {
        for(sv_index_type r = 0; r < rdim; ++r) {
          lookup_table.emplace_back(rank(r, Mp, Np));
        }
      }
    }

    // Fill popcount_table
    popcount_table.reserve(rdim);
    for(sv_index_type r = 0; r < rdim; ++r) {
      popcount_table.emplace_back(detail::popcount(r));
    }
  }

  // Rank a fermionic many-body state
  inline sv_index_type operator()(sv_index_type index) const {
    index = index ^ mask;
    sv_index_type i = 0;
    unsigned int Np = 0;
    for(unsigned int mp = 0; index != 0; ++mp) {
      sv_index_type const r = index & (rdim - 1);
      i += lookup_table[mp_block_heads[mp] + Np * rdim + r];
      Np += popcount_table[r];
      index >>= R;
    }
    return i;
  }
};

// Trie-based ranking algorithm
template <unsigned int R> class trie_ranking {

  static_assert((R > 0) && (R < 64), "Chunk size R must belong to [1; 63]");

  // Sizes of bit-string chunks (<=R for the last chunk, R for the rest)
  std::vector<unsigned int> chunk_sizes;

  // Masks used to extract bits from each chunk
  std::vector<sv_index_type> pext_masks;

  // Packed linearized representation of the trie
  std::vector<std::int64_t> trie;

  // XOR-mask to be applied to input states
  sv_index_type mask;

  inline static sv_index_type next_state(sv_index_type state) {
    sv_index_type b = state | (state - 1);
    ++b;
    sv_index_type l = state & (~b);
    l >>= detail::count_trailing_zeros(state) + 1;
    return b | l;
  }

  std::tuple<sv_index_type, sv_index_type, sv_index_type>
  make_trie(unsigned int level,
            sv_index_type rank,
            sv_index_type st,
            sv_index_type k,
            unsigned int M_subtree) {
    unsigned int chunk_size = chunk_sizes[level];

    sv_index_type const P = ~(detail::pow2(M_subtree) - 1);
    sv_index_type const p = st & P;

    sv_index_type const C = detail::pow2(chunk_size) - 1;
    sv_index_type const kp = k;
    k += detail::pow2(chunk_size);

    sv_index_type c = C;
    sv_index_type const c0 = (st >> (M_subtree - chunk_size)) & C;

    while((st & P) == p) {
      c = (st >> (M_subtree - chunk_size)) & C;
      if(level == pext_masks.size() - 1) { // Leaf
        trie[kp + c] = -static_cast<std::int64_t>(rank);
        if(st == 0x0) break; // There is exactly one state in an N = 0 sector
        ++rank;
        st = next_state(st);
      } else { // Subtree

        // Packing: Remove unused elements at the front of the child index
        sv_index_type const C_subtree =
            detail::pow2(chunk_sizes[level + 1]) - 1;
        sv_index_type const f =
            (st >> (M_subtree - chunk_size - chunk_sizes[level + 1])) &
            C_subtree;
        if(k >= f) {
          k -= f;
        }
        assert(k > kp);

        trie[kp + c] = static_cast<std::int64_t>(k);
        std::tie(rank, st, k) =
            make_trie(level + 1, rank, st, k, M_subtree - chunk_size);

        // There is exactly one state in an N = 0 sector
        if(st == 0x0) break;
      }
    }

    // Packing: Remove unused elements at the back of the child index
    sv_index_type const l = C - c;
    sv_index_type const table_head = kp + c0;
    sv_index_type const table_size = detail::pow2(chunk_size) - c0 - l;

    k -= l;
    if((level <= pext_masks.size() - 1) && (l > 0)) {
      // Remove the unused elements and move the subtree backwards by l
      for(sv_index_type kpp = table_head + table_size; kpp < k; ++kpp) {
        trie[kpp] = trie[kpp + l];
      }
      // Update the non-leaf entrees of the subtree
      for(sv_index_type kpp = table_head; kpp < k; ++kpp) {
        if(trie[kpp] > 0) {
          trie[kpp] -= l;
        }
      }
    }

    return {rank, st, k};
  }

public:
  explicit trie_ranking(detail::n_fermion_sector_params_t const& params)
    : chunk_sizes(params.M / R, R),
      trie(detail::pow2(params.M + 1), 0),
      mask((params.N <= params.M / 2) ? sv_index_type(0) :
                                        (detail::pow2(params.M) - 1)) {

    if(params.M % R != 0) chunk_sizes.push_back(params.M % R);

    pext_masks.reserve(chunk_sizes.size());
    unsigned int pext_mask_shift = params.M;
    for(unsigned int size : chunk_sizes) {
      sv_index_type pext_mask = detail::pow2(size) - 1;
      pext_mask_shift -= size;
      pext_masks.emplace_back(pext_mask << pext_mask_shift);
    }

    if(params.M > 0) {
      unsigned int const N_counted =
          (params.N <= params.M / 2) ? params.N : (params.M - params.N);
      sv_index_type k_final = 0;
      std::tie(std::ignore, std::ignore, k_final) =
          make_trie(0, 0x0, detail::pow2(N_counted) - 1, 0, params.M);
      trie.resize(k_final);
      std::for_each(trie.begin(), trie.end(), [](std::int64_t& theta) {
        if(theta < 0) theta = -theta;
      });
    }
  }

  inline sv_index_type operator()(sv_index_type a) const {
    a = a ^ mask;
    sv_index_type i = 0;
    for(auto const& pext_mask : pext_masks) {
      sv_index_type r = detail::extract_bits(a, pext_mask);
      i = trie[i + r];
    }
    return i;
  }
};

// Generator of unranked basis states
class unranking_generator {

  // Number of counted modes (either occupied or unoccupied)
  unsigned int N_counted;

  // Total number of generated unranked states
  sv_index_type total;

  // Index of current unranked state
  sv_index_type mutable i = 0;

  // Current unranked state
  sv_index_type mutable unranked;

  // XOR-mask to be applied to generated states
  sv_index_type mask;

public:
  explicit unranking_generator(detail::n_fermion_sector_params_t const& params)
    : N_counted(params.N <= params.M / 2 ? params.N : params.M - params.N),
      total(detail::binomial(params.M, N_counted)),
      unranked(detail::pow2(N_counted) - 1),
      mask((params.N == N_counted) ? sv_index_type(0) :
                                     (detail::pow2(params.M) - 1)) {}

  // Initialize iteration
  inline void init() const {
    unranked = detail::pow2(N_counted) - 1;
    i = 0;
  }

  // Return the next state
  inline sv_index_type next() const {
    if(i == 0) {
      ++i;
      return unranked ^ mask;
    }

    sv_index_type b = unranked | (unranked - 1);
    ++b;
    sv_index_type l = unranked & (~b);
    l >>= detail::count_trailing_zeros(unranked) + 1;
    unranked = b | l;
    ++i;
    return unranked ^ mask;
  }

  // Are there still states to be returned by next()?
  inline bool done() const { return i == total; }

  // Total number of generated unranked states
  inline sv_index_type size() const { return total; }
};

//
// N-fermion sector
//

// Size of the fermionic sector with N particles
template <typename HSType>
inline sv_index_type n_fermion_sector_size(HSType const& hs, unsigned int N) {
  if(hs.is_sparse())
    throw std::runtime_error(
        "n_fermion_sector_size: sparse Hilbert spaces are not supported");

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
      f_mask(detail::pow2(sector_params.M) - 1),
      M_nonfermion(hs.total_n_bits() - sector_params.M) {
    if(hs.is_sparse())
      throw std::runtime_error(
          "n_fermion_sector_view: sparse Hilbert spaces are not supported");
  }

  sv_index_type map_index(sv_index_type index) const {
    sv_index_type const ranked = ranking(index & f_mask);
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
template <typename StateVector,
          bool Ref,
          typename RankingAlgorithm,
          typename T,
          typename =
              typename std::enable_if<!std::is_const<StateVector>::value>::type>
inline void update_add_element(
    n_fermion_sector_view<StateVector, Ref, RankingAlgorithm>& view,
    sv_index_type n,
    T&& value) {
  update_add_element(view.state_vector,
                     view.map_index(n),
                     std::forward<T>(value));
}

// Set all amplitudes stored in the adapted StateVector object to zero
template <typename StateVector,
          bool Ref,
          typename RankingAlgorithm,
          typename =
              typename std::enable_if<!std::is_const<StateVector>::value>::type>
inline void
set_zeros(n_fermion_sector_view<StateVector, Ref, RankingAlgorithm>& view) {
  set_zeros(view.state_vector);
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

  sv_index_type const dim_nonfermion = detail::pow2(view.M_nonfermion);
  sv_index_type sector_index = 0;

  view.unranking.init();
  while(!view.unranking.done()) {
    sv_index_type const unranked = view.unranking.next();
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
  if(hs.is_sparse())
    throw std::runtime_error("n_fermion_sector_basis_states: sparse Hilbert "
                             "spaces are not supported");

  detail::n_fermion_sector_params_t sector_params(hs, N);

  unsigned int const M_nonfermion = hs.total_n_bits() - sector_params.M;
  if(sector_params.M == 0 && M_nonfermion == 0) return {};

  unranking_generator unranking(sector_params);

  sv_index_type const dim_nonfermion = detail::pow2(M_nonfermion);

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

// Make a non-constant N-fermion sector view
template <typename StateVector, typename HSType>
auto make_nfs_view(StateVector&& sv, HSType const& hs, unsigned int N)
    -> n_fermion_sector_view<remove_cvref_t<StateVector>,
                             std::is_lvalue_reference<StateVector>::value> {
  return {std::forward<StateVector>(sv), hs, N};
}

// Make a constant N-fermion sector view
template <typename StateVector, typename HSType>
auto make_const_nfs_view(StateVector&& sv, HSType const& hs, unsigned int N)
    -> n_fermion_sector_view<remove_cvref_t<StateVector> const,
                             std::is_lvalue_reference<StateVector>::value> {
  return {std::forward<StateVector>(sv), hs, N};
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
  mask = detail::pow2(hs.total_n_bits()) - 1;
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

      if(s < static_cast<int>(sector_gens.size()) - 1)
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

  if(hs.is_sparse())
    throw std::runtime_error(
        "n_fermion_multisector_size: sparse Hilbert spaces are not supported");

  detail::validate_sectors(hs, sectors);

  auto const total_n_bits = hs.total_n_bits();
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
    if(hs.is_sparse())
      throw std::runtime_error("n_fermion_multisector_view: sparse Hilbert "
                               "spaces are not supported");
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
template <typename StateVector,
          bool Ref,
          typename RankingAlgorithm,
          typename T,
          typename =
              typename std::enable_if<!std::is_const<StateVector>::value>::type>
inline void update_add_element(
    n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm>& view,
    sv_index_type n,
    T&& value) {
  update_add_element(view.state_vector,
                     view.map_index(n),
                     std::forward<T>(value));
}

// Set all amplitudes stored in the adapted StateVector object to zero
template <typename StateVector,
          bool Ref,
          typename RankingAlgorithm,
          typename =
              typename std::enable_if<!std::is_const<StateVector>::value>::type>
inline void set_zeros(
    n_fermion_multisector_view<StateVector, Ref, RankingAlgorithm>& view) {
  set_zeros(view.state_vector);
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

  sv_index_type const dim_nonmultisector = detail::pow2(view.M_nonmultisector);
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

  if(hs.is_sparse())
    throw std::runtime_error("n_fermion_multisector_basis_states: sparse "
                             "Hilbert spaces are not supported");

  auto sector_params = detail::make_sector_params(hs, sectors);
  if(hs.total_n_bits() == 0 && sector_params.empty()) return {};

  detail::multisector_unranking_generator unranking(sector_params);

  sv_index_type nonmultisector_bits_mask = 0;
  unsigned int const M_nonmultisector =
      make_nonmultisector_mask(hs, sector_params, nonmultisector_bits_mask);

  sv_index_type const dim_nonmultisector = detail::pow2(M_nonmultisector);

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

// Make a non-constant N-fermion sector view
template <typename StateVector, typename HSType>
auto make_nfms_view(StateVector&& sv,
                    HSType const& hs,
                    std::vector<sector_descriptor<HSType>> const& sectors)
    -> n_fermion_multisector_view<
        remove_cvref_t<StateVector>,
        std::is_lvalue_reference<StateVector>::value> {
  return {std::forward<StateVector>(sv), hs, sectors};
}

// Make a constant N-fermion sector view
template <typename StateVector, typename HSType>
auto make_const_nfms_view(StateVector&& sv,
                          HSType const& hs,
                          std::vector<sector_descriptor<HSType>> const& sectors)
    -> n_fermion_multisector_view<
        remove_cvref_t<StateVector> const,
        std::is_lvalue_reference<StateVector>::value> {
  return {std::forward<StateVector>(sv), hs, sectors};
}

} // namespace libcommute

#endif
