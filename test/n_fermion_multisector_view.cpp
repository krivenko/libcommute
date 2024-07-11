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

#include <catch.hpp>

#include <libcommute/expression/factories.hpp>
#include <libcommute/loperator/elementary_space_boson.hpp>
#include <libcommute/loperator/elementary_space_fermion.hpp>
#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/loperator.hpp>
#include <libcommute/loperator/n_fermion_sector_view.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

using namespace libcommute;

using hs_type = hilbert_space<int>;
using sd_type = sector_descriptor<hs_type>;
using state_vector = std::vector<double>;
using view_type = n_fermion_multisector_view<state_vector>;

// Compute the number of set bits in the given bit range of 'n'
unsigned int popcount(sv_index_type n,
                      std::pair<unsigned int, unsigned int> const& bit_range) {
  unsigned int count = 0;
  for(unsigned int b = 0; b <= bit_range.second; ++b) {
    if(b >= bit_range.first) count += n & sv_index_type(1);
    n >>= 1;
  }
  return count;
}

// Check that values returned by map_index() form a continuous sequence
// [0; expected_multisector_size)
template <typename StateSelector>
void check_map_index(view_type const& view,
                     hs_type const& hs,
                     StateSelector&& selector,
                     sv_index_type expected_multisector_size) {
  std::vector<sv_index_type> mapped_indices;
  for(sv_index_type index = 0; index < hs.dim(); ++index) {
    if(std::forward<StateSelector>(selector)(index)) {
      mapped_indices.push_back(view.map_index(index));
    }
  }
  std::sort(mapped_indices.begin(), mapped_indices.end());
  std::vector<sv_index_type> mapped_indices_ref(expected_multisector_size);
  std::iota(mapped_indices_ref.begin(), mapped_indices_ref.end(), 0);
  CHECK(mapped_indices == mapped_indices_ref);
}

// Build a reference list of basis states for
// n_fermion_multisector_basis_states() by repeatedly calling map_index()
template <typename StateSelector>
std::vector<sv_index_type>
build_basis_states_ref(hs_type const& hs,
                       std::vector<sd_type> const& sectors,
                       StateSelector&& selector) {
  state_vector st;
  auto view = view_type(st, hs, sectors);
  std::vector<sv_index_type> basis_states(
      n_fermion_multisector_size(hs, sectors));
  for(sv_index_type index = 0; index < hs.dim(); ++index) {
    if(std::forward<StateSelector>(selector)(index)) {
      basis_states[view.map_index(index)] = index;
    }
  }
  return basis_states;
}

TEST_CASE("Implementation details", "[detail]") {
  using namespace detail;

  SECTION("deposit_bits") {
    CHECK(deposit_bits(0x3F, 0x07020408) == 0x07020408);
    CHECK(deposit_bits(0x38, 0x07020408) == 0x07000000);
    CHECK(deposit_bits(0x0F, 0x07020408) == 0x01020408);
  }

  SECTION("extract_bits") {
    CHECK(extract_bits(0x00FF00FF, 0x07020408) == 5 /* 0b101 */);
    CHECK(extract_bits(0xFF00FF00, 0x07020408) == 58 /* 0b111010 */);
    CHECK(extract_bits(0xFF00000000000000, 0x0100000000000000) == 1);
  }

  SECTION("multisector_unranking_generator") {
    using namespace static_indices;

    using gen_t = multisector_unranking_generator;
    using params_t = detail::n_fermion_sector_params_t;

    auto check_output = [](gen_t const& g,
                           std::vector<sv_index_type> const& unranked_ref) {
      for(int i = 0; i < 2; ++i) {
        std::vector<sv_index_type> unranked;
        g.init();
        while(!g.done()) {
          unranked.emplace_back(g.next());
        }
        CHECK(unranked == unranked_ref);
      }
    };

    hs_type hs;

    // Empty Hilbert space
    check_output(gen_t({}), {0});

    // Empty sectors
    check_output(gen_t({params_t(hs, sd_type{{}, 0})}), {0});
    check_output(
        gen_t({params_t(hs, sd_type{{}, 0}), params_t(hs, sd_type{{}, 0})}),
        {0});

    // Two sectors
    for(unsigned int i = 0; i < 12; ++i)
      hs.add(make_space_fermion(int(i)));

    auto sda = sd_type{{{1}, {2}, {6}, {7}}, 2};
    auto sdb = sd_type{{{3}, {4}, {9}}, 2};

    // Sector A (N = 2, count_occupied == true)
    // 7 6 2 1
    //
    // 0 0 1 1
    // 0 1 0 1
    // 0 1 1 0
    // 1 0 0 1
    // 1 0 1 0
    // 1 1 0 0
    std::vector<sv_index_type> ref_a = {(1 << 1) + (1 << 2),
                                        (1 << 1) + (1 << 6),
                                        (1 << 2) + (1 << 6),
                                        (1 << 1) + (1 << 7),
                                        (1 << 2) + (1 << 7),
                                        (1 << 6) + (1 << 7)};
    // Sector B (N = 2, count_occupied == false)
    // 9 4 3
    //
    // 1 1 0
    // 1 0 1
    // 0 1 1
    std::vector<sv_index_type> ref_b = {(1 << 4) + (1 << 9),
                                        (1 << 3) + (1 << 9),
                                        (1 << 3) + (1 << 4)};
    std::vector<sv_index_type> ref;
    for(auto i : ref_a) {
      for(auto j : ref_b) {
        ref.emplace_back(i + j); // cppcheck-suppress useStlAlgorithm
      }
    }
    check_output(gen_t({params_t(hs, sda), params_t(hs, sdb)}), ref);
  }
}

TEST_CASE("View of a state vector projected on a direct product of "
          "multiple N-fermion sectors",
          "[n_fermion_multisector_view]") {

  using namespace static_indices;

  using detail::pow2;
  using detail::binomial;

  unsigned int const M_total = 11;

  unsigned int const N5_max = 1;
  unsigned int const Na_max = 4;
  unsigned int const Nb_max = 4;

  auto sde = [](unsigned int N) { return sd_type{{}, N}; };
  auto sd0 = [](unsigned int N) { return sd_type{{0}, N}; };
  auto sd5 = [](unsigned int N) { return sd_type{{5}, N}; };
  auto sda = [](unsigned int N) { return sd_type{{{1}, {2}, {6}, {7}}, N}; };
  auto sdb = [](unsigned int N) { return sd_type{{{3}, {4}, {8}, {9}}, N}; };

  auto sd5_index_selector = [](unsigned int N) {
    return [N](sv_index_type index) { return popcount(index, {5, 5}) == N; };
  };

  auto sdab_index_selector = [](unsigned int N1, unsigned int N2) {
    return [N1, N2](sv_index_type index) {
      return (popcount(index, {1, 2}) + popcount(index, {6, 7}) == N1) &&
             (popcount(index, {3, 4}) + popcount(index, {8, 9}) == N2);
    };
  };

  auto sda5b_index_selector =
      [](unsigned int N1, unsigned int N2, unsigned int N3) {
        return [N1, N2, N3](sv_index_type index) {
          return (popcount(index, {1, 2}) + popcount(index, {6, 7}) == N1) &&
                 (popcount(index, {5, 5}) == N2) &&
                 (popcount(index, {3, 4}) + popcount(index, {8, 9}) == N3);
        };
      };

  std::vector<sv_index_type> ab_size_ref = {1, 4, 6, 4, 1};

  SECTION("n_fermion_multisector_size") {
    hs_type hs;

    // Empty Hilbert space
    CHECK(n_fermion_multisector_size(hs, {}) == 0);
    CHECK_THROWS_AS(n_fermion_multisector_size(hs, {sda(1)}),
                    std::runtime_error);

    // Empty sectors
    CHECK(n_fermion_multisector_size(hs, {sde(0)}) == 0);
    CHECK(n_fermion_multisector_size(hs, {sde(1)}) == 0);
    CHECK(n_fermion_multisector_size(hs, {sde(0), sde(0)}) == 0);
    CHECK(n_fermion_multisector_size(hs, {sde(1), sde(1)}) == 0);

    for(unsigned int i = 0; i < M_total; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      // Empty sectors
      CHECK(n_fermion_multisector_size(hs, {sde(0)}) == 2048);
      CHECK(n_fermion_multisector_size(hs, {sde(1)}) == 0);
      CHECK(n_fermion_multisector_size(hs, {sde(0), sde(0)}) == 2048);
      CHECK(n_fermion_multisector_size(hs, {sde(1), sde(1)}) == 0);

      CHECK_THROWS_AS(
          n_fermion_multisector_size(hs, {sda(2), sd_type{{{6}}, 1}}),
          std::runtime_error);

      CHECK(n_fermion_multisector_size(hs, {}) == 2048); // No sectors

      // One sector
      for(unsigned int N = 0; N <= Na_max; ++N)
        CHECK(n_fermion_multisector_size(hs, {sda(N)}) == ab_size_ref[N] * 128);
      for(unsigned int N = 0; N <= N5_max; ++N)
        CHECK(n_fermion_multisector_size(hs, {sd5(N)}) == 1024);
      for(unsigned int N = 0; N <= Nb_max; ++N)
        CHECK(n_fermion_multisector_size(hs, {sdb(N)}) == ab_size_ref[N] * 128);

      // Two sectors
      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= Nb_max; ++N2) {
          CHECK(n_fermion_multisector_size(hs, {sda(N1), sdb(N2)}) ==
                ab_size_ref[N1] * ab_size_ref[N2] * 8);
        }
      }

      // Three sectors
      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= N5_max; ++N2) {
          for(unsigned int N3 = 0; N3 <= Nb_max; ++N3) {
            CHECK(n_fermion_multisector_size(hs, {sda(N1), sd5(N2), sdb(N3)}) ==
                  ab_size_ref[N1] * ab_size_ref[N3] * 4);
          }
        }
      }

      // Mixture with the empty sector
      for(unsigned int N = 0; N <= Na_max; ++N) {
        CHECK(n_fermion_multisector_size(hs, {sda(N), sde(0)}) ==
              ab_size_ref[N] * 128);
      }
      for(unsigned int N = 0; N <= N5_max; ++N)
        CHECK(n_fermion_multisector_size(hs, {sd5(N), sde(0)}) == 1024);
      for(unsigned int N = 0; N <= Nb_max; ++N) {
        CHECK(n_fermion_multisector_size(hs, {sdb(N), sde(0)}) ==
              ab_size_ref[N] * 128);
      }
    }

    hs.add(make_space_boson(2, 0));

    SECTION("Fermions and bosons") {
      CHECK(n_fermion_multisector_size(hs, {}) == 8192); // No sectors

      // One sector
      for(unsigned int N = 0; N <= Na_max; ++N)
        CHECK(n_fermion_multisector_size(hs, {sda(N)}) == ab_size_ref[N] * 512);
      for(unsigned int N = 0; N <= N5_max; ++N)
        CHECK(n_fermion_multisector_size(hs, {sd5(N)}) == 4096);
      for(unsigned int N = 0; N <= Nb_max; ++N)
        CHECK(n_fermion_multisector_size(hs, {sdb(N)}) == ab_size_ref[N] * 512);

      // Two sectors
      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= Nb_max; ++N2) {
          CHECK(n_fermion_multisector_size(hs, {sda(N1), sdb(N2)}) ==
                ab_size_ref[N1] * ab_size_ref[N2] * 32);
        }
      }

      // Three sectors
      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= N5_max; ++N2) {
          for(unsigned int N3 = 0; N3 <= Nb_max; ++N3) {
            CHECK(n_fermion_multisector_size(hs, {sda(N1), sd5(N2), sdb(N3)}) ==
                  ab_size_ref[N1] * ab_size_ref[N3] * 16);
          }
        }
      }

      // Mixture with the empty sector
      for(unsigned int N = 0; N <= Na_max; ++N) {
        CHECK(n_fermion_multisector_size(hs, {sda(N), sde(0)}) ==
              ab_size_ref[N] * 512);
      }
      for(unsigned int N = 0; N <= N5_max; ++N)
        CHECK(n_fermion_multisector_size(hs, {sd5(N), sde(0)}) == 4096);
      for(unsigned int N = 0; N <= Nb_max; ++N) {
        CHECK(n_fermion_multisector_size(hs, {sdb(N), sde(0)}) ==
              ab_size_ref[N] * 512);
      }
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));
      CHECK(n_fermion_multisector_size(hs_b, {}) == 32);
      CHECK_THROWS_AS(n_fermion_multisector_size(hs_b, {sda(1)}),
                      std::runtime_error);
    }
  }

  SECTION("Construction") {
    hs_type hs;
    state_vector st{};

    auto check_view = [](view_type const& view,
                         std::vector<unsigned int> const& M,
                         std::vector<sv_index_type> const& masks,
                         std::vector<unsigned int> const& N,
                         std::vector<sv_index_type> const& sector_strides,
                         unsigned int M_nonmultisector,
                         sv_index_type nonmultisector_bits_mask) {
      REQUIRE(N.size() == M.size());

      for(std::size_t i = 0; i < M.size(); ++i) {
        CHECK(view.sector_params[i].M == M[i]);
        CHECK(view.sector_params[i].mask == masks[i]);
        CHECK(view.sector_params[i].N == N[i]);
      }

      CHECK(view.sector_strides == sector_strides);
      CHECK(view.M_nonmultisector == M_nonmultisector);
      CHECK(view.nonmultisector_bits_mask == nonmultisector_bits_mask);
    };

    check_view(view_type(st, hs, {}), {}, {}, {}, {}, 0, 0x0);
    check_view(view_type(st, hs, {sde(0)}), {}, {}, {}, {}, 0, 0x0);
    check_view(view_type(st, hs, {sde(0), sde(0)}), {}, {}, {}, {}, 0, 0x0);

    CHECK_THROWS_AS(view_type(st, hs, {sda(0)}), std::runtime_error);
    CHECK_THROWS_AS(view_type(st, hs, {sda(1)}), std::runtime_error);

    hs.add(make_space_fermion(0));

    SECTION("One fermion") {
      check_view(view_type(st, hs, {}), {}, {}, {}, {}, 1, 0x1);
      check_view(view_type(st, hs, {sde(0)}), {}, {}, {}, {}, 1, 0x1);
      check_view(view_type(st, hs, {sde(0), sde(0)}), {}, {}, {}, {}, 1, 0x1);

      check_view(view_type(st, hs, {sd0(0)}), {1}, {0x1}, {0}, {1}, 0, 0x0);
      check_view(view_type(st, hs, {sd0(1)}), {1}, {0x1}, {1}, {1}, 0, 0x0);

      check_view(view_type(st, hs, {sd0(0), sde(0)}),
                 {1},
                 {0x1},
                 {0},
                 {1},
                 0,
                 0x0);
      check_view(view_type(st, hs, {sd0(1), sde(0)}),
                 {1},
                 {0x1},
                 {1},
                 {1},
                 0,
                 0x0);

      check_view(view_type(st, hs, {sde(0), sd0(0)}),
                 {1},
                 {0x1},
                 {0},
                 {1},
                 0,
                 0x0);
      check_view(view_type(st, hs, {sde(0), sd0(1)}),
                 {1},
                 {0x1},
                 {1},
                 {1},
                 0,
                 0x0);

      CHECK_THROWS_AS(view_type(st, hs, {sd0(2)}), std::runtime_error);

      CHECK_THROWS_AS(view_type(st, hs, {sda(0)}), std::runtime_error);
      CHECK_THROWS_AS(view_type(st, hs, {sda(1)}), std::runtime_error);
    }

    for(unsigned int i = 1; i < M_total; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Multiple fermions") {
      sv_index_type full = (1 << M_total) - 1;

      check_view(view_type(st, hs, {}), {}, {}, {}, {}, M_total, full);
      check_view(view_type(st, hs, {sde(0)}), {}, {}, {}, {}, M_total, full);
      check_view(view_type(st, hs, {sde(0), sde(0)}),
                 {},
                 {},
                 {},
                 {},
                 M_total,
                 full);

      unsigned int M = M_total - 1;

      check_view(view_type(st, hs, {sd0(0)}),
                 {1},
                 {0x1},
                 {0},
                 {1},
                 M,
                 full - 1);
      check_view(view_type(st, hs, {sd0(1)}),
                 {1},
                 {0x1},
                 {1},
                 {1},
                 M,
                 full - 1);

      check_view(view_type(st, hs, {sd0(0), sde(0)}),
                 {1},
                 {0x1},
                 {0},
                 {1},
                 M,
                 full - 1);
      check_view(view_type(st, hs, {sd0(1), sde(0)}),
                 {1},
                 {0x1},
                 {1},
                 {1},
                 M,
                 full - 1);

      check_view(view_type(st, hs, {sde(0), sd0(0)}),
                 {1},
                 {0x1},
                 {0},
                 {1},
                 M,
                 full - 1);
      check_view(view_type(st, hs, {sde(0), sd0(1)}),
                 {1},
                 {0x1},
                 {1},
                 {1},
                 M,
                 full - 1);

      check_view(view_type(st, hs, {sd5(0)}),
                 {1},
                 {0x20},
                 {0},
                 {1},
                 M,
                 full - (1 << 5));
      check_view(view_type(st, hs, {sd5(1)}),
                 {1},
                 {0x20},
                 {1},
                 {1},
                 M,
                 full - (1 << 5));

      check_view(view_type(st, hs, {sd5(0), sde(0)}),
                 {1},
                 {0x20},
                 {0},
                 {1},
                 M,
                 full - (1 << 5));
      check_view(view_type(st, hs, {sd5(1), sde(0)}),
                 {1},
                 {0x20},
                 {1},
                 {1},
                 M,
                 full - (1 << 5));

      check_view(view_type(st, hs, {sde(0), sd5(0)}),
                 {1},
                 {0x20},
                 {0},
                 {1},
                 M,
                 full - (1 << 5));
      check_view(view_type(st, hs, {sde(0), sd5(1)}),
                 {1},
                 {0x20},
                 {1},
                 {1},
                 M,
                 full - (1 << 5));

      M = M_total - 2 * 4;

      for(unsigned int N1 = 0; N1 < Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 < Nb_max; ++N2) {
          std::vector<sv_index_type> strides_r = {binomial(4, N2), 1};

          check_view(view_type(st, hs, {sda(N1), sdb(N2)}),
                     {4, 4},
                     {0xc6, 0x318},
                     {N1, N2},
                     strides_r,
                     M,
                     0x421);
          check_view(view_type(st, hs, {sde(0), sda(N1), sdb(N2)}),
                     {4, 4},
                     {0xc6, 0x318},
                     {N1, N2},
                     strides_r,
                     M,
                     0x421);
          check_view(view_type(st, hs, {sda(N1), sde(0), sdb(N2)}),
                     {4, 4},
                     {0xc6, 0x318},
                     {N1, N2},
                     strides_r,
                     M,
                     0x421);
          check_view(view_type(st, hs, {sda(N1), sdb(N2), sde(0)}),
                     {4, 4},
                     {0xc6, 0x318},
                     {N1, N2},
                     strides_r,
                     M,
                     0x421);
          check_view(view_type(st, hs, {sdb(N1), sda(N2)}),
                     {4, 4},
                     {0x318, 0xc6},
                     {N1, N2},
                     strides_r,
                     M,
                     0x421);
          check_view(view_type(st, hs, {sde(0), sdb(N1), sda(N2)}),
                     {4, 4},
                     {0x318, 0xc6},
                     {N1, N2},
                     strides_r,
                     M,
                     0x421);
          check_view(view_type(st, hs, {sdb(N1), sde(0), sda(N2)}),
                     {4, 4},
                     {0x318, 0xc6},
                     {N1, N2},
                     strides_r,
                     M,
                     0x421);
          check_view(view_type(st, hs, {sdb(N1), sda(N2), sde(0)}),
                     {4, 4},
                     {0x318, 0xc6},
                     {N1, N2},
                     strides_r,
                     M,
                     0x421);
        }

        CHECK_THROWS_AS(view_type(st, hs, {sda(N1), sda(N1)}),
                        std::runtime_error);
        CHECK_THROWS_AS(view_type(st, hs, {sde(0), sda(N1), sda(N1)}),
                        std::runtime_error);
        CHECK_THROWS_AS(view_type(st, hs, {sda(N1), sde(0), sda(N1)}),
                        std::runtime_error);
        CHECK_THROWS_AS(view_type(st, hs, {sda(N1), sda(N1), sde(0)}),
                        std::runtime_error);
      }

      CHECK_THROWS_AS(view_type(st, hs, {sd0(2)}), std::runtime_error);
    }

    hs.add(make_space_boson(1, 5));
    hs.add(make_space_boson(2, 6));

    SECTION("Fermions and bosons") {
      unsigned int M = M_total + 3;

      sv_index_type full = (1 << M) - 1;

      check_view(view_type(st, hs, {}), {}, {}, {}, {}, M, full);
      check_view(view_type(st, hs, {sde(0)}), {}, {}, {}, {}, M, full);
      check_view(view_type(st, hs, {sde(0), sde(0)}), {}, {}, {}, {}, M, full);

      check_view(view_type(st, hs, {sd0(0)}),
                 {1},
                 {0x1},
                 {0},
                 {1},
                 M - 1,
                 full - 1);
      check_view(view_type(st, hs, {sd0(1)}),
                 {1},
                 {0x1},
                 {1},
                 {1},
                 M - 1,
                 full - 1);

      check_view(view_type(st, hs, {sd0(0), sde(0)}),
                 {1},
                 {0x1},
                 {0},
                 {1},
                 M - 1,
                 full - 1);
      check_view(view_type(st, hs, {sd0(1), sde(0)}),
                 {1},
                 {0x1},
                 {1},
                 {1},
                 M - 1,
                 full - 1);

      check_view(view_type(st, hs, {sde(0), sd0(0)}),
                 {1},
                 {0x1},
                 {0},
                 {1},
                 M - 1,
                 full - 1);
      check_view(view_type(st, hs, {sde(0), sd0(1)}),
                 {1},
                 {0x1},
                 {1},
                 {1},
                 M - 1,
                 full - 1);

      check_view(view_type(st, hs, {sd5(0)}),
                 {1},
                 {0x20},
                 {0},
                 {1},
                 M - 1,
                 full - (1 << 5));
      check_view(view_type(st, hs, {sd5(1)}),
                 {1},
                 {0x20},
                 {1},
                 {1},
                 M - 1,
                 full - (1 << 5));

      check_view(view_type(st, hs, {sd5(0), sde(0)}),
                 {1},
                 {0x20},
                 {0},
                 {1},
                 M - 1,
                 full - (1 << 5));
      check_view(view_type(st, hs, {sd5(1), sde(0)}),
                 {1},
                 {0x20},
                 {1},
                 {1},
                 M - 1,
                 full - (1 << 5));

      check_view(view_type(st, hs, {sde(0), sd5(0)}),
                 {1},
                 {0x20},
                 {0},
                 {1},
                 M - 1,
                 full - (1 << 5));
      check_view(view_type(st, hs, {sde(0), sd5(1)}),
                 {1},
                 {0x20},
                 {1},
                 {1},
                 M - 1,
                 full - (1 << 5));

      M = M_total + 3 - 2 * 4;

      for(unsigned int N1 = 0; N1 < Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 < Nb_max; ++N2) {
          std::vector<sv_index_type> strides_r = {binomial(4, N2), 1};

          check_view(view_type(st, hs, {sda(N1), sdb(N2)}),
                     {4, 4},
                     {0xc6, 0x318},
                     {N1, N2},
                     strides_r,
                     M,
                     0x3c21);
          check_view(view_type(st, hs, {sde(0), sda(N1), sdb(N2)}),
                     {4, 4},
                     {0xc6, 0x318},
                     {N1, N2},
                     strides_r,
                     M,
                     0x3c21);
          check_view(view_type(st, hs, {sda(N1), sde(0), sdb(N2)}),
                     {4, 4},
                     {0xc6, 0x318},
                     {N1, N2},
                     strides_r,
                     M,
                     0x3c21);
          check_view(view_type(st, hs, {sda(N1), sdb(N2), sde(0)}),
                     {4, 4},
                     {0xc6, 0x318},
                     {N1, N2},
                     strides_r,
                     M,
                     0x3c21);

          check_view(view_type(st, hs, {sdb(N1), sda(N2)}),
                     {4, 4},
                     {0x318, 0xc6},
                     {N1, N2},
                     strides_r,
                     M,
                     0x3c21);
          check_view(view_type(st, hs, {sde(0), sdb(N1), sda(N2)}),
                     {4, 4},
                     {0x318, 0xc6},
                     {N1, N2},
                     strides_r,
                     M,
                     0x3c21);
          check_view(view_type(st, hs, {sdb(N1), sde(0), sda(N2)}),
                     {4, 4},
                     {0x318, 0xc6},
                     {N1, N2},
                     strides_r,
                     M,
                     0x3c21);
          check_view(view_type(st, hs, {sdb(N1), sda(N2), sde(0)}),
                     {4, 4},
                     {0x318, 0xc6},
                     {N1, N2},
                     strides_r,
                     M,
                     0x3c21);
        }

        CHECK_THROWS_AS(view_type(st, hs, {sda(N1), sda(N1)}),
                        std::runtime_error);
        CHECK_THROWS_AS(view_type(st, hs, {sde(0), sda(N1), sda(N1)}),
                        std::runtime_error);
        CHECK_THROWS_AS(view_type(st, hs, {sda(N1), sde(0), sda(N1)}),
                        std::runtime_error);
        CHECK_THROWS_AS(view_type(st, hs, {sda(N1), sda(N1), sde(0)}),
                        std::runtime_error);
      }

      CHECK_THROWS_AS(view_type(st, hs, {sd0(2)}), std::runtime_error);
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      check_view(view_type(st, hs_b, {}), {}, {}, {}, {}, 5, 0x1f);
      check_view(view_type(st, hs_b, {sde(0)}), {}, {}, {}, {}, 5, 0x1f);
      check_view(view_type(st, hs_b, {sde(0), sde(0)}),
                 {},
                 {},
                 {},
                 {},
                 5,
                 0x1f);

      CHECK_THROWS_AS(view_type(st, hs_b, {sd0(0)}), std::runtime_error);
      CHECK_THROWS_AS(view_type(st, hs_b, {sd0(1)}), std::runtime_error);
    }
  }

  SECTION("get_element() and update_add_element()") {
    hs_type hs{make_space_boson(2, 0)};
    for(unsigned int i = 0; i < M_total; ++i)
      hs.add(make_space_fermion(int(i)));

    unsigned int const N1 = 2;
    unsigned int const N2 = 3;

    state_vector st(n_fermion_multisector_size(hs, {sda(N1), sdb(N2)}));
    std::iota(st.begin(), st.end(), 0);

    // 7029 == 11 0 11 01 1 10 10 1
    // Sector A: 0110 (sector state 2)
    // Sector B: 1110 (sector state 0, count_occupied == false)
    // Non-multisector: 11011 == 27
    sv_index_type const index = 7029;
    sv_index_type const sector_index = 27 + ((2 * 4 + 0) << 5);

    SECTION("const") {
      auto view =
          n_fermion_multisector_view<state_vector const>(st,
                                                         hs,
                                                         {sda(N1), sdb(N2)});

      CHECK(std::is_same<element_type_t<decltype(view)>, double>::value);
      CHECK(get_element(view, index) == sector_index);
    }

    SECTION("non-const") {
      auto view =
          n_fermion_multisector_view<state_vector>(st, hs, {sda(N1), sdb(N2)});

      CHECK(std::is_same<element_type_t<decltype(view)>, double>::value);
      CHECK(get_element(view, index) == sector_index);

      state_vector st_ref(st.size());
      std::iota(st_ref.begin(), st_ref.end(), 0);
      st_ref[sector_index] = sector_index - 3;
      update_add_element(view, index, -3);
      CHECK(st == st_ref);

      set_zeros(view);
      CHECK(st == state_vector(st.size(), 0));
    }
  }

  SECTION("map_index()") {
    state_vector st{};

    hs_type hs;

    for(unsigned int i = 0; i < M_total; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      for(unsigned int N = 0; N <= N5_max; ++N) {
        auto view = view_type(st, hs, {sd5(N)});
        check_map_index(view,
                        hs,
                        sd5_index_selector(N),
                        detail::pow2(M_total - 1));
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= Nb_max; ++N2) {
          auto view = view_type(st, hs, {sda(N1), sdb(N2)});
          check_map_index(view,
                          hs,
                          sdab_index_selector(N1, N2),
                          binomial(4, N1) * binomial(4, N2) *
                              pow2(M_total - 8));
        }
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= N5_max; ++N2) {
          for(unsigned int N3 = 0; N3 <= Nb_max; ++N3) {
            auto view = view_type(st, hs, {sda(N1), sd5(N2), sdb(N3)});
            check_map_index(view,
                            hs,
                            sda5b_index_selector(N1, N2, N3),
                            binomial(4, N1) * binomial(4, N3) *
                                pow2(M_total - 9));
          }
        }
      }
    }

    hs.add(make_space_boson(2, int(M_total)));
    hs.add(make_space_boson(2, int(M_total + 1)));

    SECTION("Fermions and bosons") {
      unsigned int const M = M_total + 4;

      for(unsigned int N = 0; N <= N5_max; ++N) {
        auto view = view_type(st, hs, {sd5(N)});
        check_map_index(view, hs, sd5_index_selector(N), detail::pow2(M - 1));
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= Nb_max; ++N2) {
          auto view = view_type(st, hs, {sda(N1), sdb(N2)});
          check_map_index(view,
                          hs,
                          sdab_index_selector(N1, N2),
                          binomial(4, N1) * binomial(4, N2) * pow2(M - 8));
        }
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= N5_max; ++N2) {
          for(unsigned int N3 = 0; N3 <= Nb_max; ++N3) {
            auto view = view_type(st, hs, {sda(N1), sd5(N2), sdb(N3)});
            check_map_index(view,
                            hs,
                            sda5b_index_selector(N1, N2, N3),
                            binomial(4, N1) * binomial(4, N3) * pow2(M - 9));
          }
        }
      }
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      auto view = view_type(st, hs_b, {sde(0)});
      for(sv_index_type index = 0; index < hs_b.dim(); ++index) {
        CHECK(view.map_index(index) == index);
      }
    }
  }

  SECTION("foreach()") {
    hs_type hs;

    // foreach() implements a mapping of the sector basis states to the full
    // Hilbert space basis states. Here we check that the mapping is the inverse
    // of map_index().
    auto check_foreach = [](view_type const& view, state_vector& st) {
      // Start from 1 so that none of st's elements is vanishing
      std::iota(st.begin(), st.end(), 1);
      sv_index_type count = 0;
      foreach(view, [&](sv_index_type index, double a) {
        CHECK(sv_index_type(a) == view.map_index(index) + 1);
        ++count;
      });
      // Check that all elements of 'st' have been processed
      CHECK(count == st.size());
    };

    SECTION("Empty Hilbert space") {
      state_vector st{};
      foreach(view_type(st, hs, {}),
              [](sv_index_type, double) { CHECK(false); });
      foreach(view_type(st, hs, {sde(0)}),
              [](sv_index_type, double) { CHECK(false); });
      foreach(view_type(st, hs, {sde(0), sde(0)}),
              [](sv_index_type, double) { CHECK(false); });
    }

    for(unsigned int i = 0; i < M_total; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      for(unsigned int N = 0; N <= N5_max; ++N) {
        state_vector st(n_fermion_multisector_size(hs, {sd5(N)}));
        check_foreach(view_type(st, hs, {sd5(N)}), st);
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= Nb_max; ++N2) {
          state_vector st(n_fermion_multisector_size(hs, {sda(N1), sdb(N2)}));
          check_foreach(view_type(st, hs, {sda(N1), sdb(N2)}), st);
        }
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= N5_max; ++N2) {
          for(unsigned int N3 = 0; N3 <= Nb_max; ++N3) {
            state_vector st(
                n_fermion_multisector_size(hs, {sda(N1), sd5(N2), sdb(N3)}));
            check_foreach(view_type(st, hs, {sda(N1), sd5(N2), sdb(N3)}), st);
          }
        }
      }
    }

    hs.add(make_space_boson(2, int(M_total)));
    hs.add(make_space_boson(2, int(M_total + 1)));

    SECTION("Fermions and bosons") {
      for(unsigned int N = 0; N <= N5_max; ++N) {
        state_vector st(n_fermion_multisector_size(hs, {sd5(N)}));
        check_foreach(view_type(st, hs, {sd5(N)}), st);
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= Nb_max; ++N2) {
          state_vector st(n_fermion_multisector_size(hs, {sda(N1), sdb(N2)}));
          check_foreach(view_type(st, hs, {sda(N1), sdb(N2)}), st);
        }
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= N5_max; ++N2) {
          for(unsigned int N3 = 0; N3 <= Nb_max; ++N3) {
            state_vector st(
                n_fermion_multisector_size(hs, {sda(N1), sd5(N2), sdb(N3)}));
            check_foreach(view_type(st, hs, {sda(N1), sd5(N2), sdb(N3)}), st);
          }
        }
      }
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      state_vector st(n_fermion_multisector_size(hs_b, {}));
      check_foreach(view_type(st, hs_b, {}), st);
    }
  }

  SECTION("n_fermion_multisector_basis_states()") {
    state_vector st{};

    hs_type hs;

    std::vector<sv_index_type> no_states;

    CHECK(n_fermion_multisector_basis_states(hs, {}) == no_states);
    CHECK(n_fermion_multisector_basis_states(hs, {sde(0)}) == no_states);
    CHECK_THROWS_AS(n_fermion_multisector_basis_states(hs, {sde(1)}),
                    std::runtime_error);
    CHECK(n_fermion_multisector_basis_states(hs, {sde(0), sde(0)}) ==
          no_states);
    CHECK_THROWS_AS(n_fermion_multisector_basis_states(hs, {sde(0), sde(1)}),
                    std::runtime_error);
    CHECK_THROWS_AS(n_fermion_multisector_basis_states(hs, {sde(1), sde(0)}),
                    std::runtime_error);
    CHECK_THROWS_AS(n_fermion_multisector_basis_states(hs, {sde(1), sde(1)}),
                    std::runtime_error);
    CHECK_THROWS_AS(n_fermion_multisector_basis_states(hs, {sd0(0)}),
                    std::runtime_error);
    CHECK_THROWS_AS(n_fermion_multisector_basis_states(hs, {sd0(1)}),
                    std::runtime_error);

    for(unsigned int i = 0; i < M_total; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      for(unsigned int N = 0; N <= N5_max; ++N) {
        auto ref = build_basis_states_ref(hs, {sd5(N)}, sd5_index_selector(N));
        CHECK(n_fermion_multisector_basis_states(hs, {sd5(N)}) == ref);
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= Nb_max; ++N2) {
          auto ref = build_basis_states_ref(hs,
                                            {sda(N1), sdb(N2)},
                                            sdab_index_selector(N1, N2));
          CHECK(n_fermion_multisector_basis_states(hs, {sda(N1), sdb(N2)}) ==
                ref);
        }
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= N5_max; ++N2) {
          for(unsigned int N3 = 0; N3 <= Nb_max; ++N3) {
            auto ref = build_basis_states_ref(hs,
                                              {sda(N1), sd5(N2), sdb(N3)},
                                              sda5b_index_selector(N1, N2, N3));
            CHECK(n_fermion_multisector_basis_states(
                      hs,
                      {sda(N1), sd5(N2), sdb(N3)}) == ref);
          }
        }
      }
    }

    hs.add(make_space_boson(2, int(M_total)));
    hs.add(make_space_boson(2, int(M_total + 1)));

    SECTION("Fermions and bosons") {
      for(unsigned int N = 0; N <= N5_max; ++N) {
        auto ref = build_basis_states_ref(hs, {sd5(N)}, sd5_index_selector(N));
        CHECK(n_fermion_multisector_basis_states(hs, {sd5(N)}) == ref);
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= Nb_max; ++N2) {
          auto ref = build_basis_states_ref(hs,
                                            {sda(N1), sdb(N2)},
                                            sdab_index_selector(N1, N2));
          CHECK(n_fermion_multisector_basis_states(hs, {sda(N1), sdb(N2)}) ==
                ref);
        }
      }

      for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
        for(unsigned int N2 = 0; N2 <= N5_max; ++N2) {
          for(unsigned int N3 = 0; N3 <= Nb_max; ++N3) {
            auto ref = build_basis_states_ref(hs,
                                              {sda(N1), sd5(N2), sdb(N3)},
                                              sda5b_index_selector(N1, N2, N3));
            CHECK(n_fermion_multisector_basis_states(
                      hs,
                      {sda(N1), sd5(N2), sdb(N3)}) == ref);
          }
        }
      }
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      std::vector<sv_index_type> ref(n_fermion_multisector_size(hs_b, {}));
      std::iota(ref.begin(), ref.end(), 0);

      CHECK(n_fermion_multisector_basis_states(hs_b, {}) == ref);
      CHECK(n_fermion_multisector_basis_states(hs_b, {sde(0)}) == ref);
      CHECK_THROWS_AS(n_fermion_multisector_basis_states(hs_b, {sde(1)}),
                      std::runtime_error);
      CHECK(n_fermion_multisector_basis_states(hs_b, {sde(0), sde(0)}) == ref);
      CHECK_THROWS_AS(
          n_fermion_multisector_basis_states(hs_b, {sde(0), sde(1)}),
          std::runtime_error);
      CHECK_THROWS_AS(
          n_fermion_multisector_basis_states(hs_b, {sde(1), sde(0)}),
          std::runtime_error);
      CHECK_THROWS_AS(
          n_fermion_multisector_basis_states(hs_b, {sde(1), sde(1)}),
          std::runtime_error);
      CHECK_THROWS_AS(n_fermion_multisector_basis_states(hs_b, {sd0(0)}),
                      std::runtime_error);
      CHECK_THROWS_AS(n_fermion_multisector_basis_states(hs_b, {sd0(1)}),
                      std::runtime_error);
    }
  }

  SECTION("make_nfms_view() and make_const_nfms_view()") {
    hs_type hs{make_space_fermion(0), make_space_fermion(1)};

    auto sd = sd_type{{{0}, {1}}, 1};

    state_vector st{};

    state_vector& st_ref = st;
    state_vector const& st_cref = st;
    state_vector&& st_rref1 = state_vector{1, 2, 3, 4, 5, 6};
    state_vector&& st_rref2 = state_vector{1, 2, 3, 4, 5, 6};

    auto view_st = make_nfms_view(st, hs, {sd});
    auto view_st_ref = make_nfms_view(st_ref, hs, {sd});
    auto view_tmp = make_nfms_view(state_vector{1, 2, 3, 4, 5, 6}, hs, {sd});
    auto view_rref = make_nfms_view(std::move(st_rref1), hs, {sd});

    CHECK(std::is_same<decltype(view_st),
                       n_fermion_multisector_view<state_vector, true>>::value);
    CHECK(std::is_same<decltype(view_st_ref),
                       n_fermion_multisector_view<state_vector, true>>::value);
    CHECK(std::is_same<decltype(view_tmp),
                       n_fermion_multisector_view<state_vector, false>>::value);
    CHECK(std::is_same<decltype(view_rref),
                       n_fermion_multisector_view<state_vector, false>>::value);

    auto cview_st = make_const_nfms_view(st, hs, {sd});
    auto cview_st_ref = make_const_nfms_view(st_ref, hs, {sd});
    auto cview_st_cref = make_const_nfms_view(st_cref, hs, {sd});
    auto cview_tmp =
        make_const_nfms_view(state_vector{1, 2, 3, 4, 5, 6}, hs, {sd});
    auto cview_rref = make_const_nfms_view(std::move(st_rref2), hs, {sd});

    CHECK(std::is_same<
          decltype(cview_st),
          n_fermion_multisector_view<state_vector const, true>>::value);
    CHECK(std::is_same<
          decltype(cview_st_ref),
          n_fermion_multisector_view<state_vector const, true>>::value);
    CHECK(std::is_same<
          decltype(cview_st_cref),
          n_fermion_multisector_view<state_vector const, true>>::value);
    CHECK(std::is_same<
          decltype(cview_tmp),
          n_fermion_multisector_view<state_vector const, false>>::value);
    CHECK(std::is_same<
          decltype(cview_rref),
          n_fermion_multisector_view<state_vector const, false>>::value);
  }

  SECTION("loperator") {
    hs_type hs{make_space_boson(2, 0)};

    for(unsigned int i = 0; i < M_total; ++i)
      hs.add(make_space_fermion(int(i)));

    std::vector<unsigned int> sector_a_modes = {1, 2, 6, 7};
    std::vector<unsigned int> sector_b_modes = {3, 4, 8, 9};

    auto Ha = (n(1) + n(2) + n(6) + n(7)) * (a_dag(0) + a(0));
    auto Hb = (n(3) + n(4) + n(8) + n(9)) * (a_dag(0) + a(0));
    auto Hopa = make_loperator(Ha, hs);
    auto Hopb = make_loperator(Hb, hs);

    for(unsigned int N1 = 0; N1 <= Na_max; ++N1) {
      for(unsigned int N2 = 0; N2 <= Nb_max; ++N2) {
        auto sectors = std::vector<sd_type>{sda(N1), sdb(N2)};

        state_vector in(n_fermion_multisector_size(hs, sectors), 0);
        n_fermion_multisector_view<state_vector const> view_in(in, hs, sectors);
        state_vector out(n_fermion_multisector_size(hs, sectors));
        n_fermion_multisector_view<state_vector> view_out(out, hs, sectors);

        // Input:
        // First N1 modes of sector A are occupied
        // First N2 modes of sector B are occupied
        // 1 boson

        sv_index_type index_in_f = 0;
        for(unsigned int n1 = 0; n1 < N1; ++n1)
          index_in_f += detail::pow2(sector_a_modes[n1]);
        for(unsigned int n2 = 0; n2 < N2; ++n2)
          index_in_f += detail::pow2(sector_b_modes[n2]);

        in[view_in.map_index(index_in_f + detail::pow2(M_total))] = 1;

        Hopa(view_in, view_out);
        state_vector outa_ref(n_fermion_multisector_size(hs, sectors), 0);
        // 0 bosons
        outa_ref[view_out.map_index(index_in_f)] = N1;
        // 2 bosons
        outa_ref[view_out.map_index(index_in_f + detail::pow2(M_total + 1))] =
            N1 * std::sqrt(2);
        CHECK_THAT(out, Catch::Matchers::Approx(outa_ref));

        Hopb(view_in, view_out);
        state_vector outb_ref(n_fermion_multisector_size(hs, sectors), 0);
        // 0 bosons
        outb_ref[view_out.map_index(index_in_f)] = N2;
        // 2 bosons
        outb_ref[view_out.map_index(index_in_f + detail::pow2(M_total + 1))] =
            N2 * std::sqrt(2);
        CHECK_THAT(out, Catch::Matchers::Approx(outb_ref));
      }
    }
  }
}
