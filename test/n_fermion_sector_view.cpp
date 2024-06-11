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
#include <type_traits>
#include <utility>
#include <vector>

using namespace libcommute;

// Compute the number of set bits in 'n' among the first 'M' bits
unsigned int popcount(sv_index_type n, unsigned int M) {
  unsigned int count = 0;
  for(; M > 0; --M) {
    count += n & sv_index_type(1);
    n >>= 1;
  }
  return count;
}

TEST_CASE("Implementation details", "[detail]") {
  using namespace static_indices;
  using namespace detail;

  SECTION("2^n") {
    CHECK(pow2(0) == 1);
    CHECK(pow2(1) == 2);
    CHECK(pow2(4) == 16);
  }

  SECTION("Binomial coefficient") {
    CHECK(binomial(0, 0) == 1);
    CHECK(binomial(0, 1) == 0);
    CHECK(binomial(0, 2) == 0);
    CHECK(binomial(1, 0) == 1);
    CHECK(binomial(1, 1) == 1);
    CHECK(binomial(1, 2) == 0);

    CHECK(binomial(10, 4) == 210);
    CHECK(binomial(10, 5) == 252);
    CHECK(binomial(10, 6) == 210);
  }

  SECTION("binomial_sum_t") {

    auto check_b_sum = [](binomial_sum_t const& b_sum,
                          unsigned int M,
                          unsigned int N,
                          std::vector<sv_index_type> const& coeffs,
                          unsigned int coeffs_stride) {
      CHECK(b_sum.coeffs == coeffs);
      CHECK(b_sum.coeffs_stride == coeffs_stride);
      for(unsigned int n = 1; n < N; ++n) {
        for(unsigned int m = 1; m < M; ++m) {
          for(unsigned int lambda = 1; lambda < m; ++lambda) {
            sv_index_type ref = 0;
            for(unsigned int j = 2; j <= lambda; ++j)
              ref += binomial(m + 1 - j, n - 1);
            CHECK(b_sum(n, m, lambda) == ref);
          }
        }
      }
    };

    check_b_sum(binomial_sum_t(1, 0), 1, 0, {}, 2);

    check_b_sum(binomial_sum_t(5, 0), 5, 0, {}, 6);
    check_b_sum(binomial_sum_t(5, 1), 5, 1, {1, 1, 1, 1, 1}, 5);
    check_b_sum(binomial_sum_t(5, 2), 5, 2, {1, 1, 1, 1, 1, 2, 3, 4}, 4);

    check_b_sum(binomial_sum_t(6, 0), 6, 0, {}, 7);
    check_b_sum(binomial_sum_t(6, 1), 6, 1, {1, 1, 1, 1, 1, 1}, 6);
    check_b_sum(binomial_sum_t(6, 2), 6, 2, {1, 1, 1, 1, 1, 1, 2, 3, 4, 5}, 5);
    check_b_sum(binomial_sum_t(6, 3),
                6,
                3,
                {1, 1, 1, 1, 1, 2, 3, 4, 1, 3, 6, 10},
                4);
  }
}

TEST_CASE("Ranking and unranking algorithms", "[ranking_unranking]") {

  using namespace static_indices;

  using hs_type = hilbert_space<int>;

  SECTION("n_fermion_sector_params_t") {

    hs_type hs;

    auto check_params = [](n_fermion_sector_params_t const& params,
                           unsigned int M,
                           unsigned int N) {
      CHECK(params.M == M);
      CHECK(params.N == N);
    };

    check_params(n_fermion_sector_params_t(hs, 0), 0, 0);

    hs.add(make_space_fermion(0));

    SECTION("One fermion") {
      for(unsigned int N = 0; N < 1; ++N)
        check_params(n_fermion_sector_params_t(hs, N), 1, N);
    }

    for(unsigned int i = 1; i < 5; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Multiple fermions") {
      for(unsigned int N = 0; N < 6; ++N)
        check_params(n_fermion_sector_params_t(hs, N), 5, N);

      CHECK_THROWS_AS(n_fermion_sector_params_t(hs, 5 + 1), std::runtime_error);
    }

    hs.add(make_space_boson(2, 5));
    hs.add(make_space_boson(3, 6));

    SECTION("Fermions and bosons") {
      for(unsigned int N = 0; N < 6; ++N)
        check_params(n_fermion_sector_params_t(hs, N), 5, N);

      CHECK_THROWS_AS(n_fermion_sector_params_t(hs, 5 + 1), std::runtime_error);
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      check_params(n_fermion_sector_params_t(hs_b, 0), 0, 0);
      CHECK_THROWS_AS(n_fermion_sector_params_t(hs_b, 1), std::runtime_error);
    }
  }

  SECTION("combination_ranking") {

    hs_type hs;

    CHECK(combination_ranking(n_fermion_sector_params_t(hs, 0))(0x0) == 0);

    hs.add(make_space_fermion(0));

    SECTION("One fermion") {
      combination_ranking rank0(n_fermion_sector_params_t(hs, 0));
      CHECK(rank0(0x0) == 0);
      combination_ranking rank1(n_fermion_sector_params_t(hs, 1));
      CHECK(rank1(0x0) == 0);
    }

    for(unsigned int i = 1; i < 5; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Multiple fermions") {
      sv_index_type full = (1 << 5) - 1;

      combination_ranking rank0(n_fermion_sector_params_t(hs, 0));
      combination_ranking rank5(n_fermion_sector_params_t(hs, 5));
      CHECK(rank0(0x0) == 0);
      CHECK(rank5(full) == 0);

      combination_ranking rank1(n_fermion_sector_params_t(hs, 1));
      combination_ranking rank4(n_fermion_sector_params_t(hs, 4));
      for(unsigned int i = 0; i < 5; ++i) {
        CHECK(rank1(1 << i) == i);
        CHECK(rank4(full ^ (1 << i)) == i);
      }

      combination_ranking rank2(n_fermion_sector_params_t(hs, 2));
      combination_ranking rank3(n_fermion_sector_params_t(hs, 3));
      unsigned int i = 0;
      for(unsigned int i1 = 0; i1 < 5 - 1; ++i1) {
        for(unsigned int i2 = i1 + 1; i2 < 5; ++i2) {
          CHECK(rank2((1 << i1) + (1 << i2)) == i);
          CHECK(rank3(full ^ ((1 << i1) + (1 << i2))) == i);
          ++i;
        }
      }
    }
  }

  SECTION("unranking_generator") {
    auto check_output = [](unranking_generator const& g,
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

    check_output(unranking_generator(n_fermion_sector_params_t(hs, 0)),
                 std::vector<sv_index_type>{0});

    hs.add(make_space_fermion(0));

    SECTION("One fermion") {
      check_output(unranking_generator(n_fermion_sector_params_t(hs, 0)),
                   std::vector<sv_index_type>{0x0});
      check_output(unranking_generator(n_fermion_sector_params_t(hs, 1)),
                   std::vector<sv_index_type>{0x1});
    }

    for(unsigned int i = 1; i < 5; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Multiple fermions") {
      sv_index_type full = (1 << 5) - 1;

      check_output(unranking_generator(n_fermion_sector_params_t(hs, 0)),
                   std::vector<sv_index_type>{0x0});
      check_output(unranking_generator(n_fermion_sector_params_t(hs, 5)),
                   std::vector<sv_index_type>{0x1F});

      std::vector<sv_index_type> ref1 = {0x1, 0x2, 0x4, 0x8, 0x10};
      check_output(unranking_generator(n_fermion_sector_params_t(hs, 1)), ref1);
      std::transform(ref1.begin(),
                     ref1.end(),
                     ref1.begin(),
                     [full](sv_index_type i) { return full - i; });
      check_output(unranking_generator(n_fermion_sector_params_t(hs, 4)), ref1);

      std::vector<sv_index_type> ref2 =
          {0x3, 0x5, 0x9, 0x11, 0x6, 0xA, 0x12, 0xC, 0x14, 0x18};
      check_output(unranking_generator(n_fermion_sector_params_t(hs, 2)), ref2);
      std::transform(ref2.begin(),
                     ref2.end(),
                     ref2.begin(),
                     [full](sv_index_type i) { return full - i; });
      check_output(unranking_generator(n_fermion_sector_params_t(hs, 3)), ref2);
    }
  }
}

TEST_CASE("View of a state vector projected on a single N-fermion sector",
          "[n_fermion_sector_view]") {

  using namespace static_indices;

  using hs_type = hilbert_space<int>;
  using state_vector = std::vector<double>;
  using view_type = n_fermion_sector_view<state_vector>;

  SECTION("n_fermion_sector_size") {
    hs_type hs;

    // Empty Hilbert space
    CHECK(n_fermion_sector_size(hs, 0) == 0);
    CHECK(n_fermion_sector_size(hs, 1) == 0);

    // Purely fermionic Hilbert spaces
    hs.add(make_space_fermion(0));
    CHECK(n_fermion_sector_size(hs, 0) == 1);
    CHECK(n_fermion_sector_size(hs, 1) == 1);
    hs.add(make_space_fermion(1));
    CHECK(n_fermion_sector_size(hs, 0) == 1);
    CHECK(n_fermion_sector_size(hs, 1) == 2);
    CHECK(n_fermion_sector_size(hs, 2) == 1);
    hs.add(make_space_fermion(2));
    CHECK(n_fermion_sector_size(hs, 0) == 1);
    CHECK(n_fermion_sector_size(hs, 1) == 3);
    CHECK(n_fermion_sector_size(hs, 2) == 3);
    CHECK(n_fermion_sector_size(hs, 3) == 1);

    // Fermions and bosons
    hs.add(make_space_boson(4, 3));
    CHECK(n_fermion_sector_size(hs, 0) == 16);
    CHECK(n_fermion_sector_size(hs, 1) == 48);
    CHECK(n_fermion_sector_size(hs, 2) == 48);
    CHECK(n_fermion_sector_size(hs, 3) == 16);

    // Purely bosonic Hilbert space
    hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));
    CHECK(n_fermion_sector_size(hs_b, 0) == 32);
    CHECK(n_fermion_sector_size(hs_b, 1) == 0);
  }

  SECTION("Construction") {
    unsigned int const M = 5;
    hs_type hs;
    state_vector st{};

    auto check_view = [](view_type const& view,
                         unsigned int M,
                         unsigned int N,
                         unsigned int M_nonfermion) {
      CHECK(view.params.M == M);
      CHECK(view.params.N == N);
      CHECK(view.M_nonfermion == M_nonfermion);
    };

    check_view(view_type(st, hs, 0), 0, 0, 0);

    hs.add(make_space_fermion(0));

    SECTION("One fermion") {
      for(unsigned int N = 0; N <= 1; ++N)
        check_view(view_type(st, hs, N), 1, N, 0);
    }

    for(unsigned int i = 1; i < M; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Multiple fermions") {
      for(unsigned int N = 0; N <= M; ++N)
        check_view(view_type(st, hs, N), M, N, 0);

      CHECK_THROWS_AS(n_fermion_sector_view<state_vector>(st, hs, 6),
                      std::runtime_error);
    }

    hs.add(make_space_boson(2, 5));
    hs.add(make_space_boson(3, 6));

    SECTION("Fermions and bosons") {
      unsigned int const M_bosons = hs.algebra_bit_range(boson).second -
                                    hs.algebra_bit_range(boson).first + 1;

      for(unsigned int N = 0; N <= M; ++N)
        check_view(view_type(st, hs, N), M, N, M_bosons);

      CHECK_THROWS_AS(n_fermion_sector_view<state_vector>(st, hs, M + 1),
                      std::runtime_error);
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));
      unsigned int const M_bosons = hs_b.algebra_bit_range(boson).second -
                                    hs_b.algebra_bit_range(boson).first + 1;

      check_view(view_type(st, hs_b, 0), 0, 0, M_bosons);
      CHECK_THROWS_AS(view_type(st, hs_b, 1), std::runtime_error);
    }
  }

  SECTION("get_element() and update_add_element()") {
    unsigned int const M = 5;
    unsigned int const N = 2;

    hs_type hs{make_space_boson(2, 0)};
    for(unsigned int i = 0; i < M; ++i)
      hs.add(make_space_fermion(int(i)));

    state_vector st(n_fermion_sector_size(hs, N));
    std::iota(st.begin(), st.end(), 0);

    // 70 == 10 00110: 2 bosons + fermions in modes 1 and 2 (sector state 4)
    sv_index_type const index = 70;
    sv_index_type const sector_index = 2 + (4 << 2);

    SECTION("const") {
      auto view = n_fermion_sector_view<state_vector const>(st, hs, N);

      CHECK(std::is_same<element_type_t<decltype(view)>, double>::value);
      CHECK(get_element(view, index) == sector_index);
    }

    SECTION("non-const") {
      auto view = n_fermion_sector_view<state_vector>(st, hs, N);

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
    unsigned int const M = 8;
    state_vector st{};

    hs_type hs;

    // Check that values returned by map_index() form a continuous sequence
    // starting at 0
    auto check_map_index =
        [&hs](view_type const& view, unsigned int M, unsigned int N) {
          std::vector<sv_index_type> mapped_indices;
          for(sv_index_type index = 0; index < hs.dim(); ++index) {
            if(popcount(index, M) == N) {
              mapped_indices.push_back(view.map_index(index));
            }
          }
          std::sort(mapped_indices.begin(), mapped_indices.end());
          std::vector<sv_index_type> mapped_indices_ref(
              n_fermion_sector_size(hs, N));
          std::iota(mapped_indices_ref.begin(), mapped_indices_ref.end(), 0);
          CHECK(mapped_indices == mapped_indices_ref);
        };

    for(unsigned int i = 0; i < M; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      for(unsigned int N = 0; N <= M; ++N) {
        check_map_index(view_type(st, hs, N), M, N);
      }
    }

    hs.add(make_space_boson(2, int(M)));
    hs.add(make_space_boson(2, int(M + 1)));

    SECTION("Fermions and bosons") {
      for(unsigned int N = 0; N <= M; ++N) {
        check_map_index(view_type(st, hs, N), M, N);
      }
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      auto view = view_type(st, hs_b, 0);
      for(sv_index_type index = 0; index < hs_b.dim(); ++index) {
        CHECK(view.map_index(index) == index);
      }
    }
  }

  SECTION("foreach()") {
    unsigned int const M = 8;

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
      foreach(view_type(st, hs, 0),
              [](sv_index_type, double) { CHECK(false); });
    }

    for(unsigned int i = 0; i < M; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      for(unsigned int N = 0; N <= M; ++N) {
        state_vector st(n_fermion_sector_size(hs, N));
        check_foreach(view_type(st, hs, N), st);
      }
    }

    hs.add(make_space_boson(2, int(M)));
    hs.add(make_space_boson(2, int(M + 1)));

    SECTION("Fermions and bosons") {
      for(unsigned int N = 0; N <= M; ++N) {
        state_vector st(n_fermion_sector_size(hs, N));
        check_foreach(view_type(st, hs, N), st);
      }
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      state_vector st(n_fermion_sector_size(hs_b, 0));
      check_foreach(view_type(st, hs_b, 0), st);
    }
  }

  SECTION("n_fermion_sector_basis_states()") {
    unsigned int const M = 8;
    state_vector st{};

    hs_type hs;

    CHECK(n_fermion_sector_basis_states(hs, 0) == std::vector<sv_index_type>{});
    CHECK_THROWS_AS(n_fermion_sector_basis_states(hs, 1), std::runtime_error);

    auto build_basis_states_ref = [&st](hs_type const& hs, unsigned int N) {
      auto view = view_type(st, hs, N);
      std::vector<sv_index_type> basis_states(n_fermion_sector_size(hs, N));
      for(sv_index_type index = 0; index < hs.dim(); ++index) {
        if(popcount(index, M) == N) {
          basis_states[view.map_index(index)] = index;
        }
      }
      return basis_states;
    };

    for(unsigned int i = 0; i < M; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      for(unsigned int N = 0; N <= M; ++N) {
        auto ref = build_basis_states_ref(hs, N);
        CHECK(n_fermion_sector_basis_states(hs, N) == ref);
      }
    }

    hs.add(make_space_boson(2, int(M)));
    hs.add(make_space_boson(2, int(M + 1)));

    SECTION("Fermions and bosons") {
      for(unsigned int N = 0; N <= M; ++N) {
        auto ref = build_basis_states_ref(hs, N);
        CHECK(n_fermion_sector_basis_states(hs, N) == ref);
      }
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      CHECK_THROWS_AS(n_fermion_sector_basis_states(hs_b, 1),
                      std::runtime_error);

      std::vector<sv_index_type> ref(n_fermion_sector_size(hs_b, 0));
      std::iota(ref.begin(), ref.end(), 0);

      CHECK(n_fermion_sector_basis_states(hs_b, 0) == ref);
    }
  }

  SECTION("make_nfs_view() and make_const_nfs_view()") {
    hs_type hs{make_space_fermion(0), make_space_fermion(1)};

    state_vector st{};

    state_vector& st_ref = st;
    state_vector const& st_cref = st;
    state_vector&& st_rref1 = state_vector{1, 2, 3, 4, 5, 6};
    state_vector&& st_rref2 = state_vector{1, 2, 3, 4, 5, 6};

    auto view_st = make_nfs_view(st, hs, 1);
    auto view_st_ref = make_nfs_view(st_ref, hs, 1);
    auto view_tmp = make_nfs_view(state_vector{1, 2, 3, 4, 5, 6}, hs, 1);
    auto view_rref = make_nfs_view(std::move(st_rref1), hs, 1);

    CHECK(std::is_same<decltype(view_st),
                       n_fermion_sector_view<state_vector, true>>::value);
    CHECK(std::is_same<decltype(view_st_ref),
                       n_fermion_sector_view<state_vector, true>>::value);
    CHECK(std::is_same<decltype(view_tmp),
                       n_fermion_sector_view<state_vector, false>>::value);
    CHECK(std::is_same<decltype(view_rref),
                       n_fermion_sector_view<state_vector, false>>::value);

    auto cview_st = make_const_nfs_view(st, hs, 1);
    auto cview_st_ref = make_const_nfs_view(st_ref, hs, 1);
    auto cview_st_cref = make_const_nfs_view(st_cref, hs, 1);
    auto cview_tmp = make_const_nfs_view(state_vector{1, 2, 3, 4, 5, 6}, hs, 1);
    auto cview_rref = make_const_nfs_view(std::move(st_rref2), hs, 1);

    CHECK(std::is_same<decltype(cview_st),
                       n_fermion_sector_view<state_vector const, true>>::value);
    CHECK(std::is_same<decltype(cview_st_ref),
                       n_fermion_sector_view<state_vector const, true>>::value);
    CHECK(std::is_same<decltype(cview_st_cref),
                       n_fermion_sector_view<state_vector const, true>>::value);
    CHECK(
        std::is_same<decltype(cview_tmp),
                     n_fermion_sector_view<state_vector const, false>>::value);
    CHECK(
        std::is_same<decltype(cview_rref),
                     n_fermion_sector_view<state_vector const, false>>::value);
  }

  SECTION("loperator") {
    unsigned int const M = 5;

    hs_type hs{make_space_boson(2, 0)};

    for(unsigned int i = 0; i < M; ++i)
      hs.add(make_space_fermion(int(i)));

    auto H = (n(0) + n(1) + n(2) + n(3) + n(4)) * (a_dag(0) + a(0));
    auto Hop = make_loperator(H, hs);

    for(unsigned int N = 0; N <= M; ++N) {
      state_vector in(n_fermion_sector_size(hs, N), 0);
      auto view_in = n_fermion_sector_view<state_vector const>(in, hs, N);
      state_vector out(n_fermion_sector_size(hs, N));
      auto view_out = n_fermion_sector_view<state_vector>(out, hs, N);

      // 1 boson, fermions in the first N modes
      sv_index_type index_in_f = 0;
      for(unsigned int n = 0; n < N; ++n)
        index_in_f += detail::pow2(n);

      in[view_in.map_index(index_in_f + detail::pow2(M))] = 1;

      Hop(view_in, view_out);
      state_vector out_ref(n_fermion_sector_size(hs, N), 0);
      // 0 bosons
      out_ref[view_out.map_index(index_in_f)] = N;
      // 2 bosons
      out_ref[view_out.map_index(index_in_f + detail::pow2(M + 1))] =
          N * std::sqrt(2);
      CHECK_THAT(out, Catch::Matchers::Approx(out_ref));
    }
  }
}
