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
#include <vector>

using namespace libcommute;

unsigned int count_fermions(sv_index_type n, unsigned int M) {
  unsigned int count = 0;
  while(M > 0) {
    count += n & sv_index_type(1);
    --M;
    n >>= 1;
  }
  return count;
}

TEST_CASE("Implementation details", "[detail]") {
  using namespace static_indices;
  using namespace detail;

  using hs_type = hilbert_space<int>;

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

  SECTION("n_fermion_sector_params_t") {
    hs_type hs;

    CHECK_THROWS_AS(n_fermion_sector_params_t(hs, 0), std::runtime_error);

    auto check_params = [](n_fermion_sector_params_t const& params,
                           int M,
                           bool count_occupied,
                           int N_counted) {
      CHECK(params.M == M);
      CHECK(params.count_occupied == count_occupied);
      CHECK(params.N_counted == N_counted);
    };

    hs.add(make_space_fermion(0));

    SECTION("One fermion") {
      check_params(n_fermion_sector_params_t(hs, 0), 1, true, 0);
      check_params(n_fermion_sector_params_t(hs, 1), 1, false, 0);
    }

    for(int i : {1, 2, 3, 4})
      hs.add(make_space_fermion(i));

    SECTION("Multiple fermions") {
      check_params(n_fermion_sector_params_t(hs, 0), 5, true, 0);
      check_params(n_fermion_sector_params_t(hs, 1), 5, true, 1);
      check_params(n_fermion_sector_params_t(hs, 2), 5, true, 2);
      check_params(n_fermion_sector_params_t(hs, 3), 5, false, 2);
      check_params(n_fermion_sector_params_t(hs, 4), 5, false, 1);
      check_params(n_fermion_sector_params_t(hs, 5), 5, false, 0);

      CHECK_THROWS_AS(n_fermion_sector_params_t(hs, 6), std::runtime_error);
    }

    hs.add(make_space_boson(2, 5));
    hs.add(make_space_boson(3, 6));

    SECTION("Fermions and bosons") {
      check_params(n_fermion_sector_params_t(hs, 0), 5, true, 0);
      check_params(n_fermion_sector_params_t(hs, 1), 5, true, 1);
      check_params(n_fermion_sector_params_t(hs, 2), 5, true, 2);
      check_params(n_fermion_sector_params_t(hs, 3), 5, false, 2);
      check_params(n_fermion_sector_params_t(hs, 4), 5, false, 1);
      check_params(n_fermion_sector_params_t(hs, 5), 5, false, 0);

      CHECK_THROWS_AS(n_fermion_sector_params_t(hs, 6), std::runtime_error);
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      check_params(n_fermion_sector_params_t(hs_b, 0), 0, true, 0);
      CHECK_THROWS_AS(n_fermion_sector_params_t(hs_b, 1), std::runtime_error);
    }
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
    hs_type hs;
    state_vector st{};

    CHECK_THROWS_AS(n_fermion_sector_view<state_vector>(st, hs, 0),
                    std::runtime_error);

    auto check_view = [](view_type const& view,
                         int M,
                         bool count_occupied,
                         int N_counted,
                         int M_nonfermion) {
      CHECK(view.M == M);
      CHECK(view.count_occupied == count_occupied);
      CHECK(view.N_counted == N_counted);
      CHECK(view.M_nonfermion == M_nonfermion);
    };

    hs.add(make_space_fermion(0));

    SECTION("One fermion") {
      check_view(view_type(st, hs, 0), 1, true, 0, 0);
      check_view(view_type(st, hs, 1), 1, false, 0, 0);
    }

    for(int i : {1, 2, 3, 4})
      hs.add(make_space_fermion(i));

    SECTION("Multiple fermions") {
      check_view(view_type(st, hs, 0), 5, true, 0, 0);
      check_view(view_type(st, hs, 1), 5, true, 1, 0);
      check_view(view_type(st, hs, 2), 5, true, 2, 0);
      check_view(view_type(st, hs, 3), 5, false, 2, 0);
      check_view(view_type(st, hs, 4), 5, false, 1, 0);
      check_view(view_type(st, hs, 5), 5, false, 0, 0);

      CHECK_THROWS_AS(n_fermion_sector_view<state_vector>(st, hs, 6),
                      std::runtime_error);
    }

    hs.add(make_space_boson(2, 5));
    hs.add(make_space_boson(3, 6));

    SECTION("Fermions and bosons") {
      check_view(view_type(st, hs, 0), 5, true, 0, 5);
      check_view(view_type(st, hs, 1), 5, true, 1, 5);
      check_view(view_type(st, hs, 2), 5, true, 2, 5);
      check_view(view_type(st, hs, 3), 5, false, 2, 5);
      check_view(view_type(st, hs, 4), 5, false, 1, 5);
      check_view(view_type(st, hs, 5), 5, false, 0, 5);

      CHECK_THROWS_AS(n_fermion_sector_view<state_vector>(st, hs, 6),
                      std::runtime_error);
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      check_view(view_type(st, hs_b, 0), 0, true, 0, 5);
      CHECK_THROWS_AS(view_type(st, hs_b, 1), std::runtime_error);
    }
  }

  SECTION("get_element() and update_add_element()") {
    hs_type hs{make_space_boson(2, 0)};
    for(int i : {0, 1, 2, 3, 4})
      hs.add(make_space_fermion(i));

    unsigned int const N = 2;

    state_vector st(n_fermion_sector_size(hs, N));
    std::iota(st.begin(), st.end(), 0);

    SECTION("const") {
      auto view = n_fermion_sector_view<state_vector const>(st, hs, N);

      CHECK(std::is_same<element_type_t<decltype(view)>, double>::value);
      // 70 == 10 00110: 2 bosons + fermions in modes 1 and 2 (sector state 4)
      CHECK(get_element(view, 70) == 2 + (4 << 2));
    }

    SECTION("non-const") {
      auto view = n_fermion_sector_view<state_vector>(st, hs, N);

      CHECK(std::is_same<element_type_t<decltype(view)>, double>::value);
      // 70 == 10 00110: 2 bosons + fermions in modes 1 and 2 (sector state 4)
      CHECK(get_element(view, 70) == 2 + (4 << 2));

      state_vector st_ref(st.size());
      std::iota(st_ref.begin(), st_ref.end(), 0);
      st_ref[18] = 15;
      update_add_element(view, 70, -3);
      CHECK(st == st_ref);

      set_zeros(view);
      CHECK(st == state_vector(st.size(), 0));
    }
  }

  SECTION("map_index()") {
    unsigned int const M = 8;
    state_vector st{};

    hs_type hs;
    for(unsigned int i = 0; i < M; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      for(unsigned int N = 0; N <= M; ++N) {
        auto view = view_type(st, hs, N);

        std::vector<sv_index_type> mapped_indices;
        for(sv_index_type index = 0; index < hs.dim(); ++index) {
          if(count_fermions(index, M) == N) {
            mapped_indices.push_back(view.map_index(index));
          }
        }
        std::sort(mapped_indices.begin(), mapped_indices.end());
        std::vector<sv_index_type> mapped_indices_ref(
            n_fermion_sector_size(hs, N));
        std::iota(mapped_indices_ref.begin(), mapped_indices_ref.end(), 0);
        CHECK(mapped_indices == mapped_indices_ref);
      }
    }

    hs.add(make_space_boson(2, int(M)));
    hs.add(make_space_boson(2, int(M + 1)));

    SECTION("Fermions and bosons") {
      for(unsigned int N = 0; N <= M; ++N) {
        auto view = view_type(st, hs, N);

        std::vector<sv_index_type> mapped_indices;
        for(sv_index_type index = 0; index < hs.dim(); ++index) {
          if(count_fermions(index, M) == N) {
            mapped_indices.push_back(view.map_index(index));
          }
        }
        std::sort(mapped_indices.begin(), mapped_indices.end());
        std::vector<sv_index_type> mapped_indices_ref(
            n_fermion_sector_size(hs, N));
        std::iota(mapped_indices_ref.begin(), mapped_indices_ref.end(), 0);
        CHECK(mapped_indices == mapped_indices_ref);
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
    // foreach() implements a mapping of the sector basis states to the full
    // Hilbert space basis states. Here we check that the mapping is the inverse
    // of map_index().

    unsigned int const M = 8;

    hs_type hs;

    for(unsigned int i = 0; i < M; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      for(unsigned int N = 0; N <= M; ++N) {
        state_vector st(n_fermion_sector_size(hs, N));
        auto view = view_type(st, hs, N);

        for(sv_index_type index = 0; index < hs.dim(); ++index) {
          if(count_fermions(index, M) == N) {
            // Add 1 here so that none of st's elements is zero
            st[view.map_index(index)] = static_cast<double>(index) + 1;
          }
        }
        foreach(view, [](sv_index_type index, double a) {
          CHECK(sv_index_type(a) == index + 1);
        });
      }
    }

    hs.add(make_space_boson(2, int(M)));
    hs.add(make_space_boson(2, int(M + 1)));

    SECTION("Fermions and bosons") {
      for(unsigned int N = 0; N <= M; ++N) {
        state_vector st(n_fermion_sector_size(hs, N));
        auto view = view_type(st, hs, N);

        for(sv_index_type index = 0; index < hs.dim(); ++index) {
          if(count_fermions(index, M) == N) {
            // Add 1 here so that none of st's elements is zero
            st[view.map_index(index)] = static_cast<double>(index) + 1;
          }
        }
        foreach(view, [](sv_index_type index, double a) {
          CHECK(sv_index_type(a) == index + 1);
        });
      }
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      state_vector st(n_fermion_sector_size(hs_b, 0));
      auto view = view_type(st, hs_b, 0);

      for(sv_index_type index = 0; index < hs_b.dim(); ++index) {
        // Add 1 here so that none of st's elements is zero
        st[view.map_index(index)] = static_cast<double>(index) + 1;
      }
      foreach(view, [](sv_index_type index, double a) {
        CHECK(sv_index_type(a) == index + 1);
      });
    }
  }

  SECTION("n_fermion_sector_basis_states()") {
    unsigned int const M = 8;
    state_vector st;

    hs_type hs;

    CHECK_THROWS_AS(n_fermion_sector_basis_states(hs, 0), std::runtime_error);
    CHECK_THROWS_AS(n_fermion_sector_basis_states(hs, 1), std::runtime_error);

    auto build_basis_states_ref = [&st](hs_type const& hs, unsigned int N) {
      auto view = view_type(st, hs, N);
      std::vector<sv_index_type> basis_states(n_fermion_sector_size(hs, N));
      for(sv_index_type index = 0; index < hs.dim(); ++index) {
        if(count_fermions(index, M) == N) {
          basis_states[view.map_index(index)] = index;
        }
      }
      return basis_states;
    };

    for(unsigned int i = 0; i < M; ++i)
      hs.add(make_space_fermion(int(i)));

    SECTION("Purely fermionic Hilbert spaces") {
      for(unsigned int N = 0; N <= M; ++N) {
        auto basis_states_ref = build_basis_states_ref(hs, N);
        auto basis_states = n_fermion_sector_basis_states(hs, N);
        CHECK(basis_states == basis_states_ref);
      }
    }

    hs.add(make_space_boson(2, int(M)));
    hs.add(make_space_boson(2, int(M + 1)));

    SECTION("Fermions and bosons") {
      for(unsigned int N = 0; N <= M; ++N) {
        auto basis_states_ref = build_basis_states_ref(hs, N);
        auto basis_states = n_fermion_sector_basis_states(hs, N);
        CHECK(basis_states == basis_states_ref);
      }
    }

    SECTION("Purely bosonic Hilbert space") {
      hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));

      CHECK_THROWS_AS(n_fermion_sector_basis_states(hs_b, 1),
                      std::runtime_error);

      std::vector<sv_index_type> basis_states_ref(
          n_fermion_sector_size(hs_b, 0));
      std::iota(basis_states_ref.begin(), basis_states_ref.end(), 0);

      auto basis_states = n_fermion_sector_basis_states(hs_b, 0);
      CHECK(basis_states == basis_states_ref);
    }
  }

  SECTION("make_nfs_view() and make_const_nfs_view()") {
    hs_type hs{make_space_fermion(0), make_space_fermion(1)};

    state_vector st;

    state_vector& st_ref = st;
    state_vector const& st_cref = st;
    state_vector&& st_rref = state_vector{1, 2, 3, 4, 5, 6};

    auto view_st = make_nfs_view(st, hs, 1);
    auto view_st_ref = make_nfs_view(st_ref, hs, 1);
    auto view_tmp = make_nfs_view(state_vector{1, 2, 3, 4, 5, 6}, hs, 1);
    auto view_rref = make_nfs_view(std::move(st_rref), hs, 1);

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
    auto cview_rref = make_const_nfs_view(std::move(st_rref), hs, 1);

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
    hs_type hs{make_space_boson(2, 0)};
    for(int i : {0, 1, 2, 3, 4})
      hs.add(make_space_fermion(i));

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
