/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include "catch2/catch.hpp"

#include "my_complex.hpp"
#include "print_matcher.hpp"

#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/expression.hpp>
#include <libcommute/expression/factories.hpp>

using namespace libcommute;

template<typename ScalarType> void make_commutators_test_case() {

  using namespace static_indices;

  using S = ScalarType;
  using expr_t = expression<ScalarType, int>;

  auto comm = [](expr_t const& expr1, expr_t const& expr2) {
    return expr1 * expr2 - expr2 * expr1;
  };
  auto acomm = [](expr_t const& expr1, expr_t const& expr2) {
    return expr1 * expr2 + expr2 * expr1;
  };

  SECTION("Fermions") {
    const int N = 4;

    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        auto acomm_c_dag_c_dag = acomm(c_dag<S>(i), c_dag<S>(j));
        CHECK(acomm_c_dag_c_dag == expr_t());

        auto acomm_c_c = acomm(c<S>(i), c<S>(j));
        CHECK(acomm_c_c == expr_t());

        auto acomm_c_dag_c = acomm(c_dag<S>(i), c<S>(j));
        CHECK(acomm_c_dag_c == expr_t(i == j));
      }
    }
  }

  SECTION("Bosons") {
    const int N = 4;

    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        auto comm_a_dag_a_dag = comm(a_dag<S>(i), a_dag<S>(j));
        CHECK(comm_a_dag_a_dag == expr_t());

        auto comm_a_a = comm(a<S>(i), a<S>(j));
        CHECK(comm_a_a == expr_t());

        auto comm_a_a_dag = comm(a<S>(j), a_dag<S>(i));
        CHECK(comm_a_a_dag == expr_t(i == j));
      }
    }
  }

  SECTION("Spin") {
    const int N = 4;

    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        // Spin-1/2
        auto comm_Sp_Sm = comm(S_p<S>(i), S_m<S>(j));
        CHECK(comm_Sp_Sm == 2.0 * S(i == j) * S_z<S>(i));

        auto comm_Sz_Sp = comm(S_z<S>(i), S_p<S>(j));
        CHECK(comm_Sz_Sp == S(i == j) * S_p<S>(i));

        auto comm_Sz_Sm = comm(S_z<S>(i), S_m<S>(j));
        CHECK(comm_Sz_Sm == -1.0 * S(i == j) * S_m<S>(i));

        // Spin-1
        auto comm_S1p_S1m = comm(S_p<3, S>(i), S_m<3, S>(j));
        CHECK(comm_S1p_S1m == 2.0 * S(i == j) * S_z<3, S>(i));

        auto comm_S1z_S1p = comm(S_z<3, S>(i), S_p<3, S>(j));
        CHECK(comm_S1z_S1p == S(i == j) * S_p<3, S>(i));

        auto comm_S1z_S1m = comm(S_z<3, S>(i), S_m<3, S>(j));
        CHECK(comm_S1z_S1m == -1.0 * S(i == j) * S_m<3, S>(i));

        // Spin-3/2
        auto comm_S32p_S32m = comm(S_p<4, S>(i), S_m<4, S>(j));
        CHECK(comm_S32p_S32m == 2.0 * S(i == j) * S_z<4, S>(i));

        auto comm_S32z_S32p = comm(S_z<4, S>(i), S_p<4, S>(j));
        CHECK(comm_S32z_S32p == S(i == j) * S_p<4, S>(i));

        auto comm_S32z_S32m = comm(S_z<4, S>(i), S_m<4, S>(j));
        CHECK(comm_S32z_S32m == -1.0 * S(i == j) * S_m<4, S>(i));
      }
    }
  }
}

TEST_CASE("Commutation relations (double)", "[commutators_double]") {
  make_commutators_test_case<double>();
}

TEST_CASE("Commutation relations (complex)", "[commutators_complex]") {
  make_commutators_test_case<std::complex<double>>();
}

TEST_CASE("Commutation relations (my_complex)", "[commutators_my_complex]") {
  make_commutators_test_case<my_complex>();
}
