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

#include <catch.hpp>

#include <libcommute/scalars/boost_rational.hpp>

#include <libcommute/expression/expression.hpp>
#include <libcommute/expression/factories.hpp>

#include <type_traits>

using namespace libcommute;

TEST_CASE("Scalar traits of boost::rational", "[boost_rational]") {

  using rational_t = boost::rational<int>;

  CHECK(scalar_traits<rational_t>::is_zero(rational_t(0)));
  CHECK_FALSE(scalar_traits<rational_t>::is_zero(rational_t(3, 2)));
  CHECK(scalar_traits<rational_t>::make_const(0) == 0);
  CHECK(scalar_traits<rational_t>::make_const(1) == 1);
  CHECK(scalar_traits<rational_t>::make_const(var_number(2)) == 2);
  CHECK(scalar_traits<rational_t>::make_const(var_number(3, 4)) ==
        rational_t(3, 4));
  CHECK(scalar_traits<rational_t>::conj(rational_t(4)) == rational_t(4));

  SECTION("Result types of arithmetic operations") {
    CHECK(std::is_same<minus_type<rational_t>, rational_t>::value);
    CHECK(std::is_same<sum_type<rational_t, rational_t>, rational_t>::value);
    CHECK(std::is_same<diff_type<rational_t, rational_t>, rational_t>::value);
    CHECK(std::is_same<mul_type<rational_t, rational_t>, rational_t>::value);
  }

  SECTION("Detect availability of compound assignments") {
    CHECK(has_add_assign<rational_t, rational_t>::value);
    CHECK(has_sub_assign<rational_t, rational_t>::value);
    CHECK(has_mul_assign<rational_t, rational_t>::value);
  }
}

TEST_CASE("Commutators and anticommutators with boost::rational",
          "[boost_rational_commutators]") {

  using namespace static_indices;

  using S = boost::rational<int>;
  using expr_t = expression<S, int>;

  auto comm = [](expr_t const& expr1, expr_t const& expr2) {
    return expr1 * expr2 - expr2 * expr1;
  };
  auto acomm = [](expr_t const& expr1, expr_t const& expr2) {
    return expr1 * expr2 + expr2 * expr1;
  };

  SECTION("Fermions") {
    int const N = 4;

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
    int const N = 4;

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
    int const N = 4;

    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        // Spin-1/2
        auto comm_Sp_Sm = comm(S_p<S>(i), S_m<S>(j));
        CHECK(comm_Sp_Sm == S(2 * (i == j)) * S_z<S>(i));

        auto comm_Sz_Sp = comm(S_z<S>(i), S_p<S>(j));
        CHECK(comm_Sz_Sp == S(i == j) * S_p<S>(i));

        auto comm_Sz_Sm = comm(S_z<S>(i), S_m<S>(j));
        CHECK(comm_Sz_Sm == S(-(i == j)) * S_m<S>(i));

        // Spin-1
        auto comm_S1p_S1m = comm(S_p<3, S>(i), S_m<3, S>(j));
        CHECK(comm_S1p_S1m == S(2 * (i == j)) * S_z<3, S>(i));

        auto comm_S1z_S1p = comm(S_z<3, S>(i), S_p<3, S>(j));
        CHECK(comm_S1z_S1p == S(i == j) * S_p<3, S>(i));

        auto comm_S1z_S1m = comm(S_z<3, S>(i), S_m<3, S>(j));
        CHECK(comm_S1z_S1m == S(-(i == j)) * S_m<3, S>(i));

        // Spin-3/2
        auto comm_S32p_S32m = comm(S_p<4, S>(i), S_m<4, S>(j));
        CHECK(comm_S32p_S32m == S(2 * (i == j)) * S_z<4, S>(i));

        auto comm_S32z_S32p = comm(S_z<4, S>(i), S_p<4, S>(j));
        CHECK(comm_S32z_S32p == S(i == j) * S_p<4, S>(i));

        auto comm_S32z_S32m = comm(S_z<4, S>(i), S_m<4, S>(j));
        CHECK(comm_S32z_S32m == S(-(i == j)) * S_m<4, S>(i));
      }
    }
  }
}
