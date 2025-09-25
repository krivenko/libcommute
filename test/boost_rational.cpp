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

#include "commutators.hpp"

#include <type_traits>

using namespace libcommute;

TEST_CASE("Scalar traits of boost::rational", "[boost_rational]") {

  using rational_t = boost::rational<int>;

  CHECK(scalar_traits<rational_t>::is_zero(rational_t(0)));
  CHECK_FALSE(scalar_traits<rational_t>::is_zero(rational_t(3, 2)));
  CHECK(scalar_traits<rational_t>::make_const(0) == 0);
  CHECK(scalar_traits<rational_t>::make_const(1) == 1);
  CHECK(scalar_traits<rational_t>::make_const(3.0) == 3);
  CHECK(scalar_traits<rational_t>::make_const(var_number(2)) == 2);
  CHECK(scalar_traits<rational_t>::make_const(var_number(3, 4)) ==
        rational_t(3, 4));
  CHECK(scalar_traits<rational_t>::conj(rational_t(4)) == rational_t(4));

  SECTION("Result types of arithmetic operations") {
    CHECK(std::is_same<uminus_res_t<rational_t>, rational_t>::value);
    CHECK(std::is_same<sum_res_t<rational_t, rational_t>, rational_t>::value);
    CHECK(std::is_same<diff_res_t<rational_t, rational_t>, rational_t>::value);
    CHECK(std::is_same<mul_res_t<rational_t, rational_t>, rational_t>::value);
  }

  SECTION("Detect availability of compound assignments") {
    CHECK(has_add_assign<rational_t, rational_t>::value);
    CHECK(has_sub_assign<rational_t, rational_t>::value);
    CHECK(has_mul_assign<rational_t, rational_t>::value);
  }
}

TEST_CASE("Commutators and anticommutators with boost::rational",
          "[boost_rational_commutators]") {

  make_commutators_test_case_fermions<boost::rational<int>>();
  make_commutators_test_case_bosons<boost::rational<int>>();
  make_commutators_test_case_spins<boost::rational<int>>();
}
