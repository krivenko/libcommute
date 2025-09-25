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

#include <libcommute/scalars/gmpxx.hpp>

#include "commutators.hpp"

#include <type_traits>

using namespace libcommute;

TEST_CASE("Scalar traits of GMP C++ classes", "[gmpxx_classes]") {

  SECTION("mpz_class") {
    CHECK(scalar_traits<mpz_class>::is_zero(mpz_class(0)));
    CHECK_FALSE(scalar_traits<mpz_class>::is_zero(mpz_class(4)));
    CHECK(scalar_traits<mpz_class>::make_const(0) == mpz_class(0));
    CHECK(scalar_traits<mpz_class>::make_const(1) == mpz_class(1));
    CHECK(scalar_traits<mpz_class>::make_const(2.0) == mpz_class(2));
    CHECK(scalar_traits<mpz_class>::make_const(var_number(2)) == mpz_class(2));
    CHECK(scalar_traits<mpz_class>::conj(mpz_class(4)) == mpz_class(4));

    SECTION("Result types of arithmetic operations") {
      CHECK(std::is_same<uminus_res_t<mpz_class>, mpz_class>::value);
      CHECK(std::is_same<sum_res_t<mpz_class, mpz_class>, mpz_class>::value);
      CHECK(std::is_same<diff_res_t<mpz_class, mpz_class>, mpz_class>::value);
      CHECK(std::is_same<mul_res_t<mpz_class, mpz_class>, mpz_class>::value);
    }

    SECTION("Detect availability of compound assignments") {
      CHECK(has_add_assign<mpz_class, mpz_class>::value);
      CHECK(has_sub_assign<mpz_class, mpz_class>::value);
      CHECK(has_mul_assign<mpz_class, mpz_class>::value);
    }

    SECTION("Commutation relations") {
      make_commutators_test_case_fermions<mpz_class>();
      make_commutators_test_case_bosons<mpz_class>();
    }
  }

  SECTION("mpq_class") {
    CHECK(scalar_traits<mpq_class>::is_zero(mpq_class(0)));
    CHECK_FALSE(scalar_traits<mpq_class>::is_zero(mpq_class(3, 4)));
    CHECK(scalar_traits<mpq_class>::make_const(0) == mpq_class(0));
    CHECK(scalar_traits<mpq_class>::make_const(1) == mpq_class(1));
    CHECK(scalar_traits<mpq_class>::make_const(2.0) == mpq_class(2));
    CHECK(scalar_traits<mpq_class>::make_const(var_number(2, 3)) ==
          mpq_class(2, 3));
    CHECK(scalar_traits<mpq_class>::conj(mpq_class(4)) == mpq_class(4));

    SECTION("Result types of arithmetic operations") {
      CHECK(std::is_same<uminus_res_t<mpq_class>, mpq_class>::value);
      CHECK(std::is_same<sum_res_t<mpq_class, mpq_class>, mpq_class>::value);
      CHECK(std::is_same<diff_res_t<mpq_class, mpq_class>, mpq_class>::value);
      CHECK(std::is_same<mul_res_t<mpq_class, mpq_class>, mpq_class>::value);
    }

    SECTION("Detect availability of compound assignments") {
      CHECK(has_add_assign<mpq_class, mpq_class>::value);
      CHECK(has_sub_assign<mpq_class, mpq_class>::value);
      CHECK(has_mul_assign<mpq_class, mpq_class>::value);
    }

    SECTION("Commutation relations") {
      make_commutators_test_case_fermions<mpq_class>();
      make_commutators_test_case_bosons<mpq_class>();
      make_commutators_test_case_spins<mpq_class>();
    }
  }

  SECTION("mpf_class") {
    CHECK(scalar_traits<mpf_class>::is_zero(mpf_class(0)));
    CHECK_FALSE(scalar_traits<mpf_class>::is_zero(mpf_class(4)));
    CHECK(scalar_traits<mpf_class>::make_const(0) == mpf_class(0));
    CHECK(scalar_traits<mpf_class>::make_const(1.3) == mpf_class(1.3));
    CHECK(scalar_traits<mpf_class>::make_const(2.0) == mpf_class(2.0));
    CHECK(scalar_traits<mpf_class>::make_const(var_number(2)) ==
          mpf_class(2.0));
    CHECK(scalar_traits<mpf_class>::conj(mpf_class(4)) == mpf_class(4));

    SECTION("Result types of arithmetic operations") {
      CHECK(std::is_same<uminus_res_t<mpf_class>, mpf_class>::value);
      CHECK(std::is_same<sum_res_t<mpf_class, mpf_class>, mpf_class>::value);
      CHECK(std::is_same<diff_res_t<mpf_class, mpf_class>, mpf_class>::value);
      CHECK(std::is_same<mul_res_t<mpf_class, mpf_class>, mpf_class>::value);
    }

    SECTION("Detect availability of compound assignments") {
      CHECK(has_add_assign<mpf_class, mpf_class>::value);
      CHECK(has_sub_assign<mpf_class, mpf_class>::value);
      CHECK(has_mul_assign<mpf_class, mpf_class>::value);
    }

    SECTION("Commutation relations") {
      make_commutators_test_case_fermions<mpf_class>();
      make_commutators_test_case_bosons<mpf_class>();
      make_commutators_test_case_spins<mpf_class>();
    }
  }
}
