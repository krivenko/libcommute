/*******************************************************************************
 *
 * This file is part of libcommute, a C++11/14/17 header-only library allowing
 * to manipulate polynomial expressions with quantum-mechanical operators.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include "catch2/catch.hpp"

#include "int_complex.hpp"
#include "utility.hpp"

#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/expression.hpp>
#include <libcommute/expression/factories.hpp>

using namespace libcommute;

TEST_CASE("Expression with static indices", "[expression]") {
  SECTION("Constructors") {
    expression_real<int, std::string> expr0;
    CHECK_THAT(expr0, Prints<decltype(expr0)>("0"));

    expression_real<int, std::string> expr_tiny_const(1e-100);
    CHECK_THAT(expr_tiny_const, Prints<decltype(expr_tiny_const)>("0"));

    expression_real<int, std::string> expr_const(2);
    CHECK_THAT(expr_const, Prints<decltype(expr_const)>("2"));

    expression_complex<int, std::string> expr_c_const(2);
    CHECK_THAT(expr_c_const, Prints<decltype(expr_c_const)>("(2,0)"));

    expression_complex<int, std::string> expr_c_from_r(expr_const);
    CHECK_THAT(expr_c_from_r, Prints<decltype(expr_c_from_r)>("(2,0)"));

    monomial<int, std::string> mon(make_fermion(true, 1, "up"),
                                   make_fermion(false, 2, "dn"));

    expression_real<int, std::string> expr_tiny_monomial(1e-100, mon);
    CHECK_THAT(expr_tiny_monomial, Prints<decltype(expr_tiny_monomial)>("0"));

    expression_real<int, std::string> expr_monomial(3, mon);
    CHECK_THAT(expr_monomial,
               Prints<decltype(expr_monomial)>("3*C+(1,up)C(2,dn)"));
  }

  SECTION("Assignment") {
    expression_real<int, std::string> expr_const(2);
    expression_real<int, std::string> expr0;
    expr0 = expr_const;
    CHECK_THAT(expr0, Prints<decltype(expr0)>("2"));
  }

  SECTION("Unary minus") {
    auto expr_r = real::c_dag(1, "up");
    CHECK_THAT(-expr_r, Prints<decltype(expr_r)>("-1*C+(1,up)"));

    auto expr_static_int = c_dag<int_complex>(1, "up");
    CHECK_THAT(-expr_static_int,
               Prints<decltype(expr_static_int)>("{-1,0}*C+(1,up)"));
  }
}
