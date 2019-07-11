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

#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/expression.hpp>

using namespace libcommute;

TEST_CASE("Expression with static indices", "[expression]") {

  using monomial_t = monomial<int, std::string>;

  SECTION("Construction") {
    expression_real<int, std::string> expr0;
    CHECK(expr0.size() == 0);
    CHECK(expr0.get_monomials().empty());

    expression_real<int, std::string> expr_tiny_const(1e-100);
    CHECK(expr_tiny_const.size() == 0);
    CHECK(expr_tiny_const.get_monomials().empty());

    expression_real<int, std::string> expr_const(2);
    CHECK(expr_const.size() == 1);
    CHECK((expr_const.get_monomials().begin()->first) == monomial_t{});
    CHECK((expr_const.get_monomials().begin()->second) == 2.0);

    expression_complex<int, std::string> expr_c_const(2);
    CHECK(expr_c_const.size() == 1);
    CHECK((expr_c_const.get_monomials().begin()->first) == monomial_t{});
    CHECK((expr_c_const.get_monomials().begin()->second) == 2.0);

    expression_complex<int, std::string> expr_c_from_r(expr_const);
    CHECK(expr_c_from_r.size() == 1);
    CHECK((expr_c_from_r.get_monomials().begin()->first) == monomial_t{});
    CHECK((expr_c_from_r.get_monomials().begin()->second) == 2.0);

    monomial_t monomial(make_fermion(true, 1, "up"),
                        make_fermion(false, 2, "dn"));

    expression_real<int, std::string> expr_tiny_monomial(1e-100, monomial);
    CHECK(expr_tiny_monomial.size() == 0);
    CHECK(expr_tiny_monomial.get_monomials().empty());

    expression_real<int, std::string> expr_monomial(3, monomial);
    CHECK(expr_monomial.size() == 1);
    CHECK((expr_monomial.get_monomials().begin()->first) == monomial);
    CHECK((expr_monomial.get_monomials().begin()->second) == 3.0);
  }

  SECTION("assignment") {
    expression_real<int, std::string> expr_const(2);
    expression_real<int, std::string> expr0;
    expr0 = expr_const;
    CHECK(expr0.size() == 1);
    CHECK((expr0.get_monomials().begin()->first) == monomial_t{});
    CHECK((expr0.get_monomials().begin()->second) == 2.0);
  }

  // TODO
}
