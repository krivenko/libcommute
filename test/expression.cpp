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

#include "static_int.hpp"

#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/expression.hpp>
#include <libcommute/expression/factories.hpp>

using namespace libcommute;

using monomial_t = monomial<int, std::string>;

template<typename E, typename S, typename... Generators>
void check_monomial(E const& expr, S ref_coeff, Generators&&... generators) {
  monomial_t ref_monomial(std::forward<Generators>(generators)...);

  CHECK(expr.get_monomials().size() == 1);
  CHECK((expr.get_monomials().begin()->first) == ref_monomial);
  CHECK((expr.get_monomials().begin()->second) == ref_coeff);
}

TEST_CASE("Expression with static indices", "[expression]") {

  using monomial_t = monomial<int, std::string>;

  SECTION("Constructors") {
    expression_real<int, std::string> expr0;
    CHECK(expr0.size() == 0);
    CHECK(expr0.get_monomials().empty());

    expression_real<int, std::string> expr_tiny_const(1e-100);
    CHECK(expr_tiny_const.size() == 0);
    CHECK(expr_tiny_const.get_monomials().empty());

    expression_real<int, std::string> expr_const(2);
    check_monomial(expr_const, 2.0);

    expression_complex<int, std::string> expr_c_const(2);
    check_monomial(expr_c_const, 2.0);

    expression_complex<int, std::string> expr_c_from_r(expr_const);
    check_monomial(expr_c_from_r, 2.0);

    monomial_t monomial(make_fermion(true, 1, "up"),
                        make_fermion(false, 2, "dn"));

    expression_real<int, std::string> expr_tiny_monomial(1e-100, monomial);
    CHECK(expr_tiny_monomial.size() == 0);
    CHECK(expr_tiny_monomial.get_monomials().empty());

    expression_real<int, std::string> expr_monomial(3, monomial);
    check_monomial(expr_monomial,
                   3.0,
                   make_fermion(true, 1, "up"),
                   make_fermion(false, 2, "dn"));
  }

  SECTION("assignment") {
    expression_real<int, std::string> expr_const(2);
    expression_real<int, std::string> expr0;
    expr0 = expr_const;
    check_monomial(expr0, 2.0);
  }

  SECTION("Unary minus") {
    auto expr_r = real::c_dag(1, "up");
    check_monomial(-expr_r, -1.0, make_fermion(true, 1, "up"));

    auto expr_static_int = c_dag<static_int<1>>(1, "up");
    check_monomial(-expr_static_int,
                   static_int<-1>(),
                   make_fermion(true, 1, "up"));
  }

  SECTION("Multiplication by scalar") {
    auto gen = make_fermion(true, 1, "up");

    // ScalarType = double
    auto expr_r = real::c_dag(1, "up");
    check_monomial(expr_r*2, 2.0, gen);
    check_monomial(2*expr_r, 2.0, gen);
    CHECK((expr_r*0).size() == 0);
    CHECK((0*expr_r).size() == 0);

    expr_r *= 2;
    check_monomial(expr_r, 2.0, gen);
    expr_r *= 0;
    CHECK(expr_r.size() == 0);

    // ScalarType = static_int<N>
    auto expr_static_int = c_dag<static_int<1>>(1, "up");
    check_monomial(expr_static_int * static_int<2>(),
                   static_int<2>(),
                   gen);
    check_monomial(static_int<2>() * expr_static_int,
                   static_int<2>(),
                   gen);
    CHECK((expr_static_int * static_int<0>()).size() == 0);
    CHECK((static_int<0>() * expr_static_int).size() == 0);

    expr_static_int *= static_int<1>();
    check_monomial(expr_static_int,
                   static_int<1>(),
                   gen);
  }
  // TODO
}
