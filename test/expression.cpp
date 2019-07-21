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
using namespace static_indices;

TEST_CASE("Expression with static indices", "[expression]") {
  SECTION("Constructors") {
    expr_real<int, std::string> expr0;
    CHECK_THAT(expr0, Prints<decltype(expr0)>("0"));

    expr_real<int, std::string> expr_tiny_const(1e-100);
    CHECK_THAT(expr_tiny_const, Prints<decltype(expr_tiny_const)>("0"));

    expr_real<int, std::string> expr_const(2);
    CHECK_THAT(expr_const, Prints<decltype(expr_const)>("2"));

    expr_complex<int, std::string> expr_c_const(2);
    CHECK_THAT(expr_c_const, Prints<decltype(expr_c_const)>("(2,0)"));

    expr_complex<int, std::string> expr_c_from_r(expr_const);
    CHECK_THAT(expr_c_from_r, Prints<decltype(expr_c_from_r)>("(2,0)"));

    monomial<int, std::string> mon(make_fermion(true, 1, "up"),
                                   make_fermion(false, 2, "dn"));

    expr_real<int, std::string> expr_tiny_monomial(1e-100, mon);
    CHECK_THAT(expr_tiny_monomial, Prints<decltype(expr_tiny_monomial)>("0"));

    expr_real<int, std::string> expr_monomial(3, mon);
    CHECK_THAT(expr_monomial,
               Prints<decltype(expr_monomial)>("3*C+(1,up)C(2,dn)"));
  }

  SECTION("Assignment") {
    expr_real<int, std::string> expr_const(2);
    expr_real<int, std::string> expr0;
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

  SECTION("const_iterator") {
    using expr_type = expression<double, int, std::string>;
    using mon_type = expr_type::monomial_t;
    using real::c_dag;
    using real::c;
    using real::a_dag;
    using real::a;

    expr_type expr0;
    CHECK(expr0.begin() == expr0.end());
    CHECK_FALSE(expr0.begin() != expr0.end());

    auto expr = 4.0 * c_dag(1, "up") * c(2, "dn") + 1.0 +
                3.0*a(0, "x") + 2.0*a_dag(0, "y");

    CHECK_FALSE(expr.begin() == expr.end());
    CHECK(expr.begin() != expr0.end());

    std::vector<mon_type> ref_mons = {
      mon_type(),
      mon_type(make_boson(true, 0, "y")),
      mon_type(make_boson(false, 0, "x")),
      mon_type(make_fermion(true, 1, "up"), make_fermion(false, 2, "dn"))
    };
    std::vector<double> ref_coeffs = {1.0, 2.0, 3.0, 4.0};

    int n = 0;
    expr_type::const_iterator it;

    SECTION("Prefix increment/decrement") {
      // Forward iteration
      for(it = expr.begin(); it != expr.end(); ++it, ++n) {
        auto val = *it;
        CHECK(val.monomial == ref_mons[n]);
        CHECK(val.coeff == ref_coeffs[n]);
        CHECK(it->monomial == ref_mons[n]);
        CHECK(it->coeff == ref_coeffs[n]);
      }
      // Backward iteration
      for(--it, --n; n >= 0; --it, --n) {
        auto val = *it;
        CHECK(val.monomial == ref_mons[n]);
        CHECK(val.coeff == ref_coeffs[n]);
        CHECK(it->monomial == ref_mons[n]);
        CHECK(it->coeff == ref_coeffs[n]);
      }
    }

    SECTION("Postfix increment/decrement") {
      // Forward iteration
      for(it = expr.begin(); it != expr.end(); it++, n++) {
        auto val = *it;
        CHECK(val.monomial == ref_mons[n]);
        CHECK(val.coeff == ref_coeffs[n]);
        CHECK(it->monomial == ref_mons[n]);
        CHECK(it->coeff == ref_coeffs[n]);
      }
      // Backward iteration
      for(it--, n--; n >= 0; it--, n--) {
        auto val = *it;
        CHECK(val.monomial == ref_mons[n]);
        CHECK(val.coeff == ref_coeffs[n]);
        CHECK(it->monomial == ref_mons[n]);
        CHECK(it->coeff == ref_coeffs[n]);
      }
    }

    using std::swap;
    auto it1 = expr.cbegin();
    auto it2 = expr.cend();
    swap(it1, it2);
    CHECK(it1 == expr.cend());
    CHECK(it2 == expr.cbegin());
  }

  SECTION("transform()") {
    using namespace real;
    auto expr = 4.0 * c_dag(1, "up") * c(2, "dn") + 1.0 +
                3.0 * a(0, "x") + 2.0 * a_dag(0, "y");
    using mon_type = decltype(expr)::monomial_t;

    // Multiply coefficients in front of bosonic operators by 2*I
    auto f = [](mon_type const& m, double c) -> std::complex<double> {
      return (m.size() > 0 && is_boson(m[0])) ?
             std::complex<double>(0, 2*c) : .0;
    };

    auto new_expr = transform(expr, f);

    std::vector<mon_type> ref_mons = {
      mon_type(make_boson(true, 0, "y")),
      mon_type(make_boson(false, 0, "x")),
    };
    std::vector<std::complex<double>> ref_coeffs =
      {std::complex<double>(0, 4.0), std::complex<double>(0, 6.0)};

    CHECK(new_expr.size() == 2);

    int n = 0;
    for(auto it = new_expr.begin(); it != new_expr.end(); ++it, ++n) {
      auto val = *it;
      CHECK(val.monomial == ref_mons[n]);
      CHECK(val.coeff == ref_coeffs[n]);
    }
  }
}
