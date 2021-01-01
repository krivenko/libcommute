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

#include "catch2/catch.hpp"

#include "my_complex.hpp"
#include "print_matcher.hpp"

#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/generator_boson.hpp>
#include <libcommute/expression/expression.hpp>
#include <libcommute/expression/factories.hpp>

#include <string>

using namespace libcommute;
using namespace static_indices;

TEST_CASE("Compound assignment/addition", "[plus_assign]") {
  SECTION("double") {
    auto expr_r = c_dag(1, "up");
    using ref_t = decltype(expr_r);

    expr_r += c(2, "dn");
    CHECK_THAT(expr_r, Prints<ref_t>("1*C+(1,up) + 1*C(2,dn)"));
    expr_r += ref_t();
    CHECK_THAT(expr_r, Prints<ref_t>("1*C+(1,up) + 1*C(2,dn)"));
    expr_r += -c_dag(1, "up");
    CHECK_THAT(expr_r, Prints<ref_t>("1*C(2,dn)"));
  }
  SECTION("complex from double") {
    auto expr_c = make_complex(c_dag(1, "up"));
    using ref_t = decltype(expr_c);

    expr_c += c(2, "dn");
    CHECK_THAT(expr_c, Prints<ref_t>("(1,0)*C+(1,up) + (1,0)*C(2,dn)"));
    expr_c += ref_t();
    CHECK_THAT(expr_c, Prints<ref_t>("(1,0)*C+(1,up) + (1,0)*C(2,dn)"));
    expr_c += -c_dag(1, "up");
    CHECK_THAT(expr_c, Prints<ref_t>("(1,0)*C(2,dn)"));
  }
  SECTION("my_complex") {
    auto expr = c_dag<my_complex>(1, "up");
    using ref_t = decltype(expr);

    expr += c<my_complex>(2, "dn");
    CHECK_THAT(expr, Prints<ref_t>("{1,0}*C+(1,up) + {1,0}*C(2,dn)"));
    expr += ref_t();
    CHECK_THAT(expr, Prints<ref_t>("{1,0}*C+(1,up) + {1,0}*C(2,dn)"));
    expr += -c_dag<my_complex>(1, "up");
    CHECK_THAT(expr, Prints<ref_t>("{1,0}*C(2,dn)"));
  }
}

TEST_CASE("Addition", "[plus]") {
  SECTION("double") {
    auto expr_r = c_dag(1, "up");
    using ref_t = decltype(expr_r);

    // Result type
    CHECK(std::is_same<decltype(expr_r + c(2, "dn")), ref_t>::value);
    CHECK(std::is_same<decltype(c(2, "dn") + expr_r), ref_t>::value);

    CHECK_THAT((ref_t() + ref_t()), Prints<ref_t>("0"));

    CHECK_THAT((expr_r + ref_t()), Prints<ref_t>("1*C+(1,up)"));
    CHECK_THAT((ref_t() + expr_r), Prints<ref_t>("1*C+(1,up)"));
    CHECK_THAT((expr_r + c(2, "dn")),
               Prints<ref_t>("1*C+(1,up) + 1*C(2,dn)"));
    CHECK_THAT((c(2, "dn") + expr_r),
               Prints<ref_t>("1*C+(1,up) + 1*C(2,dn)"));

    expr_r += c(2, "dn");

    CHECK_THAT((expr_r + ref_t()), Prints<ref_t>("1*C+(1,up) + 1*C(2,dn)"));
    CHECK_THAT((ref_t() + expr_r), Prints<ref_t>("1*C+(1,up) + 1*C(2,dn)"));
    CHECK_THAT((expr_r + a(0, "x")),
               Prints<ref_t>("1*C+(1,up) + 1*C(2,dn) + 1*A(0,x)"));
    CHECK_THAT((a(0, "x") + expr_r),
               Prints<ref_t>("1*C+(1,up) + 1*C(2,dn) + 1*A(0,x)"));
    CHECK_THAT((expr_r + (-c_dag(1, "up"))), Prints<ref_t>("1*C(2,dn)"));
    CHECK_THAT(((-c_dag(1, "up")) + expr_r), Prints<ref_t>("1*C(2,dn)"));

    CHECK_THAT(((c_dag(1, "up") + c(2, "dn")) + (c(2, "dn") + 2.0)),
               Prints<ref_t>("2 + 1*C+(1,up) + 2*C(2,dn)"));
  }
  SECTION("complex and double") {
    auto expr1 = make_complex(c_dag(1, "up"));
    auto expr2 = c(2, "dn");
    using ref1_t = decltype(expr1);
    using ref2_t = decltype(expr2);

    // Result type
    CHECK(std::is_same<decltype(expr1 + expr2), ref1_t>::value);
    CHECK(std::is_same<decltype(expr2 + expr1), ref1_t>::value);

    CHECK_THAT((ref1_t() + ref2_t()), Prints<ref1_t>("(0,0)"));
    CHECK_THAT((ref2_t() + ref1_t()), Prints<ref1_t>("(0,0)"));

    CHECK_THAT((expr1 + ref2_t()), Prints<ref1_t>("(1,0)*C+(1,up)"));
    CHECK_THAT((ref2_t() + expr1), Prints<ref1_t>("(1,0)*C+(1,up)"));
    CHECK_THAT((expr2 + ref1_t()), Prints<ref1_t>("(1,0)*C(2,dn)"));
    CHECK_THAT((ref1_t() + expr2), Prints<ref1_t>("(1,0)*C(2,dn)"));
    CHECK_THAT((expr1 + expr2),
               Prints<ref1_t>("(1,0)*C+(1,up) + (1,0)*C(2,dn)"));
    CHECK_THAT((expr2 + expr1),
               Prints<ref1_t>("(1,0)*C+(1,up) + (1,0)*C(2,dn)"));

    expr1 += expr2;

    CHECK_THAT((expr1 + ref2_t()),
               Prints<ref1_t>("(1,0)*C+(1,up) + (1,0)*C(2,dn)"));
    CHECK_THAT((ref2_t() + expr1),
               Prints<ref1_t>("(1,0)*C+(1,up) + (1,0)*C(2,dn)"));
    CHECK_THAT((expr1 + a(0, "x")),
               Prints<ref1_t>("(1,0)*C+(1,up) + (1,0)*C(2,dn) + (1,0)*A(0,x)"));
    CHECK_THAT((a(0, "x") + expr1),
               Prints<ref1_t>("(1,0)*C+(1,up) + (1,0)*C(2,dn) + (1,0)*A(0,x)"));
    CHECK_THAT((expr1 + (-make_complex(c_dag(1, "up")))),
               Prints<ref1_t>("(1,0)*C(2,dn)"));
    CHECK_THAT(((-make_complex(c_dag(1, "up"))) + expr1),
               Prints<ref1_t>("(1,0)*C(2,dn)"));

    CHECK_THAT((make_complex(c_dag(1, "up") + c(2, "dn")) +
               (c(2, "dn") + 2.0)),
               Prints<ref1_t>("(2,0) + (1,0)*C+(1,up) + (2,0)*C(2,dn)"));
  }
  SECTION("my_complex") {
    auto expr = c_dag<my_complex>(1, "up");
    using ref_t = decltype(expr);

    // Result type
    CHECK(std::is_same<decltype(expr + c<my_complex>(2, "dn")), ref_t>::value);
    CHECK(std::is_same<decltype(c<my_complex>(2, "dn") + expr), ref_t>::value);

    CHECK_THAT((ref_t() + ref_t()), Prints<ref_t>("{0,0}"));

    CHECK_THAT((expr + ref_t()), Prints<ref_t>("{1,0}*C+(1,up)"));
    CHECK_THAT((ref_t() + expr), Prints<ref_t>("{1,0}*C+(1,up)"));
    CHECK_THAT((expr + c<my_complex>(2, "dn")),
               Prints<ref_t>("{1,0}*C+(1,up) + {1,0}*C(2,dn)"));
    CHECK_THAT((c<my_complex>(2, "dn") + expr),
               Prints<ref_t>("{1,0}*C+(1,up) + {1,0}*C(2,dn)"));

    expr += c<my_complex>(2, "dn");

    CHECK_THAT((expr + ref_t()),
               Prints<ref_t>("{1,0}*C+(1,up) + {1,0}*C(2,dn)"));
    CHECK_THAT((ref_t() + expr),
               Prints<ref_t>("{1,0}*C+(1,up) + {1,0}*C(2,dn)"));
    CHECK_THAT((expr + a<my_complex>(0, "x")),
               Prints<ref_t>("{1,0}*C+(1,up) + {1,0}*C(2,dn) + {1,0}*A(0,x)"));
    CHECK_THAT((a<my_complex>(0, "x") + expr),
               Prints<ref_t>("{1,0}*C+(1,up) + {1,0}*C(2,dn) + {1,0}*A(0,x)"));
    CHECK_THAT((expr + (-c_dag<my_complex>(1, "up"))),
               Prints<ref_t>("{1,0}*C(2,dn)"));
    CHECK_THAT(((-c_dag<my_complex>(1, "up")) + expr),
               Prints<ref_t>("{1,0}*C(2,dn)"));

    CHECK_THAT(((c_dag<my_complex>(1, "up") + c<my_complex>(2, "dn")) +
                (c<my_complex>(2, "dn") + my_complex{2,0})),
               Prints<ref_t>("{2,0} + {1,0}*C+(1,up) + {2,0}*C(2,dn)"));
  }
}
