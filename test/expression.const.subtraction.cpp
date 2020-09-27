/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <string>

using namespace libcommute;
using namespace static_indices;

TEST_CASE("Compound assignment/subtraction", "[minus_assign]") {
  SECTION("double") {
    auto expr_r = c_dag(1, "up");
    using ref_t = decltype(expr_r);
    expr_r -= 4.0;
    CHECK_THAT(expr_r, Prints<ref_t>("-4 + 1*C+(1,up)"));
    expr_r -= 0.0;
    CHECK_THAT(expr_r, Prints<ref_t>("-4 + 1*C+(1,up)"));
    expr_r -= -4.0;
    CHECK_THAT(expr_r, Prints<ref_t>("1*C+(1,up)"));
  }
  SECTION("complex from double") {
    auto expr_c = make_complex(c_dag(1, "up"));
    using ref_t = decltype(expr_c);

    expr_c -= 4.0;
    CHECK_THAT(expr_c, Prints<ref_t>("(-4,0) + (1,0)*C+(1,up)"));
    expr_c -= 0.0;
    CHECK_THAT(expr_c, Prints<ref_t>("(-4,0) + (1,0)*C+(1,up)"));
    expr_c -= -4.0;
    CHECK_THAT(expr_c, Prints<ref_t>("(1,0)*C+(1,up)"));
  }
  SECTION("my_complex") {
    auto expr = c_dag<my_complex>(1, "up");
    using ref_t = decltype(expr);

    expr -= 4;
    CHECK_THAT(expr, Prints<ref_t>("{-4,0} + {1,0}*C+(1,up)"));
    expr -= 0;
    CHECK_THAT(expr, Prints<ref_t>("{-4,0} + {1,0}*C+(1,up)"));
    expr -= -4;
    CHECK_THAT(expr, Prints<ref_t>("{1,0}*C+(1,up)"));
  }
}

TEST_CASE("Subtraction", "[minus]") {
  SECTION("double") {
    auto expr_r = c_dag(1, "up");
    using ref_t = decltype(expr_r);

    // Result type
    CHECK(std::is_same<decltype(expr_r - 2), ref_t>::value);
    CHECK(std::is_same<decltype(2 - expr_r), ref_t>::value);

    CHECK_THAT((expr_r - 0), Prints<ref_t>("1*C+(1,up)"));
    CHECK_THAT((0 - expr_r), Prints<ref_t>("-1*C+(1,up)"));
    CHECK_THAT((expr_r - 2), Prints<ref_t>("-2 + 1*C+(1,up)"));
    CHECK_THAT((2 - expr_r), Prints<ref_t>("2 + -1*C+(1,up)"));

    expr_r -= 2.0;

    CHECK_THAT((expr_r - 0), Prints<ref_t>("-2 + 1*C+(1,up)"));
    CHECK_THAT((0 - expr_r), Prints<ref_t>("2 + -1*C+(1,up)"));
    CHECK_THAT((expr_r - 2), Prints<ref_t>("-4 + 1*C+(1,up)"));
    CHECK_THAT((2 - expr_r), Prints<ref_t>("4 + -1*C+(1,up)"));
    CHECK_THAT((expr_r - (-2)), Prints<ref_t>("1*C+(1,up)"));
    CHECK_THAT(((-2) - expr_r), Prints<ref_t>("-1*C+(1,up)"));
  }
  SECTION("complex and double") {
    auto expr_c = make_complex(c_dag(1, "up"));
    using ref_t = decltype(expr_c);

    // Result type
    CHECK(std::is_same<decltype(expr_c - 2.0), ref_t>::value);
    CHECK(std::is_same<decltype(2.0 - expr_c), ref_t>::value);

    CHECK_THAT((expr_c - 0.0), Prints<ref_t>("(1,0)*C+(1,up)"));
    CHECK_THAT((0.0 - expr_c), Prints<ref_t>("(-1,-0)*C+(1,up)"));
    CHECK_THAT((expr_c - 2.0), Prints<ref_t>("(-2,0) + (1,0)*C+(1,up)"));
    CHECK_THAT((2.0 - expr_c), Prints<ref_t>("(2,0) + (-1,-0)*C+(1,up)"));

    expr_c -= 2.0;

    CHECK_THAT((expr_c - 0.0), Prints<ref_t>("(-2,0) + (1,0)*C+(1,up)"));
    CHECK_THAT((0.0 - expr_c), Prints<ref_t>("(2,-0) + (-1,-0)*C+(1,up)"));
    CHECK_THAT((expr_c - 2.0), Prints<ref_t>("(-4,0) + (1,0)*C+(1,up)"));
    CHECK_THAT((2.0 - expr_c), Prints<ref_t>("(4,-0) + (-1,-0)*C+(1,up)"));
    CHECK_THAT((expr_c - (-2.0)), Prints<ref_t>("(1,0)*C+(1,up)"));
    CHECK_THAT(((-2.0) - expr_c), Prints<ref_t>("(-1,-0)*C+(1,up)"));
  }
  SECTION("double and complex") {
    auto expr_r = c_dag(1, "up");
    using ref_t = expression<std::complex<double>, int, std::string>;
    const std::complex<double> Id(1,0);
    const std::complex<double> I(0,1);

    // Result type
    CHECK(std::is_same<decltype(expr_r - 2.0*I), ref_t>::value);
    CHECK(std::is_same<decltype(2.0*I - expr_r), ref_t>::value);

    CHECK_THAT((expr_r - 0.0*I), Prints<ref_t>("(1,0)*C+(1,up)"));
    CHECK_THAT((0.0*I - expr_r), Prints<ref_t>("(-1,0)*C+(1,up)"));
    CHECK_THAT((expr_r - 2.0*I), Prints<ref_t>("(0,-2) + (1,0)*C+(1,up)"));
    CHECK_THAT((2.0*I - expr_r), Prints<ref_t>("(0,2) + (-1,0)*C+(1,up)"));

    expr_r -= 2.0;

    CHECK_THAT((expr_r - 0.0*I), Prints<ref_t>("(-2,0) + (1,0)*C+(1,up)"));
    CHECK_THAT((0.0*I - expr_r), Prints<ref_t>("(2,0) + (-1,0)*C+(1,up)"));
    CHECK_THAT((expr_r - 2.0*Id), Prints<ref_t>("(-4,0) + (1,0)*C+(1,up)"));
    CHECK_THAT((2.0*Id - expr_r), Prints<ref_t>("(4,0) + (-1,0)*C+(1,up)"));
    CHECK_THAT((expr_r - (-2.0*Id)), Prints<ref_t>("(1,0)*C+(1,up)"));
    CHECK_THAT(((-2.0*Id) - expr_r), Prints<ref_t>("(-1,0)*C+(1,up)"));
  }
  SECTION("my_complex") {
    auto expr = c_dag<my_complex>(1, "up");
    using ref_t = decltype(expr);

    // Result type
    CHECK(std::is_same<decltype(expr - 2), ref_t>::value);
    CHECK(std::is_same<decltype(2 - expr), ref_t>::value);

    CHECK_THAT((expr - 0), Prints<ref_t>("{1,0}*C+(1,up)"));
    CHECK_THAT((0 - expr), Prints<ref_t>("{-1,0}*C+(1,up)"));
    CHECK_THAT((expr - 2), Prints<ref_t>("{-2,0} + {1,0}*C+(1,up)"));
    CHECK_THAT((2 - expr), Prints<ref_t>("{2,0} + {-1,0}*C+(1,up)"));

    expr -= 2;

    CHECK_THAT((expr - 0), Prints<ref_t>("{-2,0} + {1,0}*C+(1,up)"));
    CHECK_THAT((0 - expr), Prints<ref_t>("{2,0} + {-1,0}*C+(1,up)"));
    CHECK_THAT((expr - 2), Prints<ref_t>("{-4,0} + {1,0}*C+(1,up)"));
    CHECK_THAT((2 - expr), Prints<ref_t>("{4,0} + {-1,0}*C+(1,up)"));
    CHECK_THAT((expr - (-2)), Prints<ref_t>("{1,0}*C+(1,up)"));
    CHECK_THAT(((-2) - expr), Prints<ref_t>("{-1,0}*C+(1,up)"));
  }
}
