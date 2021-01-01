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

#include <libcommute/loperator/loperator.hpp>
#include <libcommute/expression/factories.hpp>

#include "./my_complex.hpp"

#include <complex>
#include <vector>

using namespace libcommute;
using namespace static_indices;

TEST_CASE("Linear operator with constant coefficients",
          "[loperator]") {

  SECTION("empty") {
    expression<double, int> expr0;
    auto hs = make_hilbert_space(expr0);
    auto lop = make_loperator(expr0, hs);

    const std::vector<double> in;
    CHECK(lop(in) == in);
    CHECK(lop * in == in);
    std::vector<double> out;
    lop(in, out);
    CHECK(out == in);
  }

  SECTION("double") {
    auto expr1 = 3*c_dag("dn");
    auto expr2 = 3*c("up");
    auto expr = 2*expr1 - 2*expr2;

    auto hs = make_hilbert_space(expr);
    auto lop1 = make_loperator(expr1, hs);
    auto lop2 = make_loperator(expr2, hs);
    auto lop = make_loperator(expr, hs);

    using state_vector = std::vector<double>;

    const state_vector in{1, 1, 1, 1};
    state_vector out(4);

    CHECK(lop1(in) == state_vector{0, 3, 0, 3});
    CHECK(lop1 * in == state_vector{0, 3, 0, 3});
    lop1(in, out);
    CHECK(out == state_vector{0, 3, 0, 3});

    CHECK(lop2(in) == state_vector{3, -3, 0, 0});
    CHECK(lop2 * in == state_vector{3, -3, 0, 0});
    lop2(in, out);
    CHECK(out == state_vector{3, -3, 0, 0});

    CHECK(lop(in) == state_vector{-6, 12, 0, 6});
    CHECK(lop * in == state_vector{-6, 12, 0, 6});
    lop(in, out);
    CHECK(out == state_vector{-6, 12, 0, 6});

    SECTION("complex state vector") {
      using state_vector = std::vector<std::complex<double>>;

      const state_vector in{1, 1, 1, 1};
      state_vector out(4);

      CHECK(lop(in) == state_vector{-6, 12, 0, 6});
      CHECK(lop *in == state_vector{-6, 12, 0, 6});
      lop(in, out);
      CHECK(out == state_vector{-6, 12, 0, 6});
    }
  }

  SECTION("complex") {
    const std::complex<double> I(0,1);

    auto expr1 = make_complex(3.0*c_dag("dn"));
    auto expr2 = make_complex(3.0*c("up"));
    auto expr = 2.0*I*expr1 - 2.0*I*expr2;

    auto hs = make_hilbert_space(expr);
    auto lop1 = make_loperator(expr1, hs);
    auto lop2 = make_loperator(expr2, hs);
    auto lop = make_loperator(expr, hs);

    using state_vector = std::vector<std::complex<double>>;

    const state_vector in{1, 1, 1, 1};
    state_vector out(4);

    CHECK(lop1(in) == state_vector{0, 3, 0, 3});
    CHECK(lop1 * in == state_vector{0, 3, 0, 3});
    lop1(in, out);
    CHECK(out == state_vector{0, 3, 0, 3});

    CHECK(lop2(in) == state_vector{3, -3, 0, 0});
    CHECK(lop2 *in == state_vector{3, -3, 0, 0});
    lop2(in, out);
    CHECK(out == state_vector{3, -3, 0, 0});

    CHECK(lop(in) == state_vector{-6.0*I, 12.0*I, .0, 6.0*I});
    CHECK(lop *in == state_vector{-6.0*I, 12.0*I, .0, 6.0*I});
    lop(in, out);
    CHECK(out == state_vector{-6.0*I, 12.0*I, .0, 6.0*I});
  }

  SECTION("my_complex") {
    const my_complex I(0, 1);

    auto expr1 = 3*c_dag<my_complex>("dn");
    auto expr2 = 3*c<my_complex>("up");
    auto expr = 2*I*expr1 - 2*I*expr2;

    auto hs = make_hilbert_space(expr);
    auto lop1 = make_loperator(expr1, hs);
    auto lop2 = make_loperator(expr2, hs);
    auto lop = make_loperator(expr, hs);

    using state_vector = std::vector<my_complex>;

    const state_vector in{1, 1, 1, 1};
    state_vector out(4, 0);

    lop1(in, out);
    CHECK(out == state_vector{0, 3, 0, 3});
    lop2(in, out);
    CHECK(out == state_vector{3, -3, 0, 0});
    lop(in, out);
    CHECK(out == state_vector{-6*I, 12*I, 0, 6*I});
  }
}

TEST_CASE("Linear operator with parameter-dependent coefficients",
          "[parametric_loperator]") {

  auto expr1 = 3*c_dag<my_complex>("dn");
  auto expr2 = 3*c<my_complex>("up");
  auto expr = 2*expr1 - 2*expr2;

  auto hs = make_hilbert_space(expr);
  auto lop1 = make_param_loperator(expr1, hs);
  auto lop2 = make_param_loperator(expr2, hs);
  auto lop = make_param_loperator(expr, hs);

  using state_vector = std::vector<my_complex>;

  SECTION("1 argument") {
    const state_vector in{1, 1, 1, 1};
    state_vector out(4, 0);

    lop1(in, out, 5);
    CHECK(out == state_vector{0, 15, 0, 15});
    lop2(in, out, 5);
    CHECK(out == state_vector{15, -15, 0, 0});
    lop(in, out, 5);
    CHECK(out == state_vector{-30, 60, 0, 30});

    SECTION("act_and_store_coeffs()") {
      std::vector<float> coeffs1, coeffs2, coeffs;
      coeffs1.reserve(1);
      coeffs2.reserve(1);
      coeffs.reserve(2);

      lop1.act_and_store_coeffs(in, out, coeffs1, 5);
      CHECK(coeffs1 == std::vector<float>{15});
      CHECK(out == state_vector{0, 15, 0, 15});
      lop2.act_and_store_coeffs(in, out, coeffs2, 5);
      CHECK(coeffs2 == std::vector<float>{15});
      CHECK(out == state_vector{15, -15, 0, 0});
      lop.act_and_store_coeffs(in, out, coeffs, 5);
      CHECK(coeffs == std::vector<float>{30, -30});
      CHECK(out == state_vector{-30, 60, 0, 30});
    }

    SECTION("at()") {
      auto lop1at = lop1.at(5);
      lop1at(in, out);
      CHECK(out == state_vector{0, 15, 0, 15});
      auto lop2at = lop2.at(5);
      lop2at(in, out);
      CHECK(out == state_vector{15, -15, 0, 0});
      auto lopat = lop.at(5);
      lopat(in, out);
      CHECK(out == state_vector{-30, 60, 0, 30});
    }
  }

  SECTION("2 arguments") {
    const state_vector in{1, 1, 1, 1};
    state_vector out(4, 0);

    lop1(in, out, 5, 5);
    CHECK(out == state_vector{0, 30, 0, 30});
    lop2(in, out, 5, 5);
    CHECK(out == state_vector{30, -30, 0, 0});
    lop(in, out, 5, 5);
    CHECK(out == state_vector{-60, 120, 0, 60});

    SECTION("act_and_store_coeffs()") {
      std::vector<my_complex> coeffs1, coeffs2, coeffs;
      coeffs1.reserve(1);
      coeffs2.reserve(1);
      coeffs.reserve(2);

      lop1.act_and_store_coeffs(in, out, coeffs1, 5, 5);
      CHECK(coeffs1 == std::vector<my_complex>{{30, 0}});
      CHECK(out == state_vector{0, 30, 0, 30});
      lop2.act_and_store_coeffs(in, out, coeffs2, 5, 5);
      CHECK(coeffs2 == std::vector<my_complex>{{30, 0}});
      CHECK(out == state_vector{30, -30, 0, 0});
      lop.act_and_store_coeffs(in, out, coeffs, 5, 5);
      CHECK(coeffs == std::vector<my_complex>{{60, 0}, {-60, 0}});
      CHECK(out == state_vector{-60, 120, 0, 60});
    }

    SECTION("at()") {
      auto lop1at = lop1.at(5, 5);
      lop1at(in, out);
      CHECK(out == state_vector{0, 30, 0, 30});
      auto lop2at = lop2.at(5, 5);
      lop2at(in, out);
      CHECK(out == state_vector{30, -30, 0, 0});
      auto lopat = lop.at(5, 5);
      lopat(in, out);
      CHECK(out == state_vector{-60, 120, 0, 60});
    }
  }
}
