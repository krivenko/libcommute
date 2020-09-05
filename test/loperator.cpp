/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
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

TEST_CASE("Linear operator with constant coefficients",
          "[loperator]") {

  SECTION("empty") {
    expression<double, int> expr0;
    auto hs = make_hilbert_space(expr0);
    auto qop = make_loperator(expr0, hs);

    const std::vector<double> in;
    CHECK(qop(in) == in);
    CHECK(qop * in == in);
    std::vector<double> out;
    qop(in, out);
    CHECK(out == in);
  }

  SECTION("double") {
    using namespace static_indices::real;

    auto expr1 = 3*c_dag("dn");
    auto expr2 = 3*c("up");
    auto expr = 2*expr1 - 2*expr2;

    auto hs = make_hilbert_space(expr);
    auto qop1 = make_loperator(expr1, hs);
    auto qop2 = make_loperator(expr2, hs);
    auto qop = make_loperator(expr, hs);

    using state_vector = std::vector<double>;

    const state_vector in{1, 1, 1, 1};
    state_vector out(4);

    CHECK(qop1(in) == state_vector{0, 3, 0, 3});
    CHECK(qop1 * in == state_vector{0, 3, 0, 3});
    qop1(in, out);
    CHECK(out == state_vector{0, 3, 0, 3});

    CHECK(qop2(in) == state_vector{3, -3, 0, 0});
    CHECK(qop2 * in == state_vector{3, -3, 0, 0});
    qop2(in, out);
    CHECK(out == state_vector{3, -3, 0, 0});

    CHECK(qop(in) == state_vector{-6, 12, 0, 6});
    CHECK(qop * in == state_vector{-6, 12, 0, 6});
    qop(in, out);
    CHECK(out == state_vector{-6, 12, 0, 6});

    SECTION("complex state vector") {
      using state_vector = std::vector<std::complex<double>>;

      const state_vector in{1, 1, 1, 1};
      state_vector out(4);

      CHECK(qop(in) == state_vector{-6, 12, 0, 6});
      CHECK(qop *in == state_vector{-6, 12, 0, 6});
      qop(in, out);
      CHECK(out == state_vector{-6, 12, 0, 6});
    }
  }

  SECTION("complex") {
    using namespace static_indices::complex;

    const std::complex<double> I(0,1);

    auto expr1 = 3.0*c_dag("dn");
    auto expr2 = 3.0*c("up");
    auto expr = 2.0*I*expr1 - 2.0*I*expr2;

    auto hs = make_hilbert_space(expr);
    auto qop1 = make_loperator(expr1, hs);
    auto qop2 = make_loperator(expr2, hs);
    auto qop = make_loperator(expr, hs);

    using state_vector = std::vector<std::complex<double>>;

    const state_vector in{1, 1, 1, 1};
    state_vector out(4);

    CHECK(qop1(in) == state_vector{0, 3, 0, 3});
    CHECK(qop1 * in == state_vector{0, 3, 0, 3});
    qop1(in, out);
    CHECK(out == state_vector{0, 3, 0, 3});

    CHECK(qop2(in) == state_vector{3, -3, 0, 0});
    CHECK(qop2 *in == state_vector{3, -3, 0, 0});
    qop2(in, out);
    CHECK(out == state_vector{3, -3, 0, 0});

    CHECK(qop(in) == state_vector{-6.0*I, 12.0*I, .0, 6.0*I});
    CHECK(qop *in == state_vector{-6.0*I, 12.0*I, .0, 6.0*I});
    qop(in, out);
    CHECK(out == state_vector{-6.0*I, 12.0*I, .0, 6.0*I});
  }

  SECTION("my_complex") {
    using namespace static_indices;

    const my_complex I(0, 1);

    auto expr1 = 3*c_dag<my_complex>("dn");
    auto expr2 = 3*c<my_complex>("up");
    auto expr = 2*I*expr1 - 2*I*expr2;

    auto hs = make_hilbert_space(expr);
    auto qop1 = make_loperator(expr1, hs);
    auto qop2 = make_loperator(expr2, hs);
    auto qop = make_loperator(expr, hs);

    using state_vector = std::vector<my_complex>;

    const state_vector in{1, 1, 1, 1};
    state_vector out(4, 0);

    qop1(in, out);
    CHECK(out == state_vector{0, 3, 0, 3});
    qop2(in, out);
    CHECK(out == state_vector{3, -3, 0, 0});
    qop(in, out);
    CHECK(out == state_vector{-6*I, 12*I, 0, 6*I});
  }
}

TEST_CASE("Linear operator with parameter-dependent coefficients",
          "[parametric_loperator]") {

  using namespace static_indices;

  auto expr1 = 3*c_dag<my_complex>("dn");
  auto expr2 = 3*c<my_complex>("up");
  auto expr = 2*expr1 - 2*expr2;

  auto hs = make_hilbert_space(expr);
  auto qop1 = make_param_loperator(expr1, hs);
  auto qop2 = make_param_loperator(expr2, hs);
  auto qop = make_param_loperator(expr, hs);

  using state_vector = std::vector<my_complex>;

  SECTION("1 argument") {
    const state_vector in{1, 1, 1, 1};
    state_vector out(4, 0);

    qop1(in, out, 5);
    CHECK(out == state_vector{0, 15, 0, 15});
    qop2(in, out, 5);
    CHECK(out == state_vector{15, -15, 0, 0});
    qop(in, out, 5);
    CHECK(out == state_vector{-30, 60, 0, 30});

    SECTION("act_and_store_coeffs()") {
      std::vector<float> coeffs1, coeffs2, coeffs;
      coeffs1.reserve(1);
      coeffs2.reserve(1);
      coeffs.reserve(2);

      qop1.act_and_store_coeffs(in, out, coeffs1, 5);
      CHECK(coeffs1 == std::vector<float>{15});
      CHECK(out == state_vector{0, 15, 0, 15});
      qop2.act_and_store_coeffs(in, out, coeffs2, 5);
      CHECK(coeffs2 == std::vector<float>{15});
      CHECK(out == state_vector{15, -15, 0, 0});
      qop.act_and_store_coeffs(in, out, coeffs, 5);
      CHECK(coeffs == std::vector<float>{30, -30});
      CHECK(out == state_vector{-30, 60, 0, 30});
    }

    SECTION("at()") {
      auto qop1at = qop1.at(5);
      qop1at(in, out);
      CHECK(out == state_vector{0, 15, 0, 15});
      auto qop2at = qop2.at(5);
      qop2at(in, out);
      CHECK(out == state_vector{15, -15, 0, 0});
      auto qopat = qop.at(5);
      qopat(in, out);
      CHECK(out == state_vector{-30, 60, 0, 30});
    }
  }

  SECTION("2 arguments") {
    const state_vector in{1, 1, 1, 1};
    state_vector out(4, 0);

    qop1(in, out, 5, 5);
    CHECK(out == state_vector{0, 30, 0, 30});
    qop2(in, out, 5, 5);
    CHECK(out == state_vector{30, -30, 0, 0});
    qop(in, out, 5, 5);
    CHECK(out == state_vector{-60, 120, 0, 60});

    SECTION("act_and_store_coeffs()") {
      std::vector<my_complex> coeffs1, coeffs2, coeffs;
      coeffs1.reserve(1);
      coeffs2.reserve(1);
      coeffs.reserve(2);

      qop1.act_and_store_coeffs(in, out, coeffs1, 5, 5);
      CHECK(coeffs1 == std::vector<my_complex>{{30, 0}});
      CHECK(out == state_vector{0, 30, 0, 30});
      qop2.act_and_store_coeffs(in, out, coeffs2, 5, 5);
      CHECK(coeffs2 == std::vector<my_complex>{{30, 0}});
      CHECK(out == state_vector{30, -30, 0, 0});
      qop.act_and_store_coeffs(in, out, coeffs, 5, 5);
      CHECK(coeffs == std::vector<my_complex>{{60, 0}, {-60, 0}});
      CHECK(out == state_vector{-60, 120, 0, 60});
    }

    SECTION("at()") {
      auto qop1at = qop1.at(5, 5);
      qop1at(in, out);
      CHECK(out == state_vector{0, 30, 0, 30});
      auto qop2at = qop2.at(5, 5);
      qop2at(in, out);
      CHECK(out == state_vector{30, -30, 0, 0});
      auto qopat = qop.at(5, 5);
      qopat(in, out);
      CHECK(out == state_vector{-60, 120, 0, 60});
    }
  }
}