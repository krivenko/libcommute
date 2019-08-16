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

#include <libcommute/qoperator/qoperator.hpp>
#include <libcommute/expression/factories.hpp>

#include "./int_complex.hpp"

#include <complex>
#include <vector>

using namespace libcommute;

TEST_CASE("Quantum-mechanical operator with constant coefficients",
          "[qoperator]") {

  SECTION("empty") {
    expression<double, int> expr0;
    auto hs = make_hilbert_space(expr0);
    auto qop = make_qoperator(expr0, hs);

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
    auto qop1 = make_qoperator(expr1, hs);
    auto qop2 = make_qoperator(expr2, hs);
    auto qop = make_qoperator(expr, hs);

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
    auto qop1 = make_qoperator(expr1, hs);
    auto qop2 = make_qoperator(expr2, hs);
    auto qop = make_qoperator(expr, hs);

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

  SECTION("int_complex") {
    using namespace static_indices;

    const int_complex I(0, 1);

    auto expr1 = 3*c_dag<int_complex>("dn");
    auto expr2 = 3*c<int_complex>("up");
    auto expr = 2*I*expr1 - 2*I*expr2;

    auto hs = make_hilbert_space(expr);
    auto qop1 = make_qoperator(expr1, hs);
    auto qop2 = make_qoperator(expr2, hs);
    auto qop = make_qoperator(expr, hs);

    using state_vector = std::vector<int_complex>;

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

TEST_CASE("Quantum-mechanical operator with parameter-dependent coefficients",
          "[parametric_qoperator]") {

  using namespace static_indices;

  auto expr1 = 3*c_dag<int_complex>("dn");
  auto expr2 = 3*c<int_complex>("up");
  auto expr = 2*expr1 - 2*expr2;

  auto hs = make_hilbert_space(expr);
  auto qop1 = make_param_qoperator(expr1, hs);
  auto qop2 = make_param_qoperator(expr2, hs);
  auto qop = make_param_qoperator(expr, hs);

  using state_vector = std::vector<int_complex>;

  SECTION("1 argument") {
    const state_vector in{1, 1, 1, 1};
    state_vector out(4, 0);

    qop1(in, out, 5);
    CHECK(out == state_vector{0, 15, 0, 15});
    qop2(in, out, 5);
    CHECK(out == state_vector{15, -15, 0, 0});
    qop(in, out, 5);
    CHECK(out == state_vector{-30, 60, 0, 30});
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
  }
}
