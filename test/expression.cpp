/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <catch.hpp>

#include "my_complex.hpp"
#include "print_matcher.hpp"

#include <libcommute/expression/expression.hpp>
#include <libcommute/expression/factories.hpp>
#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/generator_spin.hpp>

#include <array>

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

    expr_monomial.clear();
    CHECK(expr_monomial.size() == 0);

    SECTION("S_z products") {
      monomial<int> mon_sz(make_fermion(true, 1),
                           make_spin(spin_component::z, 1),
                           make_spin(spin_component::z, 1),
                           make_spin(spin_component::z, 1),
                           make_spin(spin_component::z, 1),
                           make_fermion(true, 2),
                           make_spin(spin_component::z, 2),
                           make_spin(spin_component::z, 2),
                           make_spin(spin_component::z, 2),
                           make_spin(spin_component::z, 2),
                           make_spin(spin_component::z, 3),
                           make_spin(spin_component::z, 3),
                           make_spin(spin_component::z, 3),
                           make_fermion(true, 4),
                           make_spin(spin_component::z, 4),
                           make_spin(spin_component::z, 4),
                           make_spin(spin_component::z, 4));
      expr_real<int> expr_sz(1.0, mon_sz);
      CHECK(expr_sz == c_dag(1) * 0.25 * 0.25 * c_dag(2) * 0.25 * 0.25 * 0.25 *
                           S_z(3) * c_dag(4) * 0.25 * S_z(4));
    }
  }

  SECTION("Assignment") {
    expr_real<int, std::string> expr_const(2);
    expr_real<int, std::string> expr0, expr1;
    expr0 = expr_const;
    CHECK_THAT(expr0, Prints<decltype(expr0)>("2"));
    expr1 = std::move(expr_const);
    CHECK_THAT(expr1, Prints<decltype(expr1)>("2"));
  }

  SECTION("Unary minus") {
    auto expr_r = c_dag(1, "up");
    CHECK_THAT(-expr_r, Prints<decltype(expr_r)>("-1*C+(1,up)"));

    auto expr_static_int = c_dag<my_complex>(1, "up");
    CHECK_THAT(-expr_static_int,
               Prints<decltype(expr_static_int)>("{-1,0}*C+(1,up)"));
  }

  SECTION("const_iterator") {
    using expr_type = expression<double, int, std::string>;
    using mon_type = expr_type::monomial_t;

    expr_type expr0;
    CHECK(expr0.begin() == expr0.end());
    CHECK_FALSE(expr0.begin() != expr0.end());

    auto expr = 4.0 * c_dag(1, "up") * c(2, "dn") + 1.0 + 3.0 * a(0, "x") +
                2.0 * a_dag(0, "y");

    CHECK_FALSE(expr.begin() == expr.end());
    CHECK(expr.begin() != expr.end());

    std::vector<mon_type> ref_mons = {
        mon_type(),
        mon_type(make_boson(true, 0, "y")),
        mon_type(make_boson(false, 0, "x")),
        mon_type(make_fermion(true, 1, "up"), make_fermion(false, 2, "dn"))};
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
      // cppcheck-suppress postfixOperator
      for(it = expr.begin(); it != expr.end(); it++, n++) {
        auto val = *it;
        CHECK(val.monomial == ref_mons[n]);
        CHECK(val.coeff == ref_coeffs[n]);
        CHECK(it->monomial == ref_mons[n]);
        CHECK(it->coeff == ref_coeffs[n]);
      }
      // Backward iteration
      // cppcheck-suppress postfixOperator
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
    auto expr = 4.0 * c_dag(1, "up") * c(2, "dn") + 1.0 + 3.0 * a(0, "x") +
                2.0 * a_dag(0, "y");
    using mon_type = decltype(expr)::monomial_t;

    // Multiply coefficients in front of bosonic operators by 2*I
    auto f = [](mon_type const& m, double c) -> std::complex<double> {
      return (m.size() > 0 && is_boson(m[0])) ? std::complex<double>(0, 2 * c) :
                                                .0;
    };

    auto new_expr = transform(expr, f);

    std::vector<mon_type> ref_mons = {
        mon_type(make_boson(true, 0, "y")),
        mon_type(make_boson(false, 0, "x")),
    };
    std::vector<std::complex<double>> ref_coeffs = {
        std::complex<double>(0, 4.0),
        std::complex<double>(0, 6.0)};

    CHECK(new_expr.size() == 2);

    int n = 0;
    for(auto it = new_expr.begin(); it != new_expr.end(); ++it, ++n) {
      auto val = *it;
      CHECK(val.monomial == ref_mons[n]);
      CHECK(val.coeff == ref_coeffs[n]);
    }
  }

  SECTION("conj()") {
    auto expr = 4.0 * c_dag(1, "up") * c(2, "dn") + 1.0 + 3.0 * a(0, "x") +
                std::complex<double>(0, 2) * a_dag(0, "y");
    auto ref = 4.0 * c_dag(2, "dn") * c(1, "up") + 1.0 + 3.0 * a_dag(0, "x") +
               std::complex<double>(0, -2) * a(0, "y");

    CHECK(conj(expr) == ref);
  }

  SECTION("Products of Spin-1/2 operators") {
    CHECK(S_z() * S_z() == expr_real<>(0.25));
    CHECK(S_p() * S_p() == expr_real<>(0));
    CHECK(S_m() * S_m() == expr_real<>(0));

    CHECK(S_p() * S_z() == -0.5 * S_p());
    CHECK(S_z() * S_m() == -0.5 * S_m());
    CHECK(S_z() * S_p() == 0.5 * S_p());
    CHECK(S_m() * S_z() == 0.5 * S_m());
    CHECK(S_p() * S_m() == 0.5 + S_z());
    CHECK(S_m() * S_p() == 0.5 - S_z());
  }

  SECTION("Powers of S_z") {
    expr_real<> sum;
    for(int n = 1; n <= 11; ++n) {
      auto p = c_dag();
      for(int i = 0; i < n; ++i)
        p *= S_z();
      p *= a();
      sum += p;
    }
    CHECK(sum == (c_dag() * (341.0 / 1024 + (1365.0 / 1024) * S_z()) * a()));
  }

  SECTION("Heisenberg chain") {
    using expr_t = expression<std::complex<double>, int>;
    using vec_expr_t = std::array<expr_t, 3>;

    // Addition of vectors
    auto add = [](vec_expr_t const& S1, vec_expr_t const& S2) -> vec_expr_t {
      return {S1[0] + S2[0], S1[1] + S2[1], S1[2] + S2[2]};
    };
    // Dot-product of vectors
    auto dot = [](vec_expr_t const& S1, vec_expr_t const& S2) -> expr_t {
      return S1[0] * S2[0] + S1[1] * S2[1] + S1[2] * S2[2];
    };
    // Cross-product of vectors
    auto cross = [](vec_expr_t const& S1, vec_expr_t const& S2) -> vec_expr_t {
      return {S1[1] * S2[2] - S1[2] * S2[1],
              S1[2] * S2[0] - S1[0] * S2[2],
              S1[0] * S2[1] - S1[1] * S2[0]};
    };

    int const N = 6;
    std::vector<vec_expr_t> S;
    S.reserve(N);
    for(int i = 0; i < N; ++i)
      S.push_back({S_x(i), S_y(i), S_z(i)});

    expr_t H;
    for(int i = 0; i < N; ++i) {
      H += dot(S[i], S[(i + 1) % N]);
    }

    vec_expr_t S_tot;
    for(int i = 0; i < N; ++i)
      S_tot = add(S_tot, S[i]);

    // H must commute with the total spin
    CHECK((H * S_tot[0] - S_tot[0] * H).size() == 0);
    CHECK((H * S_tot[1] - S_tot[1] * H).size() == 0);
    CHECK((H * S_tot[2] - S_tot[2] * H).size() == 0);

    expr_t Q3;
    for(int i = 0; i < N; ++i) {
      Q3 += dot(cross(S[i], S[(i + 1) % N]), S[(i + 2) % N]);
    }

    // Q3 is a higher-order integral of motion
    CHECK((H * Q3 - Q3 * H).size() == 0);
  }
}
