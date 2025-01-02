/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2025 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <catch.hpp>

#include "./gamma.hpp"

#include <libcommute/algebra_ids.hpp>
#include <libcommute/expression/expression.hpp>

using namespace libcommute;

TEST_CASE("Gamma matrices", "[gamma]") {
  using mon_type = monomial<int>;
  using expr_type = expression<std::complex<double>, int>;
  std::complex<double> I(0, 1);

  // Metric
  auto eta = [](int mu, int nu) -> double {
    return (mu == nu) * (mu == 0 ? 1 : -1);
  };

  std::vector<expr_type> Gamma;  // Gamma matrices
  std::vector<expr_type> Gammac; // Covariant Gamma matrices
  for(int mu : {0, 1, 2, 3}) {
    Gamma.emplace_back(1.0, mon_type(generator_gamma(mu)));
    Gammac.emplace_back(eta(mu, mu), mon_type(generator_gamma(mu)));
  }

  SECTION("Commutation relations") {
    for(int mu = 0; mu < 4; ++mu) {
      for(int nu = 0; nu < 4; ++nu) {
        CHECK(Gamma[mu] * Gamma[nu] + Gamma[nu] * Gamma[mu] ==
              expr_type(2 * eta(mu, nu)));
      }
    }
  }

  SECTION("Hermitian conjugate") {
    CHECK(conj(Gamma[0]) == Gamma[0]);
    for(int k = 1; k < 4; ++k) {
      CHECK(conj(Gamma[k]) == -Gamma[k]);
    }
    for(int mu = 0; mu < 4; ++mu) {
      CHECK(conj(Gamma[mu]) == Gamma[0] * Gamma[mu] * Gamma[0]);
    }
  }

  auto Gamma5 = I * Gamma[0] * Gamma[1] * Gamma[2] * Gamma[3];

  SECTION("Gamma^5") {
    for(int mu = 0; mu < 4; ++mu) {
      CHECK(Gamma5 * Gamma[mu] + Gamma[mu] * Gamma5 == expr_type());
    }
    CHECK(Gamma5 * Gamma5 == expr_type(1));
    CHECK(conj(Gamma5) == Gamma5);
  }

  //
  // https://en.wikipedia.org/wiki/Gamma_matrices#Miscellaneous_identities
  //

  SECTION("Identity 1") {
    expr_type s;
    for(int mu = 0; mu < 4; ++mu)
      s += Gamma[mu] * Gammac[mu];
    CHECK(s == expr_type(4));
  }

  SECTION("Identity 2") {
    for(int nu = 0; nu < 4; ++nu) {
      expr_type s;
      for(int mu = 0; mu < 4; ++mu)
        s += Gamma[mu] * Gamma[nu] * Gammac[mu];
      CHECK(s == -2.0 * Gamma[nu]);
    }
  }

  SECTION("Identity 3") {
    for(int nu = 0; nu < 4; ++nu) {
      for(int rho = 0; rho < 4; ++rho) {
        expr_type s;
        for(int mu = 0; mu < 4; ++mu)
          s += Gamma[mu] * Gamma[nu] * Gamma[rho] * Gammac[mu];
        CHECK(s == expr_type(4 * eta(nu, rho)));
      }
    }
  }

  SECTION("Identity 4") {
    for(int nu = 0; nu < 4; ++nu) {
      for(int rho = 0; rho < 4; ++rho) {
        for(int sigma = 0; sigma < 4; ++sigma) {
          expr_type s;
          for(int mu = 0; mu < 4; ++mu)
            s += Gamma[mu] * Gamma[nu] * Gamma[rho] * Gamma[sigma] * Gammac[mu];
          CHECK(s == -2.0 * Gamma[sigma] * Gamma[rho] * Gamma[nu]);
        }
      }
    }
  }

  // 4D Levi-Civita symbol
  auto eps = [](int i1, int i2, int i3, int i4) -> double {
    return (i4 - i3) * (i4 - i2) * (i4 - i1) * (i3 - i2) * (i3 - i1) *
           (i2 - i1) / 12;
  };

  SECTION("Identity 5") {
    for(int mu = 0; mu < 4; ++mu) {
      for(int nu = 0; nu < 4; ++nu) {
        for(int rho = 0; rho < 4; ++rho) {
          auto lhs = Gamma[mu] * Gamma[nu] * Gamma[rho];
          auto rhs = eta(mu, nu) * Gamma[rho] + eta(nu, rho) * Gamma[mu] -
                     eta(mu, rho) * Gamma[nu];
          for(int sigma = 0; sigma < 4; ++sigma) {
            rhs += -I * eps(sigma, mu, nu, rho) * Gammac[sigma] * Gamma5;
          }
          CHECK(lhs == rhs);
        }
      }
    }
  }
}
