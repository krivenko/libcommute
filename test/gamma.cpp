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

#include <libcommute/expression/expression.hpp>

//
// Implement algebra of 4-dimensional \Gamma-matrices
//

namespace libcommute {

class generator_gamma : public generator<int> {

  using base = generator<int>;
  using linear_function_t = typename base::linear_function_t;

public:

  // 0 is the lowest algebra ID available to user-defined algebras
  virtual int algebra_id() const override { return 0; }

  // Value semantics
  generator_gamma(int index) : base(index) {}
  generator_gamma(generator_gamma const&) = default;
  generator_gamma(generator_gamma&&) noexcept = default;
  generator_gamma& operator=(generator_gamma const&) = default;
  generator_gamma& operator=(generator_gamma&&) noexcept = default;
  virtual ~generator_gamma() {}

  // Make a smart pointer that manages a copy of this generator
  virtual std::unique_ptr<base> clone() const override {
    return make_unique<generator_gamma>(*this);
  }

  // c = -1, f(g) = 2\eta(g1, g2)
  virtual double commute(base const& g2, linear_function_t & f) const override {
    assert(*this > g2);
    bool diag = base::equal(g2);
    f.set(diag * (std::get<0>(indices_) == 0 ? 2 : -2));
    return -1;
  }

  // (\Gamma^0)^2 = I_4
  // (\Gamma^k)^2 = -I_4 for k=1,2,3
  virtual bool collapse_power(int power, linear_function_t & f) const override {
    assert(power >= 2);
    auto coeff = [this](int p) {
      return std::get<0>(indices_) == 0 ? 1. : std::pow(-1, p/2);
    };
    if(power%2 == 0) {
      f.set(coeff(power));
    } else {
      f.set(0, clone(), coeff(power-1));
    }
    return true;
  }

  virtual bool has_vanishing_power(int power) const {
    assert(power >= 2);
    return false;
  }

  // Gamma^0 is Hermitian and Gamma^k are anti-Hermitian
  virtual void conj(linear_function_t & f) const override {
    f.set(0, clone(), std::get<0>(indices_) == 0 ? 1 : -1);
  }
};

} // namespace libcommute

using namespace libcommute;

TEST_CASE("Gamma matrices", "[gamma]") {
  using mon_type = monomial<int>;
  using expr_type = expression<std::complex<double>, int>;
  std::complex<double> I(0,1);

  // Metric
  auto eta = [](int mu, int nu) -> double {
    return (mu == nu) * (mu == 0 ? 1 : -1);
  };

  std::vector<expr_type> Gamma; // Gamma matrices
  std::vector<expr_type> Gammac; // Covariant Gamma matrices
  for(int mu : {0, 1, 2, 3}) {
    Gamma.emplace_back(expr_type(1.0, mon_type(generator_gamma(mu))));
    Gammac.emplace_back(expr_type(eta(mu, mu),
                                  mon_type(generator_gamma(mu))));
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
      CHECK(s == -2.0*Gamma[nu]);
    }
  }

  SECTION("Identity 3") {
    for(int nu = 0; nu < 4; ++nu) {
      for(int rho = 0; rho < 4; ++rho) {
        expr_type s;
        for(int mu = 0; mu < 4; ++mu)
          s += Gamma[mu] * Gamma[nu] * Gamma[rho] * Gammac[mu];
        CHECK(s == expr_type(4*eta(nu, rho)));
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
    return (i4 - i3)*(i4 - i2)*(i4 - i1)*(i3 - i2)*(i3 - i1)*(i2 - i1)/12;
  };

  SECTION("Identity 5") {
    for(int mu = 0; mu < 4; ++mu) {
      for(int nu = 0; nu < 4; ++nu) {
        for(int rho = 0; rho < 4; ++rho) {
          auto lhs = Gamma[mu] * Gamma[nu] * Gamma[rho];
          auto rhs = eta(mu, nu) * Gamma[rho] +
                     eta(nu, rho) * Gamma[mu] -
                     eta(mu, rho) * Gamma[nu];
          for(int sigma = 0; sigma < 4; ++sigma) {
            rhs += -I*eps(sigma, mu, nu, rho) * Gammac[sigma] * Gamma5;
          }
          CHECK(lhs == rhs);
        }
      }
    }
  }
}
