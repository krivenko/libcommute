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
#include <libcommute/loperator/elementary_space.hpp>
#include <libcommute/loperator/loperator.hpp>

#include <algorithm>
#include <array>
#include <complex>
#include <memory>
#include <vector>

std::complex<double> const I(0, 1);

namespace libcommute {

//
// 4-dimensional elementary space generated by gamma matrices
//

class elementary_space_gamma final : public elementary_space<int> {

  using base = elementary_space<int>;

public:
  // Value semantics
  elementary_space_gamma() : base(0) {}
  elementary_space_gamma(elementary_space_gamma const&) = default;
  elementary_space_gamma(elementary_space_gamma&&) noexcept = default;
  elementary_space_gamma& operator=(elementary_space_gamma const&) = default;
  elementary_space_gamma&
  operator=(elementary_space_gamma&&) noexcept = default;
  ~elementary_space_gamma() final = default;

  // cppcheck-suppress duplInheritedMember
  std::unique_ptr<base> clone() const final {
    return make_unique<elementary_space_gamma>(*this);
  }

  int algebra_id() const final { return libcommute::gamma; }
  int n_bits() const final { return 2; }
};

//
// Action of a gamma-matrix monomial
//

template <> class monomial_action<libcommute::gamma> {

  std::vector<int> index_sequence;

public:
  monomial_action(detail::monomial_range_t<int> const& m_range,
                  hilbert_space<int> const& hs) {
    for(auto it = m_range.first; it != m_range.second; ++it)
      index_sequence.push_back(std::get<0>(it->indices()));
    std::reverse(index_sequence.begin(), index_sequence.end());
  }

  inline bool act(sv_index_type& index, std::complex<double>& coeff) const {
    for(int i : index_sequence) {
      switch(i) {
      case 0: coeff *= std::array<double, 4>{1, 1, -1, -1}[index]; break;
      case 1:
        coeff *= std::array<double, 4>{-1, -1, 1, 1}[index];
        index = std::array<sv_index_type, 4>{3, 2, 1, 0}[index];
        break;
      case 2:
        coeff *= std::array<std::complex<double>, 4>{-I, I, I, -I}[index];
        index = std::array<sv_index_type, 4>{3, 2, 1, 0}[index];
        break;
      case 3:
        coeff *= std::array<std::complex<double>, 4>{-1, 1, 1, -1}[index];
        index = std::array<sv_index_type, 4>{2, 3, 0, 1}[index];
        break;
      }
    }
    return true;
  }
};

} // namespace libcommute

using namespace libcommute;

template <typename ExprType>
void check_loperator(ExprType const& expr1, ExprType const& expr2) {
  hilbert_space<int> hs{elementary_space_gamma()};

  using libcommute::gamma;
  auto lop1 = loperator<std::complex<double>, gamma>(expr1, hs);
  auto lop2 = loperator<std::complex<double>, gamma>(expr2, hs);

  std::vector<std::complex<double>> in(4), out1(4), out2(4);
  std::fill(in.begin(), in.end(), .0);

  for(sv_index_type in_index : {0, 1, 2, 3}) {
    in[in_index] = 1;
    CHECK(lop1(in) == lop2(in));
    CHECK(lop1 * in == lop2 * in);
    lop1(in, out1);
    lop2(in, out2);
    CHECK(out1 == out2);
    in[in_index] = 0; // cppcheck-suppress unreadVariable
  }
}

TEST_CASE("loperator for expressions with gamma-matrices",
          "[loperator_gamma]") {
  using mon_type = monomial<int>;
  using expr_type = expression<std::complex<double>, int>;

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
        auto lhs = Gamma[mu] * Gamma[nu] + Gamma[nu] * Gamma[mu];
        auto rhs = expr_type(2 * eta(mu, nu));
        check_loperator(lhs, rhs);
      }
    }
  }

  SECTION("Hermitian conjugate") {
    check_loperator(conj(Gamma[0]), Gamma[0]);
    for(int k = 1; k < 4; ++k) {
      check_loperator(conj(Gamma[k]), -Gamma[k]);
    }
    for(int mu = 0; mu < 4; ++mu) {
      check_loperator(conj(Gamma[mu]), Gamma[0] * Gamma[mu] * Gamma[0]);
    }
  }

  auto Gamma5 = I * Gamma[0] * Gamma[1] * Gamma[2] * Gamma[3];

  SECTION("Gamma^5") {
    for(int mu = 0; mu < 4; ++mu) {
      check_loperator(Gamma5 * Gamma[mu] + Gamma[mu] * Gamma5, expr_type());
    }
    check_loperator(Gamma5 * Gamma5, expr_type(1));
    check_loperator(conj(Gamma5), Gamma5);
  }

  //
  // https://en.wikipedia.org/wiki/Gamma_matrices#Miscellaneous_identities
  //

  SECTION("Identity 1") {
    expr_type s;
    for(int mu = 0; mu < 4; ++mu)
      s += Gamma[mu] * Gammac[mu];
    check_loperator(s, expr_type(4));
  }

  SECTION("Identity 2") {
    for(int nu = 0; nu < 4; ++nu) {
      expr_type s;
      for(int mu = 0; mu < 4; ++mu)
        s += Gamma[mu] * Gamma[nu] * Gammac[mu];
      check_loperator(s, -2.0 * Gamma[nu]);
    }
  }

  SECTION("Identity 3") {
    for(int nu = 0; nu < 4; ++nu) {
      for(int rho = 0; rho < 4; ++rho) {
        expr_type s;
        for(int mu = 0; mu < 4; ++mu)
          s += Gamma[mu] * Gamma[nu] * Gamma[rho] * Gammac[mu];
        check_loperator(s, expr_type(4 * eta(nu, rho)));
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
          check_loperator(s, -2.0 * Gamma[sigma] * Gamma[rho] * Gamma[nu]);
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
          check_loperator(lhs, rhs);
        }
      }
    }
  }
}
