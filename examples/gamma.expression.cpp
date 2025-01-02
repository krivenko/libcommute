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

//
// Custom algebra of Dirac gamma matrices: polynomial expressions
//

#include "gamma.hpp"

//
// Check that expressions with gamma-matrices behave as expected
//

int main() {

  using namespace libcommute;

  // Minkowski metric tensor
  auto eta = [](int mu, int nu) -> double {
    if(mu != nu) return 0;
    return mu == 0 ? 1.0 : -1.0;
  };

  // Check the commutation relations
  for(int mu = 0; mu < 4; ++mu) {
    for(int nu = 0; nu < 4; ++nu) {
      auto gamma_mu = make_gamma(mu);
      auto gamma_nu = make_gamma(nu);

      std::cout << "{" << gamma_mu << ", " << gamma_nu << "}"
                << " - 2\\eta(" << mu << ", " << nu << "} = "
                << (gamma_mu * gamma_nu + gamma_nu * gamma_mu - 2 * eta(mu, nu))
                << '\n';
    }
  }

  // \gamma^5
  std::complex<double> const I(0, 1);
  auto gamma5 =
      I * make_gamma(0) * make_gamma(1) * make_gamma(2) * make_gamma(3);

  // \gamma^5 is Hermitian ...
  std::cout << "gamma5 - conj(gamma5) = " << (gamma5 - conj(gamma5)) << '\n';

  // ... and anti-commutes with \gamma^\mu.
  for(int mu = 0; mu < 4; ++mu) {
    auto gamma_mu = make_gamma(mu);
    std::cout << "{gamma5, " << gamma_mu
              << "} = " << (gamma5 * gamma_mu + gamma_mu * gamma5) << '\n';
  }
}
