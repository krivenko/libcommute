/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2024 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

//
// Holstein model on a square lattice.
//

#include <libcommute/libcommute.hpp>

#include <cmath>
#include <iostream>
#include <string>

int main() {

  //
  // Let us define Hamiltonian of an electronic tight-binding model
  // on a square lattice.
  //

  // Number of lattice sites in each direction
  // (the total number of sites is N * N)
  int const N = 10;

  // Electron hopping constant - energy parameter of the TB model
  double const t = 2.0;

  // TB Hamiltonian as an expression with real coefficients and
  // statically-typed indices. In this case, the indices are a pair of
  // integers - lattice site coordinates - and a spin label ("up" or "down").
  libcommute::static_indices::expr_real<int, int, std::string> H_e;

  // Use functions c_dag() and c() that return fermionic creation/annihilation
  // operators.
  using libcommute::static_indices::c_dag;
  using libcommute::static_indices::c;

  // Are two sites neighbors along an axis with periodicity?
  auto neighbors = [](int i, int j) {
    return std::abs(i - j) == 1 || std::abs(i - j) == N;
  };

  // Iterate over spin projections
  for(auto spin : {"up", "down"}) {
    // Iterate over all lattice sites with coordinates i = (ix, iy)
    for(int ix = 0; ix < N; ++ix) {
      for(int iy = 0; iy < N; ++iy) {
        // Iterate over all lattice sites with coordinates j = (jx, jy)
        for(int jx = 0; jx < N; ++jx) {
          for(int jy = 0; jy < N; ++jy) {
            // Skip all pairs of lattice sites i and j that are not
            // nearest-neighbors. The modulus operation accounts for
            // periodic boundary conditions on the lattice.
            if((neighbors(ix, jx) && iy == jy) ||
               (ix == jx && neighbors(iy, jy))) {
              // Add a hopping term
              H_e += -t * c_dag(ix, iy, spin) * c(jx, jy, spin);
            }
          }
        }
      }
    }
  }

  // Frequency of the localized phonon
  double const w0 = 0.5;

  //
  // Hamiltonian of phonons localized at lattice sites.
  // We want this object to have the same type as H_e.
  //
  decltype(H_e) H_ph;

  // Use functions a_dag() and a() that return bosonic creation/annihilation
  // operators.
  using libcommute::static_indices::a_dag;
  using libcommute::static_indices::a;

  // Iterate over all lattice sites
  for(int ix = 0; ix < N; ++ix) {
    for(int iy = 0; iy < N; ++iy) {
      // Energy of the localized phonon at site (ix, iy)
      H_ph += w0 * a_dag(ix, iy, "") * a(ix, iy, "");
    }
  }

  // Electron-phonon coupling constant
  double const g = 0.1;

  //
  // Hamiltonian of electron-phonon coupling.
  //
  decltype(H_e) H_e_ph;

  // Use function n() that returns the fermionic number operator n = c_dag * c
  using libcommute::static_indices::n;

  // Iterate over spin projections
  for(auto spin : {"up", "down"}) {
    // Iterate over all lattice sites
    for(int ix = 0; ix < N; ++ix) {
      for(int iy = 0; iy < N; ++iy) {
        // Electron-phonon coupling at site (ix, iy)
        H_e_ph += g * n(ix, iy, spin) * (a_dag(ix, iy, "") + a(ix, iy, ""));
      }
    }
  }

  // Holstein Hamiltonian.
  auto H_H = H_e + H_ph + H_e_ph;

  // Print H_H. There will be quite a lot of terms for the 100-site lattice!
  std::cout << "H_H = " << H_H << '\n';

  // Check hermiticity of H_H
  std::cout << "H_H - H_H^\\dagger = " << (H_H - conj(H_H)) << '\n';

  // Check that H_H commutes with the total number of electrons
  decltype(H_H) N_e;
  for(auto spin : {"up", "down"}) {
    for(int ix = 0; ix < N; ++ix) {
      for(int iy = 0; iy < N; ++iy) {
        N_e += n(ix, iy, spin);
      }
    }
  }

  std::cout << "[H_H, N_e] = " << (H_H * N_e - N_e * H_H) << '\n';

  return 0;
}
