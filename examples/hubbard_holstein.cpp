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

#include <cstdlib>
#include <iostream>

#include <libcommute/libcommute.hpp>

// This small program illustrates how to use libcommute to construct Hamiltonian
// of a lattice model.

// The system considered here is the Hubbard-Holstein model on a square 10x10
// lattice with nearest-neighbor hopping.
//
// \hat H = -\sum_{\langle i,j \rangle, \sigma}
//            t_{ij} c^\dagger_{i,\sigma} c_{j,\sigma}
//          + U \sum_i n_{i,\uparrow} n_{i,\downarrow}
//          + \sum_i a^\dagger_i a_i
//          + g \sum_{i,\sigma} n_{i,\sigma} (a^\dagger_i + a_i)
//
using namespace libcommute;

int main() {

  //
  // Parameters of the system
  //

  // Linear sizes of the lattice
  const int Nx = 10;
  const int Ny = 10;

  // Hopping constant
  const double t = 0.5;
  // Coulomb repulsion
  const double U = 2.0;
  // Electron-phonon coupling constant
  const double g = 0.1;

  // Expression with real coefficients 'H' will represent the Hamiltonian.
  // It is initially set to zero by its default-constructor.
  // Every creation and annihilation operator met in the expression must
  // carry two integer (coordinates of a lattice site) and one string
  // index.
  expression<double,               // type of expression's coefficients
             int, int, std::string // types of operator indices
             > H;

  // The following 'factory' functions make quantum operators with
  // statically typed indices and real coefficients.
  using static_indices::real::c_dag; // Create an electron
  using static_indices::real::c;     // Destroy an electron
  using static_indices::real::n;     // Number of electrons
  using static_indices::real::a_dag; // Create a phonon
  using static_indices::real::a;     // Destroy a phonon

  // Hopping terms of H
  for(auto spin : {"up", "down"}) {
    for(int ix = 0; ix < Nx; ++ix) {
      for(int iy = 0; iy < Ny; ++iy) {
        for(int jx = 0; jx < Nx; ++jx) {
          for(int jy = 0; jy < Ny; ++jy) {
            // Skip all pairs of lattice sites (ix,iy) and (jx,jy) that are
            // not nearest-neighbors.
            if((std::abs(ix - jx) % Nx == 1 && iy == jy) ||
               (ix == jx && std::abs(iy - jy) % Ny == 1)
            ) {
              // Add a hopping term
              H += -t * c_dag(ix, iy, spin) * c(jx, jy, spin);
            }
          }
        }
      }
    }
  }

  // Coulomb repulsion terms
  for(int ix = 0; ix < Nx; ++ix)
    for(int iy = 0; iy < Ny; ++iy) {
      H += U * n(ix, iy, "up") * n(ix, iy, "down");
    }

  // Energy of phonons
  for(int ix = 0; ix < Nx; ++ix)
    for(int iy = 0; iy < Ny; ++iy) {
      // The spin index is left blank for bosonic operators
      H += a_dag(ix, iy, "") * a(ix, iy, "");
    }

  // Electron-phonon coupling
  for(auto spin : {"up", "down"}) {
    for(int ix = 0; ix < Nx; ++ix) {
      for(int iy = 0; iy < Ny; ++iy) {
        H += g * n(ix, iy, spin) * (a_dag(ix, iy, "") + a(ix, iy, ""));
      }
    }
  }

  // Total number of terms (monomials) in 'H'.
  std::cout << "Total number of terms in H: " << H.size() << std::endl;
  // Is H Hermitian?
  std::cout << "H^\\dagger - H = " << (conj(H) - H) << std::endl;

  // Does H commute with N and S_z?
  decltype(H) N, S_z;
  for(int ix = 0; ix < Nx; ++ix) {
    for(int iy = 0; iy < Ny; ++iy) {
      N += n(ix, iy, "up") + n(ix, iy, "down");
      S_z += 0.5 * (n(ix, iy, "up") - n(ix, iy, "down"));
    }
  }
  std::cout << "[H, N] = " << (H * N - N * H) << std::endl;
  std::cout << "[H, S_z] = " << (H * S_z - S_z * H) << std::endl;

  // Iterate over all terms in 'H' and print those of degree 3.
  //
  // Monomials of degree 3 come from the electron-phonon coupling and
  // are products of two fermionic and one bosonic operators.
  for(auto const& term: H) {
    if(term.monomial.size() == 3) {
      // term.coeff is coefficient in front of the monomial
      std::cout << term.monomial << " => " << term.coeff << std::endl;
    }
  }

  return 0;
}
