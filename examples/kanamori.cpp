/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSl and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPl was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

//
// Kanamori interaction Hamiltonian and its expression in terms of
// integrals of motion N, S^2 and L^2.
//
//   "Strong Correlations from Hund's Coupling",
//   A. Georges, L. de' Medici and J. Mravlje,
//   Annu. Rev. Condens. Matter Phys. 2013. 4:137–78,
//   https://doi.org/10.1146/annurev-conmatphys-020911-125045
//

#include <cmath>
#include <iostream>
#include <string>

#include <libcommute/libcommute.hpp>

using namespace libcommute;

int main() {

  // For functions c_dag(), c() and n().
  using namespace static_indices;

  // Orbital degeneracy of the shell (t_{2g} triplet)
  int const n_orbs = 3;
  // Coulomb integrals U and J
  double const U = 4.0;
  double const J = 0.2;

  // Kanamori Hamiltonian (Eq. (2))
  static_indices::expr_real<int /* orbital index */,
                            std::string /* spin index */>
      H_K;

  // Intraorbital density-density interaction terms
  for(int m = 0; m < n_orbs; ++m) {
    H_K += U * n(m, "up") * n(m, "down");
  }

  // Interorbital density-density interaction terms (different spins)
  for(int m1 = 0; m1 < n_orbs; ++m1) {
    for(int m2 = 0; m2 < n_orbs; ++m2) {
      if(m1 == m2) continue;
      H_K += (U - 2 * J) * n(m1, "up") * n(m2, "down");
    }
  }

  // Interorbital density-density interaction terms (same spin)
  for(int m1 = 0; m1 < n_orbs; ++m1) {
    for(int m2 = 0; m2 < n_orbs; ++m2) {
      if(m1 >= m2) continue;
      H_K += (U - 3 * J) * n(m1, "up") * n(m2, "up");
      H_K += (U - 3 * J) * n(m1, "down") * n(m2, "down");
    }
  }

  // Spin-flip terms
  for(int m1 = 0; m1 < n_orbs; ++m1) {
    for(int m2 = 0; m2 < n_orbs; ++m2) {
      if(m1 == m2) continue;
      H_K += -J * c_dag(m1, "up") * c(m1, "down") * c_dag(m2, "down") *
             c(m2, "up");
    }
  }

  // Pair-hopping terms
  for(int m1 = 0; m1 < n_orbs; ++m1) {
    for(int m2 = 0; m2 < n_orbs; ++m2) {
      if(m1 == m2) continue;
      H_K +=
          J * c_dag(m1, "up") * c_dag(m1, "down") * c(m2, "down") * c(m2, "up");
    }
  }

  // Print the Hamiltonian
  std::cout << "H_K = " << H_K << std::endl;

  //
  // Integrals of motion N, S^2 and L^2 (Eq. (4)).
  //

  std::complex<double> const I(0, 1);

  // Total number of particles
  decltype(H_K) N;
  for(int m = 0; m < n_orbs; ++m) {
    N += n(m, "up") + n(m, "down");
  }

  // Total spin operators S_x, S_y, S_z.
  static_indices::expr_complex< // to allow for complex coefficients in S_x, S_y
      int,
      std::string>
      Sx, Sy, Sz;
  for(int m = 0; m < n_orbs; ++m) {
    Sx += 0.5 * (c_dag(m, "up") * c(m, "down") + c_dag(m, "down") * c(m, "up"));
    Sy += 0.5 * I *
          (c_dag(m, "down") * c(m, "up") - c_dag(m, "up") * c(m, "down"));
    Sz += 0.5 * (n(m, "up") - n(m, "down"));
  }
  // Operator S^2 = S_x S_x + S_y S_y + S_z S_z
  auto S2 = Sx * Sx + Sy * Sy + Sz * Sz;

  // Levi-Civita symbol \epsilon_{ijk}
  auto eps = [](int i, int j, int k) -> double {
    return (j - i) * (k - j) * (k - i) / 2;
  };

  // Orbital isospin generators L_x, L_y, L_z.
  static_indices::expr_complex< // to allow for complex coefficients
      int,
      std::string>
      Lx, Ly, Lz;

  for(std::string spin : {"up", "down"}) {
    for(int m1 = 0; m1 < n_orbs; ++m1) {
      for(int m2 = 0; m2 < n_orbs; ++m2) {
        Lx += I * eps(0, m1, m2) * c_dag(m1, spin) * c(m2, spin);
        Ly += I * eps(1, m1, m2) * c_dag(m1, spin) * c(m2, spin);
        Lz += I * eps(2, m1, m2) * c_dag(m1, spin) * c(m2, spin);
      }
    }
  }
  // Operator L^2 = L_x L_x + L_y L_y + L_z L_z
  auto L2 = Lx * Lx + Ly * Ly + Lz * Lz;

  // Hamiltonian as a function of N, S^2 and L^2 (Eq. (7)).
  auto H_t2g = (U - 3 * J) / 2 * N * (N - 1.0) - 2 * J * S2 - J / 2 * L2 +
               (5.0 / 2) * J * N;

  // Must be zero
  std::cout << "H_K - H_t2g = " << (H_K - H_t2g) << std::endl;

  return 0;
}
