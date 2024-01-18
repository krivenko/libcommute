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

//
// Spectrum of the Jaynes–Cummings model
//
//   "Comparison of quantum and semiclassical radiation theories
//   with application to the beam maser",
//   E. T. Jaynes and F.W. Cummings,
//   Proc. IEEE, Vol. 51, issue 1, pp. 89-109 (1963),
//   https://doi.org/10.1109/PROC.1963.1664
//

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <libcommute/libcommute.hpp>

using namespace libcommute;

int main() {

  // Planck constant
  double const hbar = 1.0;
  // Angular frequency of the cavity mode
  double const w_c = 1.0;
  // Atomic transition frequency
  double const w_a = 1.1;
  // Atom-field coupling constant
  double const Omega = 0.2;
  // Detuning
  double const delta = w_a - w_c;

  // For a_dag(), a(), S_p(), S_m() and S_z()
  using namespace static_indices;

  // The Jaynes-Cummings Hamiltonian
  auto H = hbar * w_c * a_dag() * a() + hbar * w_a * S_z() +
           hbar * 0.5 * Omega * (a() * S_p() + a_dag() * S_m());

  // Check that the number operator commutes with the Hamiltonian
  auto N = a_dag() * a() + S_z();
  std::cout << "[H, N] = " << (H * N - N * H) << std::endl;

  // Complete Hilbert space of the system is a product of the isospin-1/2 space
  // and the bosonic space truncated to 2^10 - 1 excitations.
  auto es_spin = static_indices::make_space_spin(0.5);
  auto es_boson = static_indices::make_space_boson(10);

  hilbert_space</* Our operators have no indices */> hs(es_spin, es_boson);

  // Ground state energy
  double const E_gs = -hbar * w_a / 2;

  // Excited energy levels E_{\pm}(n)
  auto energy = [=](int pm, int n) {
    return hbar * w_c * (n - 0.5) +
           pm * 0.5 * hbar * std::sqrt(delta * delta + n * Omega * Omega);
  };

  // Mixing angle \alpha_n
  auto alpha = [&](int n) { return std::atan(Omega * std::sqrt(n) / delta); };

  // Two state vectors in the complete Hilbert space
  std::vector<double> phi(hs.dim()), psi(hs.dim());

  // L^2 norm of |\spi> - E|\phi>
  // Will be used to verify correctness of the eigenpairs
  auto diff = [&](std::vector<double> const& psi,
                  std::vector<double> const& phi,
                  double E) {
    double d = 0;
    for(unsigned int i = 0; i < hs.dim(); ++i) {
      d += std::pow(psi[i] - E * phi[i], 2);
    }
    return std::sqrt(d);
  };

  // loperator object corresponding to H
  auto Hop = make_loperator(H, hs);

  //
  // Verify that |0, g> is the ground state
  //

  // Index of the |0, g> basis vector
  auto gs_index = hs.basis_state_index(es_boson, 0) + // n = 0
                  hs.basis_state_index(es_spin, 0);   // 0 -> g

  // Set |\phi> = |0, g>
  std::fill(phi.begin(), phi.end(), 0);
  phi[gs_index] = 1;

  psi = Hop(phi);

  // Check |\psi> \approx E_{gs} |\phi>
  std::cout << "|||0,g> - E_g|0,g>||_{L^2} = " << diff(psi, phi, E_gs)
            << std::endl;

  //
  // Now, do a similar check for the Jaynes-Cummings ladder doublets.
  //

  // Check the first 20 doublets
  for(int n = 1; n <= 20; ++n) {

    // Index of basis state |n - 1, e>
    auto n1_e_index = hs.basis_state_index(es_boson, n - 1) + // n - 1
                      hs.basis_state_index(es_spin, 1);       // 1 -> e
    // Index of basis state |n, g>
    auto n_g_index = hs.basis_state_index(es_boson, n) + // n
                     hs.basis_state_index(es_spin, 0);   // 0 -> g

    // Set \phi to be the "+" dressed state
    std::fill(phi.begin(), phi.end(), 0);
    phi[n1_e_index] = std::cos(alpha(n) / 2);
    phi[n_g_index] = std::sin(alpha(n) / 2);

    psi = Hop(phi);

    std::cout << "|||" << n << ",+> - E_+(" << n << ")|" << n
              << ",+>||_{L^2} = " << diff(psi, phi, energy(1, n)) << std::endl;

    // Set \phi to be the "-" dressed state
    std::fill(phi.begin(), phi.end(), 0);
    phi[n1_e_index] = std::sin(alpha(n) / 2);
    phi[n_g_index] = -std::cos(alpha(n) / 2);

    psi = Hop(phi);

    std::cout << "|||" << n << ",-> - E_-(" << n << ")|" << n
              << ",->||_{L^2} = " << diff(psi, phi, energy(-1, n)) << std::endl;
  }

  return 0;
}
