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
// Periodic spin-1/2 Heisenberg chain and its integrals of motion
//
// Expressions for the integrals of motion are taken from
//
//   "Quantum Integrals of Motion for the Heisenberg Spin Chain",
//   M. P. Grabowski and P. Mathieu
//   Mod. Phys. Lett. A, Vol. 09, No. 24, pp. 2197-2206 (1994),
//   https://doi.org/10.1142/S0217732394002057
//

#include <array>
#include <iostream>
#include <vector>

#include <libcommute/libcommute.hpp>

using namespace libcommute;

// Type of the spin-1/2 Heisenberg Hamiltonian and its charges.
// All spin operators are going to carry one integer index - the site index.
using expr_t = static_indices::expr_complex<int>;

// Spins S_i are operator-valued 3-dimensional vectors. Here, we use std::array
// as a simplistic implementation of such vectors.
using vector_expr_t = std::array<expr_t, 3>;

// Addition of two vectors.
vector_expr_t operator+(vector_expr_t const& S1, vector_expr_t const& S2) {
  return {S1[0] + S2[0], S1[1] + S2[1], S1[2] + S2[2]};
}

// Scalar product of two vectors.
expr_t dot(vector_expr_t const& S1, vector_expr_t const& S2) {
  return S1[0] * S2[0] + S1[1] * S2[1] + S1[2] * S2[2];
}

// Cross product of two vectors
vector_expr_t cross(vector_expr_t const& S1, vector_expr_t const& S2) {
  return {S1[1] * S2[2] - S1[2] * S2[1],
          S1[2] * S2[0] - S1[0] * S2[2],
          S1[0] * S2[1] - S1[1] * S2[0]};
}

int main() {

  // For functions S_x(), S_y() and S_z().
  // Note that the x/y functions exist only for the complex expressions.
  using namespace static_indices;

  // Number of spins in the chain
  int const N = 20;
  // Heisenberg exchange constant
  double const g = 2;

  // List of spin operators {S_0, S_1, ..., S_{N-1}}
  std::vector<vector_expr_t> S;
  S.reserve(N);
  for(int i = 0; i < N; ++i)
    S.push_back({S_x(i), S_y(i), S_z(i)});

  // Hamiltonian of the spin-1/2 Heisenberg chain.
  expr_t H;
  for(int i = 0; i < N; ++i) {
    // Index shift modulo N ensures periodic boundary conditions.
    H += g * dot(S[i], S[(i + 1) % N]);
  }

  // Total spin of the chain
  vector_expr_t S_tot;
  for(int i = 0; i < N; ++i)
    S_tot = S_tot + S[i];

  // All three components of S commute with the Hamiltonian.
  std::cout << "[H, S_x] = " << (H * S_tot[0] - S_tot[0] * H) << '\n';
  std::cout << "[H, S_y] = " << (H * S_tot[1] - S_tot[1] * H) << '\n';
  std::cout << "[H, S_z] = " << (H * S_tot[2] - S_tot[2] * H) << '\n';

  // Higher charge Q_3 (1st line of Eq. (10)).
  expr_t Q3;
  for(int i = 0; i < N; ++i) {
    Q3 += dot(cross(S[i], S[(i + 1) % N]), S[(i + 2) % N]);
  }
  std::cout << "[H, Q3] = " << (H * Q3 - Q3 * H) << '\n';

  // Higher charge Q_4 (2nd line of Eq. (10)).
  expr_t Q4;
  for(int i = 0; i < N; ++i) {
    Q4 += 4.0 * dot(cross(cross(S[i], S[(i + 1) % N]), S[(i + 2) % N]),
                    S[(i + 3) % N]);
    Q4 += dot(S[i], S[(i + 2) % N]);
  }
  std::cout << "[H, Q4] = " << (H * Q4 - Q4 * H) << '\n';

  // Higher charge Q_5 (3rd line of Eq. (10)).
  expr_t Q5;
  for(int i = 0; i < N; ++i) {
    Q5 += 4.0 * dot(cross(cross(cross(S[i], S[(i + 1) % N]), S[(i + 2) % N]),
                          S[(i + 3) % N]),
                    S[(i + 4) % N]);
    Q5 += dot(cross(S[i], S[(i + 2) % N]), S[(i + 3) % N]);
    Q5 += dot(cross(S[i], S[(i + 1) % N]), S[(i + 3) % N]);
  }
  std::cout << "[H, Q5] = " << (H * Q5 - Q5 * H) << '\n';

  // Check that the higher charges pairwise commute.
  std::cout << "[Q3, Q4] = " << (Q3 * Q4 - Q4 * Q3) << '\n';
  std::cout << "[Q3, Q5] = " << (Q3 * Q5 - Q5 * Q3) << '\n';
  std::cout << "[Q4, Q5] = " << (Q4 * Q5 - Q5 * Q4) << '\n';

  return 0;
}
