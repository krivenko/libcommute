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

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <libcommute/libcommute.hpp>

// In this program we show how to construct a matrix representation of
// an expression. In this particular case, of the Heisenberg exchange
// interaction between two spins 1/2,
//
// \hat H = \mathbf{S}^1 \cdot \mathbf{S}^2 =
//    \frac{1}{2}\left[ S^1_+ S^2_- + S^1_- S^2_+ \right] + S^1_z S^2_z
//
using namespace libcommute;

int main() {

  // The following 'factory' functions make spin operators with statically typed
  // indices and real coefficients.
  using static_indices::S_p;  // Spin-1/2 raising operator S_+
  using static_indices::S_m;  // Spin-1/2 lowering operator S_-
  using static_indices::S_z;  // Spin-1/2 operator S_z

  // Expression 'H' will represent the exchange interaction term.
  // Our spin operators will carry one integer index (site 1 or 2).
  auto H = 0.5 * (S_p(1) * S_m(2) + S_m(1) * S_p(2)) + S_z(1) * S_z(2);

  // Print 'H'
  std::cout << "H = " << H << std::endl;

  // Automatically analyze structure of 'H' and construct a 4-dimensional
  // Hilbert space (direct product of two spin-1/2 spaces).
  auto hs = make_hilbert_space(H);
  std::cout << "dim(hs) = " << hs.dim() << std::endl;

  // Construct a 'loperator' object that represents action of expression 'H' on
  // state vectors in the Hilbert space 'hs'.
  auto Hop = make_loperator(H, hs);

  // Here, we will act with 'Hop' on each of the 4 basis states |\psi> in 'hs',
  // |\phi> = Hop |\psi>, and print components of |\phi>. In other words,
  // we are going to construct the matrix representation <\phi|Hop|\psi>.

  // Preallocate state vectors.
  // Other containers, such as Eigen::Vector could be used instead.
  std::vector<double> phi(4), psi(4);
  // Iterate over basis states
  for(int i = 0; i < 4; ++i) {
    // Reset vectors |\phi> and |\psi> to zero.
    std::fill(phi.begin(), phi.end(), 0);
    std::fill(psi.begin(), psi.end(), 0);
    psi[i] = 1; // 'psi' is i-th basis vector now

    phi = Hop(psi);
    // NB.: It is generally recommended to use the in-place syntax
    //  Hop(psi, phi);
    // as it eliminates a memory allocation needed to store the result.

    // Print the result
    std::cout << "H|" << i << "> = ";
    for(int j = 0; j < 4; ++j) {
      std::cout << "+(" << phi[j] << ")|" << j << ">";
    }
    std::cout << std::endl;
  }

  return 0;
}
