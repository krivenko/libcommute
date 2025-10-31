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
// Sparse Hilbert spaces and the 'compressed_state_view' functionality
// demonstrated by the example of the spin-1 Heisenberg ring.
//
// The Hamiltonian of the model is
//  \hat H = \sum_{i=0}^{N-1} \mathbf{S}_i \cdot \mathbf{S}_{(i+1)\bmod N} =
//           \sum_{\substack{i=0,\\j = (i+1)\bmod N}}^{N-1} S_{z,i} S_{z,j} +
//              \frac{1}{2} \left[ S_{+,i} S_{-,j} + S_{-,i} S_{+,j} \right].
//

// 'state_vector_eigen3.hpp' must be included before 'libcommute.hpp'.
// This way libcommute will know how to deal with Eigen's vector types.
#include <libcommute/loperator/state_vector_eigen3.hpp>

#include <libcommute/libcommute.hpp>

#include <Eigen/Eigenvalues>

using namespace libcommute;

int main() {
  int const N = 6; // Number of spins

  // Expression of the Heisenberg Hamiltonian to be constructed.
  // Every creation and annihilation operator found in the expression carries
  // one integer index (ring site).
  static_indices::expr_real<int> H;

  // The following functions return operators with statically typed indices.
  using static_indices::S_z; // Spin projection operator S_z
  using static_indices::S_p; // Spin-raising operator S_+
  using static_indices::S_m; // Spin-lowering operator S_-

  // Heisenberg interaction of the S=1 spins
  for(int i = 0; i < N; ++i) {
    int j = (i + 1) % N;
    H += S_z<3>(i) * S_z<3>(j) +
         0.5 * (S_p<3>(i) * S_m<3>(j) + S_m<3>(i) * S_p<3>(j));
  }

  // Automatically analyse structure of 'H' and construct a Hilbert space.
  auto hs = make_hilbert_space(H);

  // The Hilbert space 'hs' is sparse: The minimal required container size for
  // state vectors exceeds the dimension of the space.
  std::cout << "Is the Hilbert space sparse? " << hs.is_sparse() << '\n';
  std::cout << "Hilbert space dimension: " << hs.dim() << '\n';
  std::cout << "Minimal container size for state vectors: " << hs.vec_size()
            << '\n';

  auto const dim = hs.dim(); // 3^N

  // Construct an 'loperator' object that represents action of the expression
  // 'H' on state vectors in the Hilbert space 'hs'.
  auto Hop = make_loperator(H, hs);

  // Preallocate an input state vector of the length 3^N.
  // This container is compressed, i.e. it stores quantum amplitudes of
  // the *physical* basis states only!
  Eigen::VectorXd psi(dim);

  // Preallocate a Hamiltonian matrix to be filled and diagonalized.
  // This matrix is also compressed.
  Eigen::MatrixXd Hmat(dim, dim);

  // Create a view of the compressed state 'psi'.
  auto psi_view = make_const_comp_state_view(psi, hs);

  // Fill 'Hmat'
  // Important: k takes on values corresponding to the *physical* basis states
  // only.
  foreach(hs, [dim, &hs, &Hop, &psi, &psi_view, &Hmat](sv_index_type k) {
    // Translate k into the compressed (continuous) range [0; 3^N-1]
    auto const k_comp = psi_view.map_index(k);

    // Reset vector |\psi> to zero.
    psi = Eigen::VectorXd::Zero((int)dim);

    // Now, 'psi' is the basis vector corresponding to 'k'.
    psi((int)k_comp) = 1;

    // Create a view of the column of 'Hmat' that corresponds to 'k'.
    auto HMat_col_view = make_comp_state_view(Hmat.col((int)k_comp), hs);

    // Act with 'Hop' on the basis state 'k' and write the result into
    // the column of 'Hmat' that also corresponds to 'k'.
    Hop(psi_view, HMat_col_view);
  });

  // Diagonalize 'Hmat'.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> diag(Hmat,
                                                      Eigen::EigenvaluesOnly);

  std::cout << "10 lowest eigenvalues of H = "
            << diag.eigenvalues().head(10).transpose() << '\n';

  return 0;
}
