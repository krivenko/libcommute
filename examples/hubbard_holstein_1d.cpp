/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

//
// Partial diagonalization of a one-dimensional Hubbard-Holstein model
// on a chain with L sites,
//
// \hat H =
//   - t \sum_{\sigma}\sum_i (c^\dagger_{i,\sigma} c_{i+1,\sigma} + h.c.)
//   - \mu \sum_i n_{i,\sigma}
//   + U \sum_i n_{i,\uparrow} n_{i,\downarrow}
//   + \omega \sum_i a^\dagger_i a_i
//   + g \sum_i (n_{i,\uparrow} + n_{i,\downarrow})(a^\dagger_i + a_i).
//
// Diagonalization is performed within the 2-electron sector by means of Eigen's
// 'SelfAdjointEigenSolver' and the 'n_fermion_sector_view' objects.
//

#include <libcommute/loperator/state_vector_eigen3.hpp>
#include <libcommute/libcommute.hpp>

#include <Eigen/Eigenvalues>

using namespace libcommute;

int main() {
  int const L = 4;       // Length of the chain
  int const b = 1;       // 2^b phonon states per chain site
  double const t = 0.5;  // Electron hopping constant
  double const mu = 1.0; // Chemical potential
  double const U = 2.0;  // Coulomb repulsion
  double const w = 2.0;  // Phonon frequency
  double const g = 0.1;  // Electron-phonon coupling

  // Expression of the Hubbard-Holstein Hamiltonian to be constructed.
  // Every creation and annihilation operator found in the expression carries
  // an integer index (chain site) and a string index (spin).
  static_indices::expr_real<int, std::string> H;

  // The following functions return operators with statically typed indices.
  using static_indices::c_dag; // Create an electron
  using static_indices::c;     // Destroy an electron
  using static_indices::n;     // Number of electrons
  using static_indices::a_dag; // Create a phonon
  using static_indices::a;     // Destroy a phonon

  // Electron hopping terms of H.
  for(auto spin : {"up", "down"}) {
    for(int i = 0; i < L - 1; ++i) {
      H += -t * (c_dag(i, spin) * c(i + 1, spin) + hc);
    }
  }
  // Chemical potential terms.
  for(int i = 0; i < L; ++i)
    H += -mu * (n(i, "up") + n(i, "down"));
  // Coulomb repulsion terms.
  for(int i = 0; i < L; ++i)
    H += U * n(i, "up") * n(i, "down");
  // Energy of the localized phonons.
  for(int i = 0; i < L; ++i)
    H += w * a_dag(i, "") * a(i, "");
  // Electron-phonon coupling.
  for(int i = 0; i < L; ++i)
    H += g * (n(i, "up") + n(i, "down")) * (a_dag(i, "") + a(i, ""));

  // Automatically analyse structure of 'H' and construct a Hilbert space.
  // Only the lowest 2^b states will be accounted for for each localized phonon.
  auto hs = make_hilbert_space(H, boson_es_constructor(b));
  std::cout << "Full Hilbert space dimension: " << hs.dim() << '\n';

  // Construct an 'loperator' object that represents action of the expression
  // 'H' on state vectors in the Hilbert space 'hs'.
  auto Hop = make_loperator(H, hs);

  // Diagonalize the 2-electron sector of the model.
  // using Eigen 3 and the 'n_fermion_sector_view' class.
  int const N_el = 2;

  auto sec_size = n_fermion_sector_size(hs, N_el);
  std::cout << "Size of the N_el = 2 sector: " << sec_size << '\n';

  // Preallocate a Hamiltonian matrix to be filled and diagonalized. Note that
  // we have to store matrix elements only within the small sector.
  Eigen::MatrixXd Hmat(sec_size, sec_size);

  // Preallocate an input state vector.
  Eigen::VectorXd psi(sec_size);

  // Create a constant 2-fermion sector view of the source vector 'psi'.
  // The view makes the small vector 'psi' compatible with 'Hop', which wants
  // to act in the full Hilbert space.
  auto psi_view = make_const_nfs_view(psi, hs, N_el);

  // Fill 'Hmat'.
  for(int k = 0; k < sec_size; ++k) {
    // Reset vector |\psi> to zero.
    psi = Eigen::VectorXd::Zero((int)sec_size);
    psi[k] = 1; // 'psi' is the k-th basis vector now.

    // Create a 2-fermion sector view of the k-th column of 'Hmat'.
    auto HMat_col_view = make_nfs_view(Hmat.col(k), hs, N_el);

    // Act with 'Hop' on the k-th sector basis state and write the result into
    // the k-th column of 'Hmat'.
    Hop(psi_view, HMat_col_view);
  }

  // Diagonalize 'Hmat'.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> diag(Hmat,
                                                      Eigen::EigenvaluesOnly);

  std::cout << "10 lowest eigenvalues of the N_el = " << N_el
            << " sector: " << diag.eigenvalues().head(10).transpose() << '\n';

  return 0;
}
