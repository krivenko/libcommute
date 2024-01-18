/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSl and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPl was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

//
// Fermi-Hubbard model on a large square lattice and diagonalization of
// its moderately sized sectors.
//
// \hat H = -\sum_{\langle i,j \rangle, \sigma}
//            t_{ij} c^\dagger_{i,\sigma} c_{j,\sigma}
//          -\mu \sum_{i, \sigma} n_{i,\sigma}
//          + U \sum_i n_{i,\uparrow} n_{i,\downarrow}
//

// 'state_vector_eigen3.hpp' must be included before 'libcommute.hpp'.
// This way libcommute will know how to deal with Eigen's vector types.
#include <libcommute/loperator/state_vector_eigen3.hpp>

#include <libcommute/libcommute.hpp>

// For Eigen::SelfAdjointEigenSolver
#include <Eigen/Eigenvalues>

using namespace libcommute;

int main() {

  //
  // Parameters of the system
  //

  // Linear sizes of the lattice
  int const Nx = 4;
  int const Ny = 4;

  // Hopping constant
  double const t = 0.5;
  // Chemical potential
  double const mu = 1.0;
  // Coulomb repulsion
  double const U = 2.0;

  // Fermi-Hubbard Hamiltonian to be constructed.
  // Every creation and annihilation operator met in the Hamiltonian
  // carry two integer indices (coordinates of a lattice site) and
  // one string index (spin).
  static_indices::expr_real<int, int, std::string> H;

  // The following functions make quantum operators with
  // statically typed indices.
  using static_indices::c_dag; // Create an electron
  using static_indices::c;     // Destroy an electron
  using static_indices::n;     // Number of electrons

  // Are two sites neighbors along the x-axis with periodicity?
  auto neighbors_x = [](int ix, int jx) {
    return std::abs(ix - jx) == 1 || std::abs(ix - jx) == Nx - 1;
  };
  // Are two sites neighbors along the y-axis with periodicity?
  auto neighbors_y = [](int iy, int jy) {
    return std::abs(iy - jy) == 1 || std::abs(iy - jy) == Ny - 1;
  };

  // Hopping terms of H
  for(auto spin : {"up", "down"}) {
    for(int ix = 0; ix < Nx; ++ix) {
      for(int iy = 0; iy < Ny; ++iy) {
        for(int jx = 0; jx < Nx; ++jx) {
          for(int jy = 0; jy < Ny; ++jy) {
            // Skip all pairs of lattice sites (ix,iy) and (jx,jy)
            // that are not nearest neighbors.
            if((neighbors_x(ix, jx) && iy == jy) ||
               (ix == jx && neighbors_y(iy, jy))) {
              // Add a hopping term
              H += -t * c_dag(ix, iy, spin) * c(jx, jy, spin);
            }
          }
        }
      }
    }
  }

  // Chemical potential terms
  for(int ix = 0; ix < Nx; ++ix)
    for(int iy = 0; iy < Ny; ++iy) {
      H += -mu * (n(ix, iy, "up") + n(ix, iy, "down"));
    }

  // Coulomb repulsion terms
  for(int ix = 0; ix < Nx; ++ix)
    for(int iy = 0; iy < Ny; ++iy) {
      H += U * n(ix, iy, "up") * n(ix, iy, "down");
    }

  // Automatically analyze structure of 'H' and construct a Hilbert space
  auto hs = make_hilbert_space(H);
  std::cout << "Full Hilbert space dimension: " << hs.dim() << std::endl;

  // Construct a 'loperator' object that represents action of expression 'H' on
  // state vectors in the Hilbert space 'hs'.
  auto Hop = make_loperator(H, hs);

  //
  // Diagonalize the N = 2 sector of the model using Eigen 3 and the
  // 'n_fermion_sector_view' class.
  //
  {
    int const N = 2;

    auto const sector_size = n_fermion_sector_size(hs, N);
    std::cout << "Size of the N = 2 sector: " << sector_size << std::endl;

    // Preallocate a Hamiltonian matrix to be filled and diagonalized.
    // Note that we have to store matrix elements only within the small sector.
    Eigen::MatrixXd Hmat(sector_size, sector_size);

    // Preallocate an input state vector
    Eigen::VectorXd psi(sector_size);

    // Create a constant 2-fermion sector view of the source vector 'psi'
    auto psi_view = make_const_nfs_view(psi, hs, N);

    // Fill 'Hmat'
    for(int i = 0; i < int(sector_size); ++i) {
      // Reset vector |\psi> to zero.
      psi = Eigen::VectorXd::Zero(int(sector_size));
      psi[i] = 1; // 'psi' is the i-th basis vector now

      // Create a 2-fermion sector view of the i-th column of 'Hmat'.
      auto HMat_col_view = make_nfs_view(Hmat.col(i), hs, N);

      // Act with Hop on the i-th sector basis state and write the result into
      // the i-th column of Hmat.
      Hop(psi_view, HMat_col_view);
    }

    // Diagonalize 'Hmat'
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> diag(Hmat,
                                                        Eigen::EigenvaluesOnly);

    std::cout << "10 lowest eigenvalues of the N = 2 sector: "
              << diag.eigenvalues().head(10).transpose() << std::endl;
  }

  //
  // Diagonalize the N_up = 1, N_down = 1 multisector of the model using
  // Eigen 3 and the 'n_fermion_multisector_view' class.
  //
  {
    int const N_up = 1, N_down = 1;

    // Define sectors contributing to the multisector
    sector_descriptor<decltype(hs)> sector_up{{}, N_up};
    sector_descriptor<decltype(hs)> sector_down{{}, N_down};

    // Fill index sets {S_up} and {S_down}
    for(int ix = 0; ix < Nx; ++ix) {
      for(int iy = 0; iy < Ny; ++iy) {
        sector_up.indices.emplace(ix, iy, "up");
        sector_down.indices.emplace(ix, iy, "down");
      }
    }

    std::vector<sector_descriptor<decltype(hs)>> sectors = {sector_up,
                                                            sector_down};

    auto const multisector_size = n_fermion_multisector_size(hs, sectors);
    std::cout << "Size of the N_up = 1, N_down = 1 multisector: "
              << multisector_size << std::endl;

    // Preallocate a Hamiltonian matrix to be filled and diagonalized.
    // Note that we have to store matrix elements only within the multisector.
    Eigen::MatrixXd Hmat(multisector_size, multisector_size);

    // Preallocate an input state vector
    Eigen::VectorXd psi(multisector_size);

    // Create a constant multisector view of the source vector 'psi'
    auto psi_view = make_const_nfms_view(psi, hs, sectors);

    // Fill 'Hmat'
    for(int i = 0; i < int(multisector_size); ++i) {
      // Reset vector |\psi> to zero.
      psi = Eigen::VectorXd::Zero(int(multisector_size));
      psi[i] = 1; // 'psi' is the i-th basis vector now

      // Create a multisector view of the i-th column of 'Hmat'.
      auto HMat_col_view = make_nfms_view(Hmat.col(i), hs, sectors);

      // Act with Hop on the i-th multisector basis state and write the result
      // into the i-th column of Hmat.
      Hop(psi_view, HMat_col_view);
    }

    // Diagonalize 'Hmat'
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> diag(Hmat,
                                                        Eigen::EigenvaluesOnly);

    std::cout << "10 lowest eigenvalues of the multisector: "
              << diag.eigenvalues().head(10).transpose() << std::endl;
  }

  return 0;
}
