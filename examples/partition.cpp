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
// Hilbert space partition into disjoint subspaces
//

#include <libcommute/expression/factories.hpp>
#include <libcommute/loperator/mapped_basis_view.hpp>
#include <libcommute/loperator/space_partition.hpp>

#include <iostream>
#include <vector>

using namespace libcommute;

int main() {

  using namespace static_indices; // For c(), c_dag() and n()

  //
  // Build Hamiltonian of the 3-orbital Hubbard-Kanamori atom
  //

  int const n_orbs = 3;
  double const mu = 0.7; // Chemical
  double const U = 3.0;  // Coulomb repulsion
  double const J = 0.3;  // Hund interaction

  // Hamiltonian
  expression<double, std::string, int> H;

  // Chemical potential terms
  for(int o = 0; o < n_orbs; ++o) {
    H += -mu * (n("up", o) + n("dn", o));
  }

  // Intraorbital interactions
  for(int o = 0; o < n_orbs; ++o) {
    H += U * n("up", o) * n("dn", o);
  }
  // Interorbital interactions, different spins
  for(int o1 = 0; o1 < n_orbs; ++o1) {
    for(int o2 = 0; o2 < n_orbs; ++o2) {
      if(o1 == o2) continue;
      H += (U - 2 * J) * n("up", o1) * n("dn", o2);
    }
  }
  // Interorbital interactions, equal spins
  for(int o1 = 0; o1 < n_orbs; ++o1) {
    for(int o2 = 0; o2 < n_orbs; ++o2) {
      if(o2 >= o1) continue;
      H += (U - 3 * J) * n("up", o1) * n("up", o2);
      H += (U - 3 * J) * n("dn", o1) * n("dn", o2);
    }
  }
  // Spin flip and pair hopping terms
  for(int o1 = 0; o1 < n_orbs; ++o1) {
    for(int o2 = 0; o2 < n_orbs; ++o2) {
      if(o1 == o2) continue;
      H += -J * c_dag("up", o1) * c_dag("dn", o1) * c("up", o2) * c("dn", o2);
      H += -J * c_dag("up", o1) * c_dag("dn", o2) * c("up", o2) * c("dn", o1);
    }
  }

  //
  // Hilbert space of the problem
  //

  auto hs = make_hilbert_space(H);

  //
  // Linear operator form of the Hamiltonian
  //

  auto Hop = make_loperator(H, hs);

  //
  // Reveal invariant subspaces of H
  //

  auto sp1 = space_partition(Hop, hs);

  std::cout << "Total dimension of the Hilbert space is " << sp1.dim() << '\n';

  std::cout << "H has " << sp1.n_subspaces() << " invariant subspaces\n";

  //
  // Once again, this time saving the matrix elements of H
  //

  matrix_elements_map<double> H_elements;
  auto sp2 = space_partition(Hop, hs, H_elements);

  std::cout << "H has " << H_elements.size() << " non-vanishing matrix elements"
            << '\n';

  //
  // Now merge some invariant subspaces to make sure that all electron
  // creation and annihilation operators connect one subspace to one subspace.
  //

  for(std::string spin : {"up", "dn"}) {
    for(int o = 0; o < n_orbs; ++o) {
      auto Cdagop = make_loperator(c_dag(spin, o), hs);
      auto Cop = make_loperator(c(spin, o), hs);

      auto matrix_elements = sp2.merge_subspaces(Cdagop, Cop, hs);

      std::cout << c_dag(spin, o) << " has " << matrix_elements.first.size()
                << " non-vanishing matrix elements\n";
      std::cout << c(spin, o) << " has " << matrix_elements.first.size()
                << " non-vanishing matrix elements\n";
    }
  }

  std::cout << "Number of disjoint subspaces after merging is "
            << sp2.n_subspaces() << '\n';

  //
  // foreach()
  //

  std::cout << "Distribution of basis states over the subspaces:\n";
  foreach(sp2, [](sv_index_type index, sv_index_type subspace) {
    std::cout << index << " => " << subspace << '\n';
  });

  //
  // Bonus: Using mapped_basis_view to act on a vector in one of the invariant
  // subspaces.
  //

  // Collect all basis state indices from one subspace
  std::vector<sv_index_type> basis_states_in_subspace_24 =
      sp2.subspace_basis(24);

  auto sp24_dim = basis_states_in_subspace_24.size();
  std::cout << "Subspace 24 is " << sp24_dim << "-dimensional\n";

  // Make a mapper object
  basis_mapper mapper(basis_states_in_subspace_24);

  // 'psi' and 'phi' are small vectors entirely lying in the 24-th subspace
  std::vector<double> psi(sp24_dim), phi(sp24_dim);

  //
  // These views have dimension of the full Hilbert space!
  //
  auto psi_view = mapper.make_const_view(psi);
  auto phi_view = mapper.make_view(phi);

  //
  // Act on 'psi_view' as if it was a vector in the full space
  // No exception here because 'Hop' connects subspace 24 only to itself.
  //
  Hop(psi_view, phi_view);

  return 0;
}
