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

#include <catch.hpp>

#include <libcommute/expression/factories.hpp>
#include <libcommute/loperator/space_partition.hpp>
#include <libcommute/loperator/sparse_state_vector.hpp>

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

using namespace libcommute;

using subspaces_t = std::set<std::set<sv_index_type>>;

template <typename HSType>
subspaces_t collect_subspaces(space_partition<HSType> const& sp) {
  // Sets are used to neglect order of subspaces and of states within a subspace
  std::vector<std::set<sv_index_type>> v_cl(sp.n_subspaces());
  foreach(sp, [&](int i, int subspace) { v_cl[subspace].insert(i); });
  subspaces_t cl{v_cl.cbegin(), v_cl.cend()};
  return cl;
}

TEST_CASE("Automatic Hilbert space partition", "[space_partition]") {
  using namespace static_indices;

  // 3 orbital Hubbard-Kanamori atom
  int const n_orbs = 3;
  double const mu = 0.7;
  double const U = 3.0;
  double const J = 0.3;

  // Hamiltonian
  expression<double, std::string, int> H;

  // Chemical potential
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
  // Spin flips and pair hoppings
  for(int o1 = 0; o1 < n_orbs; ++o1) {
    for(int o2 = 0; o2 < n_orbs; ++o2) {
      if(o1 == o2) continue;
      H += -J * c_dag("up", o1) * c_dag("dn", o1) * c("up", o2) * c("dn", o2);
      H += -J * c_dag("up", o1) * c_dag("dn", o2) * c("up", o2) * c("dn", o1);
    }
  }

  // Hilbert space
  auto hs = make_hilbert_space(H);
  // Linear operator form of the Hamiltonian
  auto Hop = make_loperator(H, hs);

  using es_type = elementary_space_fermion<std::string, int>;
  sv_index_type d0 = 1 << hs.bit_range(es_type("dn", 0)).first;
  sv_index_type d1 = 1 << hs.bit_range(es_type("dn", 1)).first;
  sv_index_type d2 = 1 << hs.bit_range(es_type("dn", 2)).first;
  sv_index_type u0 = 1 << hs.bit_range(es_type("up", 0)).first;
  sv_index_type u1 = 1 << hs.bit_range(es_type("up", 1)).first;
  sv_index_type u2 = 1 << hs.bit_range(es_type("up", 2)).first;

  SECTION("space_partition()") {
    auto sp = make_space_partition(Hop, hs);

    CHECK(sp.dim() == 64);
    CHECK(sp.n_subspaces() == 44);

    // Expected classification of states
    subspaces_t ref_cl{
        // N=0
        {0},
        // N=1
        {d0},
        {d1},
        {d2},
        {u0},
        {u1},
        {u2},
        // N=2, same spin
        {d0 + d1},
        {d0 + d2},
        {d1 + d2},
        {u0 + u1},
        {u0 + u2},
        {u1 + u2},
        // N=2, pair hopping
        {d0 + u0, d1 + u1, d2 + u2},
        // N=2, spin flip
        {d0 + u1, d1 + u0},
        {d0 + u2, d2 + u0},
        {d1 + u2, d2 + u1},
        // N=3
        {d0 + d1 + d2},
        {u0 + u1 + u2},
        {d0 + d1 + u0, d1 + d2 + u2},
        {d0 + d2 + u0, d1 + d2 + u1},
        {d0 + d1 + u1, d0 + d2 + u2},
        {d0 + u0 + u1, d2 + u1 + u2},
        {d1 + u0 + u1, d2 + u0 + u2},
        {d0 + u0 + u2, d1 + u1 + u2},
        {d1 + d2 + u0, d0 + d2 + u1, d0 + d1 + u2},
        {d2 + u0 + u1, d0 + u1 + u2, d1 + u0 + u2},
        // N=4, 2 holes with the same spin
        {d2 + u0 + u1 + u2},
        {d1 + u0 + u1 + u2},
        {d0 + u0 + u1 + u2},
        {d0 + d1 + d2 + u2},
        {d0 + d1 + d2 + u1},
        {d0 + d1 + d2 + u0},
        // N=4, pair hopping
        {d1 + d2 + u1 + u2, d0 + d2 + u0 + u2, d0 + d1 + u0 + u1},
        // N=4, spin flip
        {d1 + d2 + u0 + u2, d0 + d2 + u1 + u2},
        {d1 + d2 + u0 + u1, d0 + d1 + u1 + u2},
        {d0 + d2 + u0 + u1, d0 + d1 + u0 + u2},
        // N=5
        {d1 + d2 + u0 + u1 + u2},
        {d0 + d2 + u0 + u1 + u2},
        {d0 + d1 + u0 + u1 + u2},
        {d0 + d1 + d2 + u1 + u2},
        {d0 + d1 + d2 + u0 + u2},
        {d0 + d1 + d2 + u0 + u1},
        // N=6
        {d0 + d1 + d2 + u0 + u1 + u2}};

    CHECK(collect_subspaces(sp) == ref_cl);

    SECTION("subspace_bases()") {
      auto bases = sp.subspace_bases();
      CHECK(bases.size() == sp.n_subspaces());
      for(auto const& basis : bases) {
        std::set<sv_index_type> basis_set{basis.cbegin(), basis.cend()};
        CHECK(ref_cl.count(basis_set) == 1);
      }
    }

    SECTION("subspace_basis()") {
      for(int subspace = 0; subspace < (int)sp.n_subspaces(); ++subspace) {
        auto basis = sp.subspace_basis(subspace);
        std::set<sv_index_type> basis_set{basis.cbegin(), basis.cend()};
        CHECK(ref_cl.count(basis_set) == 1);
      }
      CHECK_THROWS_AS(sp.subspace_basis(sp.n_subspaces()), std::runtime_error);
    }
  }

  SECTION("Matrix elements") {
    matrix_elements_map<double> matrix_elements;
    auto sp = make_space_partition(Hop, hs, matrix_elements);

    struct melem_t {
      sv_index_type from;
      sv_index_type to;
      double val;
      bool operator<(melem_t const& me) const {
        return ((from != me.from) ? from < me.from : to < me.to);
      }
      bool operator==(melem_t const& me) const {
        return from == me.from && to == me.to && std::abs(val - me.val) < 1e-10;
      }
    };
    using melem_set_t = std::set<std::set<melem_t>>;

    std::vector<std::set<melem_t>> v_melem(sp.n_subspaces());
    for(auto const& x : matrix_elements) {
      v_melem[sp[x.first.first]].insert(
          melem_t{x.first.second, x.first.first, x.second});
    }
    melem_set_t melem(v_melem.cbegin(), v_melem.cend());

    // Expected matrix elements
    melem_set_t ref_melem{
        // N=0
        std::set<melem_t>{},
        // N=1
        {{d0, d0, -mu}},
        {{d1, d1, -mu}},
        {{d2, d2, -mu}},
        {{u0, u0, -mu}},
        {{u1, u1, -mu}},
        {{u2, u2, -mu}},
        // N=2, same spin
        {{d0 + d1, d0 + d1, -2 * mu + U - 3 * J}},
        {{d0 + d2, d0 + d2, -2 * mu + U - 3 * J}},
        {{d1 + d2, d1 + d2, -2 * mu + U - 3 * J}},
        {{u0 + u1, u0 + u1, -2 * mu + U - 3 * J}},
        {{u0 + u2, u0 + u2, -2 * mu + U - 3 * J}},
        {{u1 + u2, u1 + u2, -2 * mu + U - 3 * J}},
        // N=2, pair hopping
        {{d0 + u0, d0 + u0, -2 * mu + U},
         {d1 + u1, d1 + u1, -2 * mu + U},
         {d2 + u2, d2 + u2, -2 * mu + U},
         {d0 + u0, d1 + u1, J},
         {d0 + u0, d2 + u2, J},
         {d1 + u1, d2 + u2, J},
         {d1 + u1, d0 + u0, J},
         {d2 + u2, d0 + u0, J},
         {d2 + u2, d1 + u1, J}},
        // N=2, spin flip
        {{d0 + u1, d0 + u1, -2 * mu + U - 2 * J},
         {d1 + u0, d1 + u0, -2 * mu + U - 2 * J},
         {d0 + u1, d1 + u0, J},
         {d1 + u0, d0 + u1, J}},
        {{d0 + u2, d0 + u2, -2 * mu + U - 2 * J},
         {d2 + u0, d2 + u0, -2 * mu + U - 2 * J},
         {d0 + u2, d2 + u0, J},
         {d2 + u0, d0 + u2, J}},
        {{d1 + u2, d1 + u2, -2 * mu + U - 2 * J},
         {d2 + u1, d2 + u1, -2 * mu + U - 2 * J},
         {d1 + u2, d2 + u1, J},
         {d2 + u1, d1 + u2, J}},
        // N=3
        {{d0 + d1 + d2, d0 + d1 + d2, -3 * mu + 3 * U - 9 * J}},
        {{u0 + u1 + u2, u0 + u1 + u2, -3 * mu + 3 * U - 9 * J}},
        {{d0 + d1 + u0, d0 + d1 + u0, -3 * mu + 3 * U - 5 * J},
         {d1 + d2 + u2, d1 + d2 + u2, -3 * mu + 3 * U - 5 * J},
         {d0 + d1 + u0, d1 + d2 + u2, -J},
         {d1 + d2 + u2, d0 + d1 + u0, -J}},
        {{d0 + d2 + u0, d0 + d2 + u0, -3 * mu + 3 * U - 5 * J},
         {d1 + d2 + u1, d1 + d2 + u1, -3 * mu + 3 * U - 5 * J},
         {d0 + d2 + u0, d1 + d2 + u1, J},
         {d1 + d2 + u1, d0 + d2 + u0, J}},
        {{d0 + d1 + u1, d0 + d1 + u1, -3 * mu + 3 * U - 5 * J},
         {d0 + d2 + u2, d0 + d2 + u2, -3 * mu + 3 * U - 5 * J},
         {d0 + d1 + u1, d0 + d2 + u2, J},
         {d0 + d2 + u2, d0 + d1 + u1, J}},
        {{d0 + u0 + u1, d0 + u0 + u1, -3 * mu + 3 * U - 5 * J},
         {d2 + u1 + u2, d2 + u1 + u2, -3 * mu + 3 * U - 5 * J},
         {d0 + u0 + u1, d2 + u1 + u2, -J},
         {d2 + u1 + u2, d0 + u0 + u1, -J}},
        {{d1 + u0 + u1, d1 + u0 + u1, -3 * mu + 3 * U - 5 * J},
         {d2 + u0 + u2, d2 + u0 + u2, -3 * mu + 3 * U - 5 * J},
         {d1 + u0 + u1, d2 + u0 + u2, J},
         {d2 + u0 + u2, d1 + u0 + u1, J}},
        {{d0 + u0 + u2, d0 + u0 + u2, -3 * mu + 3 * U - 5 * J},
         {d1 + u1 + u2, d1 + u1 + u2, -3 * mu + 3 * U - 5 * J},
         {d0 + u0 + u2, d1 + u1 + u2, J},
         {d1 + u1 + u2, d0 + u0 + u2, J}},
        {{d1 + d2 + u0, d1 + d2 + u0, -3 * mu + 3 * U - 7 * J},
         {d0 + d2 + u1, d0 + d2 + u1, -3 * mu + 3 * U - 7 * J},
         {d0 + d1 + u2, d0 + d1 + u2, -3 * mu + 3 * U - 7 * J},
         {d1 + d2 + u0, d0 + d2 + u1, J},
         {d0 + d2 + u1, d1 + d2 + u0, J},
         {d1 + d2 + u0, d0 + d1 + u2, -J},
         {d0 + d1 + u2, d1 + d2 + u0, -J},
         {d0 + d2 + u1, d0 + d1 + u2, J},
         {d0 + d1 + u2, d0 + d2 + u1, J}},
        {{d2 + u0 + u1, d2 + u0 + u1, -3 * mu + 3 * U - 7 * J},
         {d0 + u1 + u2, d0 + u1 + u2, -3 * mu + 3 * U - 7 * J},
         {d1 + u0 + u2, d1 + u0 + u2, -3 * mu + 3 * U - 7 * J},
         {d2 + u0 + u1, d0 + u1 + u2, -J},
         {d0 + u1 + u2, d2 + u0 + u1, -J},
         {d2 + u0 + u1, d1 + u0 + u2, J},
         {d1 + u0 + u2, d2 + u0 + u1, J},
         {d0 + u1 + u2, d1 + u0 + u2, J},
         {d1 + u0 + u2, d0 + u1 + u2, J}},
        // N=4, 2 holes with the same spin
        {{d2 + u0 + u1 + u2, d2 + u0 + u1 + u2, -4 * mu + 6 * U - 13 * J}},
        {{d1 + u0 + u1 + u2, d1 + u0 + u1 + u2, -4 * mu + 6 * U - 13 * J}},
        {{d0 + u0 + u1 + u2, d0 + u0 + u1 + u2, -4 * mu + 6 * U - 13 * J}},
        {{d0 + d1 + d2 + u0, d0 + d1 + d2 + u0, -4 * mu + 6 * U - 13 * J}},
        {{d0 + d1 + d2 + u1, d0 + d1 + d2 + u1, -4 * mu + 6 * U - 13 * J}},
        {{d0 + d1 + d2 + u2, d0 + d1 + d2 + u2, -4 * mu + 6 * U - 13 * J}},
        // N=4, pair hopping
        {{d1 + d2 + u1 + u2, d1 + d2 + u1 + u2, -4 * mu + 6 * U - 10 * J},
         {d0 + d2 + u0 + u2, d0 + d2 + u0 + u2, -4 * mu + 6 * U - 10 * J},
         {d0 + d1 + u0 + u1, d0 + d1 + u0 + u1, -4 * mu + 6 * U - 10 * J},
         {d1 + d2 + u1 + u2, d0 + d2 + u0 + u2, J},
         {d0 + d2 + u0 + u2, d1 + d2 + u1 + u2, J},
         {d1 + d2 + u1 + u2, d0 + d1 + u0 + u1, J},
         {d0 + d1 + u0 + u1, d1 + d2 + u1 + u2, J},
         {d0 + d2 + u0 + u2, d0 + d1 + u0 + u1, J},
         {d0 + d1 + u0 + u1, d0 + d2 + u0 + u2, J}},
        // N=4, spin flip
        {{d1 + d2 + u0 + u2, d1 + d2 + u0 + u2, -4 * mu + 6 * U - 12 * J},
         {d0 + d2 + u1 + u2, d0 + d2 + u1 + u2, -4 * mu + 6 * U - 12 * J},
         {d1 + d2 + u0 + u2, d0 + d2 + u1 + u2, J},
         {d0 + d2 + u1 + u2, d1 + d2 + u0 + u2, J}},
        {{d1 + d2 + u0 + u1, d1 + d2 + u0 + u1, -4 * mu + 6 * U - 12 * J},
         {d0 + d1 + u1 + u2, d0 + d1 + u1 + u2, -4 * mu + 6 * U - 12 * J},
         {d1 + d2 + u0 + u1, d0 + d1 + u1 + u2, J},
         {d0 + d1 + u1 + u2, d1 + d2 + u0 + u1, J}},
        {{d0 + d2 + u0 + u1, d0 + d2 + u0 + u1, -4 * mu + 6 * U - 12 * J},
         {d0 + d1 + u0 + u2, d0 + d1 + u0 + u2, -4 * mu + 6 * U - 12 * J},
         {d0 + d2 + u0 + u1, d0 + d1 + u0 + u2, J},
         {d0 + d1 + u0 + u2, d0 + d2 + u0 + u1, J}},
        // N=5
        {{d1 + d2 + u0 + u1 + u2,
          d1 + d2 + u0 + u1 + u2,
          -5 * mu + 10 * U - 20 * J}},
        {{d0 + d2 + u0 + u1 + u2,
          d0 + d2 + u0 + u1 + u2,
          -5 * mu + 10 * U - 20 * J}},
        {{d0 + d1 + u0 + u1 + u2,
          d0 + d1 + u0 + u1 + u2,
          -5 * mu + 10 * U - 20 * J}},
        {{d0 + d1 + d2 + u1 + u2,
          d0 + d1 + d2 + u1 + u2,
          -5 * mu + 10 * U - 20 * J}},
        {{d0 + d1 + d2 + u0 + u2,
          d0 + d1 + d2 + u0 + u2,
          -5 * mu + 10 * U - 20 * J}},
        {{d0 + d1 + d2 + u0 + u1,
          d0 + d1 + d2 + u0 + u1,
          -5 * mu + 10 * U - 20 * J}},
        // N=6
        {{d0 + d1 + d2 + u0 + u1 + u2,
          d0 + d1 + d2 + u0 + u1 + u2,
          -6 * mu + 15 * U - 30 * J}}};

    CHECK(melem == ref_melem);
  }

  SECTION("merge_subspaces()") {

    auto check_one_to_one = [&hs](space_partition<decltype(hs)> const& sp,
                                  std::vector<decltype(Hop)> const& ops) {
      // Calculated classification of states
      auto cl = collect_subspaces(sp);

      sparse_state_vector<double> in_state(hs.vec_size());
      sparse_state_vector<double> out_state(hs.vec_size());

      for(auto const& op : ops) {
        for(auto const& i_sp : cl) {
          std::set<sv_index_type> f_sp;
          for(auto i : i_sp) {
            in_state[i] = 1.0;
            op(in_state, out_state);
            foreach(out_state, [&f_sp](sv_index_type f, double a) {
              if(std::abs(a) < 1e-10) return;
              f_sp.insert(f);
            });
            set_zeros(in_state);
          }

          // op maps i_sp onto zero
          if(f_sp.size() == 0) continue;

          // Check if op maps i_sp to only one subspace
          auto n =
              std::count_if(cl.begin(),
                            cl.end(),
                            [&f_sp](std::set<sv_index_type> const& f_sp_ref) {
                              return std::includes(f_sp_ref.cbegin(),
                                                   f_sp_ref.cend(),
                                                   f_sp.cbegin(),
                                                   f_sp.cend());
                            });
          CHECK(n == 1);
        }
      }
    };

    SECTION("Hubbard-Kanamori") {
      auto sp = make_space_partition(Hop, hs);

      std::vector<decltype(Hop)> ops1, ops2, all_ops;
      for(int o = 0; o < n_orbs; ++o) {
        ops1.emplace_back(c_dag("up", o) + c("dn", o), hs);
        ops2.emplace_back(c("up", o) + c_dag("dn", o), hs);

        all_ops.emplace_back(ops1.back());
        all_ops.emplace_back(ops2.back());

        sp.merge_subspaces(ops1.back(), ops2.back());
      }

      check_one_to_one(sp, all_ops);
    }

    SECTION("4 bath orbitals") {
      int const n_bath_orbs = 4;
      std::vector<double> const eps = {-0.2, -0.1, 0.1, 0.2};
      double const V = 0.7;

      for(std::string spin : {"dn", "up"}) {
        for(int o = 0; o < n_bath_orbs; ++o) {
          H += eps[o] * n(spin, n_orbs + o);
        }
        for(int o1 = 0; o1 < n_orbs; ++o1) {
          for(int o2 = 0; o2 < n_bath_orbs; ++o2) {
            H += V * (c_dag(spin, o1) * c(spin, n_orbs + o2) +
                      c_dag(spin, n_orbs + o2) * c(spin, o1));
          }
        }
      }

      // Hilbert space
      auto hs_bath = make_hilbert_space(H);
      // Linear operator form of the Hamiltonian
      auto H_bath_op = make_loperator(H, hs_bath);

      auto sp = make_space_partition(H_bath_op, hs_bath);

      std::vector<decltype(Hop)> Cd, C, all_ops;
      for(std::string spin : {"dn", "up"}) {
        for(int o = 0; o < n_orbs + n_bath_orbs; ++o) {
          Cd.emplace_back(c_dag(spin, o), hs_bath);
          C.emplace_back(c(spin, o), hs_bath);

          all_ops.emplace_back(Cd.back());
          all_ops.emplace_back(C.back());

          sp.merge_subspaces(Cd.back(), C.back());
        }
      }

      check_one_to_one(sp, all_ops);
    }
  }

  SECTION("find_connections") {

    std::vector<expression<double, std::string, int>> expr = {
        {},
        c_dag("up", 0) * c("dn", 1),
        c_dag("up", 0) * c("dn", 1) + c_dag("dn", 0) * c("dn", 1),
        c_dag("dn", 1) * c_dag("up", 1),
        n("up", 2)};

    auto sp = make_space_partition(Hop, hs);
    auto op = make_loperator(expr[1] + expr[2] + expr[3] + expr[4], hs);
    auto conns = sp.find_connections(op);

    connections_map conns_ref;
    std::vector<double> in_state(hs.vec_size());
    foreach(hs, [&](sv_index_type i) {
      in_state[i] = 1.0;
      auto out_state = op(in_state);
      foreach(out_state, [&](sv_index_type f, double a) {
        if(std::abs(a) < 1e-10) return;
        conns_ref.insert(std::make_pair(sp[i], sp[f]));
      });
      in_state[i] = 0;
    });

    CHECK(conns == conns_ref);

    for(auto const& expr1 : expr) {
      for(auto const& expr2 : expr) {
        auto conns1 = sp.find_connections(make_loperator(expr1, hs));
        auto conns2 = sp.find_connections(make_loperator(expr2, hs));
        auto conns12 =
            sp.find_connections(make_loperator(expr1 + expr2, hs));

        connections_map conns12_ref;
        std::set_union(conns1.begin(),
                       conns1.end(),
                       conns2.begin(),
                       conns2.end(),
                       std::inserter(conns12_ref, conns12_ref.end()));

        CHECK(conns12 == conns12_ref);
      }
    }
  }

  SECTION("Sparse Hilbert space") {
    auto H_sp = (S_p<3>(0) * S_m<3>(1) + S_m<3>(0) * S_p<3>(1)) +
                (S_p<3>(1) * S_m<3>(2) + S_m<3>(1) * S_p<3>(2)) +
                (S_p<3>(2) * S_m<3>(0) + S_m<3>(2) * S_p<3>(0));
    auto hs_sp = make_hilbert_space(H_sp);
    REQUIRE(hs_sp.is_sparse());
    auto H_sp_op = make_loperator(H_sp, hs_sp);
    auto sp = make_space_partition(H_sp_op, hs_sp);

    CHECK(sp.dim() == 27);
    CHECK(sp.n_subspaces() == 7); // S=3 septuplet

    using es_type = elementary_space_spin<int>;
    sv_index_type S0_m1 = 0 << hs_sp.bit_range(es_type(1.0, 0)).first;
    sv_index_type S0_0 = 1 << hs_sp.bit_range(es_type(1.0, 0)).first;
    sv_index_type S0_p1 = 2 << hs_sp.bit_range(es_type(1.0, 0)).first;
    sv_index_type S1_m1 = 0 << hs_sp.bit_range(es_type(1.0, 1)).first;
    sv_index_type S1_0 = 1 << hs_sp.bit_range(es_type(1.0, 1)).first;
    sv_index_type S1_p1 = 2 << hs_sp.bit_range(es_type(1.0, 1)).first;
    sv_index_type S2_m1 = 0 << hs_sp.bit_range(es_type(1.0, 2)).first;
    sv_index_type S2_0 = 1 << hs_sp.bit_range(es_type(1.0, 2)).first;
    sv_index_type S2_p1 = 2 << hs_sp.bit_range(es_type(1.0, 2)).first;

    // Expected classification of states
    std::vector<std::set<sv_index_type>> ref_cl{
        // S_z=-3
        {S0_m1 + S1_m1 + S2_m1},
        // S_z=-2
        {S0_m1 + S1_m1 + S2_0, S0_m1 + S1_0 + S2_m1, S0_0 + S1_m1 + S2_m1},
        // S_z=-1
        {S0_m1 + S1_0 + S2_0,
         S0_0 + S1_m1 + S2_0,
         S0_0 + S1_0 + S2_m1,
         S0_m1 + S1_m1 + S2_p1,
         S0_m1 + S1_p1 + S2_m1,
         S0_p1 + S1_m1 + S2_m1},
        // S_z=0
        {S0_0 + S1_0 + S2_0,
         S0_m1 + S1_p1 + S2_0,
         S0_p1 + S1_m1 + S2_0,
         S0_0 + S1_m1 + S2_p1,
         S0_0 + S1_p1 + S2_m1,
         S0_m1 + S1_0 + S2_p1,
         S0_p1 + S1_0 + S2_m1},
        // S_z=1
        {S0_p1 + S1_0 + S2_0,
         S0_0 + S1_p1 + S2_0,
         S0_0 + S1_0 + S2_p1,
         S0_p1 + S1_p1 + S2_m1,
         S0_p1 + S1_m1 + S2_p1,
         S0_m1 + S1_p1 + S2_p1},
        // S_z=2
        {S0_p1 + S1_p1 + S2_0, S0_p1 + S1_0 + S2_p1, S0_0 + S1_p1 + S2_p1},
        // S_z=3
        {S0_p1 + S1_p1 + S2_p1}};
    CHECK(collect_subspaces(sp) == subspaces_t(ref_cl.begin(), ref_cl.end()));

    auto op_pp = make_loperator(S_p<3>(0) + S_m<3>(1), hs_sp);
    auto op_mm = make_loperator(S_m<3>(0) + S_p<3>(1), hs_sp);
    sp.merge_subspaces(op_pp, op_mm);

    // Two subspaces: Even and odd S_z
    std::set<sv_index_type> ref_cl_odd, ref_cl_even;
    for(int S_z : {-3, -1, 1, 3})
      ref_cl_odd.insert(ref_cl[S_z + 3].begin(), ref_cl[S_z + 3].end());
    for(int S_z : {-2, 0, 2})
      ref_cl_even.insert(ref_cl[S_z + 3].begin(), ref_cl[S_z + 3].end());

    CHECK(sp.n_subspaces() == 2);
    CHECK(collect_subspaces(sp) == subspaces_t{ref_cl_odd, ref_cl_even});
  }
}
