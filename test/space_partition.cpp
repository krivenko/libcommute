/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include "catch2/catch.hpp"

#include <libcommute/expression/factories.hpp>
#include <libcommute/loperator/space_partition.hpp>

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

using namespace libcommute;

TEST_CASE("Automatic Hilbert space partition", "[space_partition]") {
  using namespace static_indices::real;

  // 3 orbital Hubbard-Kanamori atom
  const int n_orbs = 3;
  const double mu = 0.7;
  const double U = 3.0;
  const double J = 0.3;

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
  for (int o1 = 0; o1 < n_orbs; ++o1) {
    for (int o2 = 0; o2 < n_orbs; ++o2) {
      if (o1 == o2) continue;
      H += (U - 2 * J) * n("up", o1) * n("dn", o2);
    }
  }
  // Interorbital interactions, equal spins
  for (int o1 = 0; o1 < n_orbs; ++o1) {
    for (int o2 = 0; o2 < n_orbs; ++o2) {
      if (o2 >= o1) continue;
      H += (U - 3 * J) * n("up", o1) * n("up", o2);
      H += (U - 3 * J) * n("dn", o1) * n("dn", o2);
    }
  }
  // Spin flips and pair hoppings
  for (int o1 = 0; o1 < n_orbs; ++o1) {
    for (int o2 = 0; o2 < n_orbs; ++o2) {
      if (o1 == o2) continue;
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
    auto sp = space_partition(Hop, hs);

    CHECK(sp.dim() == 64);
    CHECK(sp.n_subspaces() == 44);

    // Calculated classification of states
    // sets are used to neglect order of subspaces and of states within
    // a subspace
    std::vector<std::set<sv_index_type>> v_cl(sp.n_subspaces());
    foreach(sp, [&](int i, int subspace) { v_cl[subspace].insert(i); });
    std::set<std::set<sv_index_type>> cl{v_cl.cbegin(), v_cl.cend()};

    // Expected classification of states
    std::set<std::set<sv_index_type>> ref_cl {
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
      {d0 + d1 + d2 + u0 + u1 + u2}
    };

    CHECK(cl == ref_cl);
  }

  SECTION("Matrix elements") {
    matrix_elements_map<double> matrix_elements;
    auto sp = space_partition(Hop, hs, matrix_elements);

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
      v_melem[sp[x.first.first]].insert(melem_t{x.first.second,
                                                x.first.first,
                                                x.second}
                                       );
    }
    melem_set_t melem(v_melem.cbegin(), v_melem.cend());

    // Expected matrix elements
    melem_set_t ref_melem {
      // N=0
      std::set<melem_t>{},
      // N=1
      {{d0, d0, -mu}}, {{d1, d1, -mu}}, {{d2, d2, -mu}},
      {{u0, u0, -mu}}, {{u1, u1, -mu}}, {{u2, u2, -mu}},
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
       {d0 + u0, d1 + u1, J}, {d0 + u0, d2 + u2, J}, {d1 + u1, d2 + u2, J},
       {d1 + u1, d0 + u0, J}, {d2 + u2, d0 + u0, J}, {d2 + u2, d1 + u1, J}},
      // N=2, spin flip
      {{d0 + u1, d0 + u1, -2 * mu + U - 2 * J},
       {d1 + u0, d1 + u0, -2 * mu + U - 2 * J},
       {d0 + u1, d1 + u0, J}, {d1 + u0, d0 + u1, J}},
      {{d0 + u2, d0 + u2, -2 * mu + U - 2 * J},
       {d2 + u0, d2 + u0, -2 * mu + U - 2 * J},
       {d0 + u2, d2 + u0, J}, {d2 + u0, d0 + u2, J}},
      {{d1 + u2, d1 + u2, -2 * mu + U - 2 * J},
       {d2 + u1, d2 + u1, -2 * mu + U - 2 * J},
       {d1 + u2, d2 + u1, J}, {d2 + u1, d1 + u2, J}},
      // N=3
      {{d0 + d1 + d2, d0 + d1 + d2, -3 * mu + 3 * U - 9 * J}},
      {{u0 + u1 + u2, u0 + u1 + u2, -3 * mu + 3 * U - 9 * J}},
      {{d0 + d1 + u0, d0 + d1 + u0, -3 * mu + 3 * U - 5 * J},
       {d1 + d2 + u2, d1 + d2 + u2, -3 * mu + 3 * U - 5 * J},
       {d0 + d1 + u0, d1 + d2 + u2, -J}, {d1 + d2 + u2, d0 + d1 + u0, -J}},
      {{d0 + d2 + u0, d0 + d2 + u0, -3 * mu + 3 * U - 5 * J},
       {d1 + d2 + u1, d1 + d2 + u1, -3 * mu + 3 * U - 5 * J},
       {d0 + d2 + u0, d1 + d2 + u1, J}, {d1 + d2 + u1, d0 + d2 + u0, J}},
      {{d0 + d1 + u1, d0 + d1 + u1, -3 * mu + 3 * U - 5 * J},
       {d0 + d2 + u2, d0 + d2 + u2, -3 * mu + 3 * U - 5 * J},
       {d0 + d1 + u1, d0 + d2 + u2, J}, {d0 + d2 + u2, d0 + d1 + u1, J}},
      {{d0 + u0 + u1, d0 + u0 + u1, -3 * mu + 3 * U - 5 * J},
       {d2 + u1 + u2, d2 + u1 + u2, -3 * mu + 3 * U - 5 * J},
       {d0 + u0 + u1, d2 + u1 + u2, -J}, {d2 + u1 + u2, d0 + u0 + u1, -J}},
      {{d1 + u0 + u1, d1 + u0 + u1, -3 * mu + 3 * U - 5 * J},
       {d2 + u0 + u2, d2 + u0 + u2, -3 * mu + 3 * U - 5 * J},
       {d1 + u0 + u1, d2 + u0 + u2, J}, {d2 + u0 + u2, d1 + u0 + u1, J}},
      {{d0 + u0 + u2, d0 + u0 + u2, -3 * mu + 3 * U - 5 * J},
       {d1 + u1 + u2, d1 + u1 + u2, -3 * mu + 3 * U - 5 * J},
       {d0 + u0 + u2, d1 + u1 + u2, J}, {d1 + u1 + u2, d0 + u0 + u2, J}},
      {{d1 + d2 + u0, d1 + d2 + u0, -3 * mu + 3 * U - 7 * J},
       {d0 + d2 + u1, d0 + d2 + u1, -3 * mu + 3 * U - 7 * J},
       {d0 + d1 + u2, d0 + d1 + u2, -3 * mu + 3 * U - 7 * J},
       {d1 + d2 + u0, d0 + d2 + u1, J},
       {d0 + d2 + u1, d1 + d2 + u0, J}, {d1 + d2 + u0, d0 + d1 + u2, -J},
       {d0 + d1 + u2, d1 + d2 + u0, -J}, {d0 + d2 + u1, d0 + d1 + u2, J},
       {d0 + d1 + u2, d0 + d2 + u1, J}},
      {{d2 + u0 + u1, d2 + u0 + u1, -3 * mu + 3 * U - 7 * J},
       {d0 + u1 + u2, d0 + u1 + u2, -3 * mu + 3 * U - 7 * J},
       {d1 + u0 + u2, d1 + u0 + u2, -3 * mu + 3 * U - 7 * J},
       {d2 + u0 + u1, d0 + u1 + u2, -J},
       {d0 + u1 + u2, d2 + u0 + u1, -J}, {d2 + u0 + u1, d1 + u0 + u2, J},
       {d1 + u0 + u2, d2 + u0 + u1, J}, {d0 + u1 + u2, d1 + u0 + u2, J},
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
        -6 * mu + 15 * U - 30 * J}}
    };

    CHECK(melem == ref_melem);
  }

  SECTION("merge_subspaces()") {

    auto sp = space_partition(Hop, hs);

    std::vector<decltype(Hop)> Cd, C, all_ops;
    for(std::string spin : {"dn", "up"}) {
      for(int o = 0; o < n_orbs; ++o) {
        Cd.emplace_back(c_dag(spin, o), hs);
        C.emplace_back(c(spin, o), hs);

        all_ops.emplace_back(Cd.back());
        all_ops.emplace_back(C.back());

        sp.merge_subspaces(Cd.back(), C.back(), hs);
      }
    }

    // Calculated classification of states
    std::vector<std::set<sv_index_type>> v_cl(sp.n_subspaces());
    foreach(sp, [&](int i, int subspace) { v_cl[subspace].insert(i); });
    std::set<std::set<sv_index_type>> cl{v_cl.cbegin(), v_cl.cend()};

    std::vector<double> in_state(sp.dim());

    for(auto const& op : all_ops) {
      for(auto const& i_sp : cl) {
        std::set<sv_index_type> f_sp;
        for(auto i : i_sp) {
          in_state[i] = 1.0;
          auto out_state = op(in_state);
          foreach(out_state, [&f_sp](sv_index_type f, double a) {
            if(std::abs(a) < 1e-10) return;
            f_sp.insert(f);
          });
          in_state[i] = 0;
        }

        // op maps i_sp onto zero
        if(f_sp.size() == 0) continue;

        // Check if op maps i_sp to only one subspace
        int n = 0;
        for(auto const& f_sp_ref : cl) {
          if(std::includes(f_sp_ref.cbegin(), f_sp_ref.cend(),
                           f_sp.cbegin(), f_sp.cend())) ++n;
        }
        CHECK(n == 1);
      }
    }
  }
}
