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
#include <libcommute/qoperator/qoperator.hpp>
#include <libcommute/qoperator/remapped_basis_view.hpp>

#include <map>
#include <vector>
#include <unordered_map>

using namespace libcommute;

TEST_CASE("Basis-remapped view of a state vector",
          "[remapped_basis_view]") {

  using namespace static_indices::real;
  using state_vector = std::vector<double>;

  // Finite system of 4 fermions: 2 orbitals, two spin projections
  // Hamiltonian: spin flips
  auto Hex = 2 * c_dag("up", 1) * c("up", 2) * c_dag("dn", 2) * c("dn", 1) +
             2 * c_dag("up", 2) * c("up", 1) * c_dag("dn", 1) * c("dn", 2);
  auto Hp = 2 * c_dag("up", 1) * c("up", 2) * c_dag("dn", 1) * c("dn", 2) +
            2 * c_dag("up", 2) * c("up", 1) * c_dag("dn", 2) * c("dn", 1);

  auto hs = make_hilbert_space(Hex + Hp);
  REQUIRE(hs.dim() == 16);

  using static_indices::make_space_fermion;
  REQUIRE(hs.bit_range(make_space_fermion("dn", 1)) == std::make_pair(0, 0));
  REQUIRE(hs.bit_range(make_space_fermion("dn", 2)) == std::make_pair(1, 1));
  REQUIRE(hs.bit_range(make_space_fermion("up", 1)) == std::make_pair(2, 2));
  REQUIRE(hs.bit_range(make_space_fermion("up", 2)) == std::make_pair(3, 3));

  // Map all basis states with 2 electrons so that their indices are contiguous
  std::unordered_map<sv_index_type, sv_index_type> map;
  // 3 = 1 + 2 -> |dn>_1 |dn>_2
  // 5 = 1 + 4 -> |dn>_1 |up>_1
  // 6 = 2 + 4 -> |dn>_2 |up>_1
  // 9 = 1 + 8 -> |dn>_1 |up>_2
  // 10 = 2 + 8 -> |dn>_2 |up>_2
  // 12 = 4 + 8 -> |up>_1 |up>_2
  for(auto i : {3, 5, 6, 9, 10, 12}) map[i] = map.size();
  REQUIRE(map.size() == 6);

  state_vector st{0, 1, 2, 3, 4, 5};

  using foreach_res_t = std::map<sv_index_type, double>;
  const foreach_res_t foreach_res_ref = {
    std::make_pair(5, 2.0),
    std::make_pair(6, 4.0),
    std::make_pair(9, 6.0),
    std::make_pair(10, 8.0),
    std::make_pair(12, 10.0)
  };

  SECTION("remapped_basis_view") {

    SECTION("const") {
      auto view = remapped_basis_view<state_vector, true>(st, map);

      CHECK(std::is_same<element_type_t<decltype(view)>, double>::value);
      CHECK(get_element(view, 10) == 4);

      foreach_res_t foreach_res;
      auto f = [&foreach_res](sv_index_type n, double a) {
        foreach_res[n] = 2*a;
      };
      foreach(view, f);
      CHECK(foreach_res == foreach_res_ref);
    }

    SECTION("non-const") {
      auto view = remapped_basis_view<state_vector>(st, map);

      CHECK(std::is_same<element_type_t<decltype(view)>, double>::value);
      CHECK(get_element(view, 10) == 4);

      foreach_res_t foreach_res;
      auto f = [&foreach_res](sv_index_type n, double a) {
        foreach_res[n] = 2*a;
      };
      foreach(view, f);
      CHECK(foreach_res == foreach_res_ref);

      update_add_element(view, 6, -3);
      CHECK(st == state_vector{0, 1, -1, 3, 4, 5});
      set_zeros(view);
      CHECK(st == state_vector{0, 0, 0, 0, 0, 0});
    }
  }

  SECTION("qoperator") {
    state_vector out(6);

    SECTION("Spin flips") {
      auto Hop = make_qoperator(Hex, hs);

      state_vector in1{1, 1, 1, 1, 1, 1};
      Hop(remapped_basis_view<state_vector, true>(in1, map),
          remapped_basis_view<state_vector, false>(out, map)
      );
      CHECK(out == state_vector{0, 0, 2, 2, 0, 0});

      state_vector in2{1, 1, 1, -1, 1, 1};
      Hop(remapped_basis_view<state_vector, true>(in2, map),
          remapped_basis_view<state_vector, false>(out, map)
      );
      CHECK(out == state_vector{0, 0, -2, 2, 0, 0});
    }

    SECTION("Pair hops") {
      auto Hop = make_qoperator(Hp, hs);

      state_vector in1{1, 1, 1, 1, 1, 1};
      Hop(remapped_basis_view<state_vector, true>(in1, map),
          remapped_basis_view<state_vector, false>(out, map)
      );
      CHECK(out == state_vector{0, 2, 0, 0, 2, 0});

      state_vector in2{1, 1, 1, 1, -1, 1};
      Hop(remapped_basis_view<state_vector, true>(in2, map),
          remapped_basis_view<state_vector, false>(out, map)
      );
      CHECK(out == state_vector{0, -2, 0, 0, 2, 0});
    }
  }
}
