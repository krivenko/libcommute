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
#include <libcommute/qoperator/mapped_basis_view.hpp>

#include <map>
#include <set>
#include <vector>
#include <unordered_map>

using namespace libcommute;

template<typename Key, typename Value>
void check_equal_maps_up_to_value_permutation(
  std::unordered_map<Key, Value> const& m1,
  std::unordered_map<Key, Value> const& m2
)
{
  std::set<Key> keys1, keys2;
  std::set<Value> values1, values2;
  for(auto const& p : m1) {
    keys1.emplace(p.first);
    values1.emplace(p.second);
  }
  for(auto const& p : m2) {
    keys2.emplace(p.first);
    values2.emplace(p.second);
  }
  CHECK(keys1 == keys2);
  CHECK(values2 == values2);
}

TEST_CASE("Basis-mapped view of a state vector",
          "[mapped_basis_view]") {

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

  // 3 = 1 + 2 -> |dn>_1 |dn>_2
  // 5 = 1 + 4 -> |dn>_1 |up>_1
  // 6 = 2 + 4 -> |dn>_2 |up>_1
  // 9 = 1 + 8 -> |dn>_1 |up>_2
  // 10 = 2 + 8 -> |dn>_2 |up>_2
  // 12 = 4 + 8 -> |up>_1 |up>_2

  // Map all basis states with 2 electrons so that their indices are contiguous
  std::unordered_map<sv_index_type, sv_index_type> map;
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

  SECTION("mapped_basis_view") {

    SECTION("const") {
      auto view = mapped_basis_view<state_vector, true>(st, map);

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
      auto view = mapped_basis_view<state_vector>(st, map);

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
      Hop(mapped_basis_view<state_vector, true>(in1, map),
          mapped_basis_view<state_vector, false>(out, map)
      );
      CHECK(out == state_vector{0, 0, 2, 2, 0, 0});

      state_vector in2{1, 1, 1, -1, 1, 1};
      Hop(mapped_basis_view<state_vector, true>(in2, map),
          mapped_basis_view<state_vector, false>(out, map)
      );
      CHECK(out == state_vector{0, 0, -2, 2, 0, 0});
    }

    SECTION("Pair hops") {
      auto Hop = make_qoperator(Hp, hs);

      state_vector in1{1, 1, 1, 1, 1, 1};
      Hop(mapped_basis_view<state_vector, true>(in1, map),
          mapped_basis_view<state_vector, false>(out, map)
      );
      CHECK(out == state_vector{0, 2, 0, 0, 2, 0});

      state_vector in2{1, 1, 1, 1, -1, 1};
      Hop(mapped_basis_view<state_vector, true>(in2, map),
          mapped_basis_view<state_vector, false>(out, map)
      );
      CHECK(out == state_vector{0, -2, 0, 0, 2, 0});
    }
  }

  SECTION("basis_mapper") {
    state_vector st{1, 1, 1, 1, 1, 1};

    SECTION("basis_state_indices") {
      std::vector<sv_index_type> basis_indices{3, 5, 6, 9, 10, 12};
      basis_mapper mapper(basis_indices);
      CHECK(mapper.size() == 6);
      SECTION("const") {
        auto view = mapper.make_const_view(st);
        CHECK(view.map == map);
      }
      SECTION("non-const") {
        auto view = mapper.make_view(st);
        CHECK(view.map == map);
      }
    }

    SECTION("O|vac>") {
      auto P = c_dag("dn", 1) * c_dag("dn", 2) +
               c_dag("dn", 1) * c_dag("up", 1) +
               c_dag("dn", 2) * c_dag("up", 1) +
               c_dag("dn", 1) * c_dag("up", 2) +
               c_dag("dn", 2) * c_dag("up", 2) +
               c_dag("up", 1) * c_dag("up", 2);
      basis_mapper mapper(make_qoperator(P, hs), hs);
      CHECK(mapper.size() == 6);
      SECTION("const") {
        auto view = mapper.make_const_view(st);
        check_equal_maps_up_to_value_permutation(view.map, map);
      }
      SECTION("non-const") {
        auto view = mapper.make_view(st);
        check_equal_maps_up_to_value_permutation(view.map, map);
      }
    }

    SECTION("Compositions") {
      using O_list_t = std::vector<qoperator<double, fermion, boson, spin>>;

      basis_mapper mapper_empty(O_list_t{}, hs, 0);
      CHECK(mapper_empty.size() == 1);

      O_list_t O_list{
        make_qoperator(c_dag("dn", 1), hs),
        make_qoperator(c_dag("dn", 2), hs),
        make_qoperator(c_dag("up", 1), hs),
        make_qoperator(c_dag("up", 2), hs)
      };

      basis_mapper mapper_N0(O_list, hs, 0);
      CHECK(mapper_N0.size() == 1);
      basis_mapper mapper(O_list, hs, 2);

      SECTION("const") {
        CHECK(mapper_empty.make_const_view(st).map.size() == 1);
        CHECK(mapper_empty.make_const_view(st).map.at(0) == 0);
        CHECK(mapper_N0.make_const_view(st).map.size() == 1);
        CHECK(mapper_N0.make_const_view(st).map.at(0) == 0);

        auto view = mapper.make_const_view(st);
        check_equal_maps_up_to_value_permutation(view.map, map);
      }
      SECTION("non-const") {
        CHECK(mapper_empty.make_view(st).map.size() == 1);
        CHECK(mapper_empty.make_view(st).map.at(0) == 0);
        CHECK(mapper_N0.make_view(st).map.size() == 1);
        CHECK(mapper_N0.make_view(st).map.at(0) == 0);

        auto view = mapper.make_view(st);
        check_equal_maps_up_to_value_permutation(view.map, map);
      }
    }

    SECTION("Compositions/bosons") {
      using O_list_t = std::vector<qoperator<double, fermion, boson, spin>>;

      auto hs = make_hilbert_space(a_dag(1) + a_dag(2) + a_dag(3) + a_dag(4),
                                   boson_bs_constructor(4));

      O_list_t O_list{
        make_qoperator(a_dag(1), hs),
        make_qoperator(a_dag(2), hs),
        make_qoperator(a_dag(3), hs),
        make_qoperator(a_dag(4), hs)
      };

      std::vector<int> map_size_ref{1, 4, 10, 20, 35, 56, 84, 120, 165, 220};
      for(int N = 0; N < 10; ++N) {
        basis_mapper mapper(O_list, hs, N);
        CHECK(mapper.size() == map_size_ref[N]);
      }
    }
  }
}
