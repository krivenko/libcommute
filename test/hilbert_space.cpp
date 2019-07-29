/*******************************************************************************
 *
 * This file is part of libcommute, a C++11/14/17 header-only library allowing
 * to manipulate polynomial expressions with quantum-mechanical operators.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include "catch2/catch.hpp"

#include <libcommute/qoperator/basis_space_fermion.hpp>
#include <libcommute/qoperator/basis_space_boson.hpp>
#include <libcommute/qoperator/basis_space_spin.hpp>
#include <libcommute/qoperator/hilbert_space.hpp>

#include <utility>

using namespace libcommute;

TEST_CASE("Hilbert space", "[hilbert_space]") {
  using namespace static_indices;

  // Setup
  using bs_type = basis_space<std::string, int>;
  using hs_type = hilbert_space<std::string, int>;

  // Fermionic generators
  auto bs_f_dn = make_space_fermion("dn", 0);
  auto bs_f_up = make_space_fermion("up", 0);
  std::vector<bs_type*> fermion_bs = {&bs_f_dn,&bs_f_up};

  // Bosonic generators
  auto bs_b_x = make_space_boson(4, "x", 0);
  auto bs_b_y = make_space_boson(4, "y", 0);
  std::vector<bs_type*> boson_bs = {&bs_b_x,&bs_b_y};

  // Spin-1/2 algebra generators
  auto bs_s_i = make_space_spin(0.5, "i", 0);
  auto bs_s_j = make_space_spin(0.5, "j", 0);
  std::vector<bs_type*> spin_bs = {&bs_s_i,&bs_s_j};

  // Spin-1 algebra generators
  auto bs_s1_i = make_space_spin(1.0, "i", 0);
  auto bs_s1_j = make_space_spin(1.0, "j", 0);
  std::vector<bs_type*> spin1_bs = {&bs_s1_i,&bs_s1_j};

  // Spin-3/2 algebra generators
  auto bs_s32_i = make_space_spin(3.0/2, "i", 0);
  auto bs_s32_j = make_space_spin(3.0/2, "j", 0);
  std::vector<bs_type*> spin32_bs = {&bs_s32_i,&bs_s32_j};

  SECTION("Equality") {
    hs_type hs_empty;
    CHECK(hs_empty == hs_empty);
    CHECK_FALSE(hs_empty != hs_empty);

    hs_type hs(bs_s32_i, bs_s32_j,
               bs_s1_i, bs_s1_j,
               bs_s_i, bs_s_j,
               bs_b_x, bs_b_y,
               bs_f_dn, bs_f_up
              );
    CHECK(hs == hs);
    CHECK_FALSE(hs != hs);
  }

  SECTION("Construction") {
    hs_type hs_empty;
    CHECK(hs_empty.size() == 0);
    CHECK(hs_empty.total_n_bits() == 0);

    hs_type hs1(bs_s32_i, bs_s32_j,
                bs_s1_i, bs_s1_j,
                bs_s_i, bs_s_j,
                bs_b_x, bs_b_y,
                bs_f_dn, bs_f_up
               );
    CHECK(hs1.size() == 10);
    CHECK(hs1.total_n_bits() == 20);

    hs_type hs2(hs1);
    CHECK(hs2 == hs1);

    hs_type hs3(std::move(hs2));
    CHECK(hs3 == hs1);

    CHECK_THROWS_AS(hs_type(bs_s32_i, bs_s_i, bs_s_i, bs_f_dn),
                    hs_type::basis_space_exists);
  }

  SECTION("Assignment") {
    hs_type hs_empty;
    hs_type hs1(bs_s32_i, bs_s32_j, bs_s1_i, bs_s1_j);
    hs_type hs2(bs_s_i, bs_s_j, bs_b_x, bs_b_y, bs_f_dn, bs_f_up);
    hs_type hs3(bs_s_i, bs_b_x, bs_b_y, bs_f_up);
    CHECK_FALSE(hs_empty == hs2);
    CHECK_FALSE(hs1 == hs2);
    SECTION("Copy") {
      hs_empty = hs2;
      CHECK(hs_empty == hs2);
      hs1 = hs2;
      CHECK(hs1 == hs2);
    }
    SECTION("Move") {
      auto hs_ref = hs2;
      hs_empty = std::move(hs2);
      CHECK(hs_empty == hs_ref);
      hs_ref = hs3;
      hs1 = std::move(hs3);
      CHECK(hs1 == hs_ref);
    }
  }

  SECTION("has() and bit_range()") {
    hs_type hs(bs_s32_i, bs_s32_j,
               bs_s1_j,
               bs_s_i, bs_s_j,
               bs_b_x,
               bs_f_dn, bs_f_up
              );
    CHECK(hs.size() == 8);
    CHECK(hs.total_n_bits() == 14);

    CHECK(hs.has(bs_s32_i));
    CHECK(hs.bit_range(bs_s32_i) == std::make_pair(0, 1));
    CHECK(hs.has(bs_s32_j));
    CHECK(hs.bit_range(bs_s32_j) == std::make_pair(2, 3));
    CHECK_FALSE(hs.has(bs_s1_i));
    CHECK_THROWS_AS(hs.bit_range(bs_s1_i), hs_type::basis_space_not_found);
    CHECK(hs.has(bs_s1_j));
    CHECK(hs.bit_range(bs_s1_j) == std::make_pair(4, 5));
    CHECK(hs.has(bs_s_i));
    CHECK(hs.bit_range(bs_s_i) == std::make_pair(6, 6));
    CHECK(hs.has(bs_s_j));
    CHECK(hs.bit_range(bs_s_j) == std::make_pair(7, 7));
    CHECK(hs.has(bs_b_x));
    CHECK(hs.bit_range(bs_b_x) == std::make_pair(8, 11));
    CHECK_FALSE(hs.has(bs_b_y));
    CHECK_THROWS_AS(hs.bit_range(bs_b_y), hs_type::basis_space_not_found);
    CHECK(hs.has(bs_f_dn));
    CHECK(hs.bit_range(bs_f_dn) == std::make_pair(12, 12));
    CHECK(hs.has(bs_f_up));
    CHECK(hs.bit_range(bs_f_up) == std::make_pair(13, 13));
  }

  SECTION("add()") {
    hs_type hs(bs_s32_i);

    auto check_hs = [&hs](bs_type const& bs,
                          int size,
                          int total_n_bits,
                          int b, int e) {
      CHECK(hs.size() == size);
      CHECK(hs.total_n_bits() == total_n_bits);
      CHECK(hs.has(bs));
      CHECK(hs.bit_range(bs) == std::make_pair(b, e));
    };

    check_hs(bs_s32_i, 1, 2, 0, 1);
    hs.add(bs_s32_j);
    check_hs(bs_s32_j, 2, 4, 2, 3);
    hs.add(bs_s1_j);
    check_hs(bs_s1_j, 3, 6, 4, 5);
    hs.add(bs_s_i);
    check_hs(bs_s_i, 4, 7, 6, 6);
    hs.add(bs_s_j);
    check_hs(bs_s_j, 5, 8, 7, 7);
    hs.add(bs_b_x);
    check_hs(bs_b_x, 6, 12, 8, 11);
    hs.add(bs_f_dn);
    check_hs(bs_f_dn, 7, 13, 12, 12);
    hs.add(bs_f_up);
    check_hs(bs_f_up, 8, 14, 13, 13);

    CHECK_THROWS_AS(hs.add(bs_s_j), hs_type::basis_space_exists);
  }

  SECTION("Construction from expression") {
    // TODO
  }

}
