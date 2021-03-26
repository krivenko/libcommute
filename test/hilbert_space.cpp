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

#include "catch2/catch.hpp"

#include <libcommute/expression/factories.hpp>
#include <libcommute/loperator/elementary_space_fermion.hpp>
#include <libcommute/loperator/elementary_space_boson.hpp>
#include <libcommute/loperator/elementary_space_spin.hpp>
#include <libcommute/loperator/hilbert_space.hpp>

#include <utility>

using namespace libcommute;

TEST_CASE("Hilbert space", "[hilbert_space]") {
  using namespace static_indices;

  // Setup
  using es_type = elementary_space<std::string, int>;
  using hs_type = hilbert_space<std::string, int>;

  // Fermionic elementary spaces
  auto es_f_dn = make_space_fermion("dn", 0);
  auto es_f_up = make_space_fermion("up", 0);
  std::vector<es_type*> fermion_es = {&es_f_dn,&es_f_up};

  // Bosonic elementary spaces
  auto es_b_x = make_space_boson(4, "x", 0);
  auto es_b_y = make_space_boson(4, "y", 0);
  std::vector<es_type*> boson_es = {&es_b_x,&es_b_y};

  // Spin-1/2 elementary spaces
  auto es_s_i = make_space_spin(0.5, "i", 0);
  auto es_s_j = make_space_spin(0.5, "j", 0);
  std::vector<es_type*> spin_es = {&es_s_i,&es_s_j};

  // Spin-1 elementary spaces
  auto es_s1_i = make_space_spin(1.0, "i", 0);
  auto es_s1_j = make_space_spin(1.0, "j", 0);
  std::vector<es_type*> spin1_es = {&es_s1_i,&es_s1_j};

  // Spin-3/2 elementary spaces
  auto es_s32_i = make_space_spin(3.0/2, "i", 0);
  auto es_s32_j = make_space_spin(3.0/2, "j", 0);
  std::vector<es_type*> spin32_es = {&es_s32_i,&es_s32_j};

  SECTION("Equality") {
    hs_type hs_empty;
    CHECK(hs_empty == hs_empty);
    CHECK_FALSE(hs_empty != hs_empty);

    hs_type hs(es_s32_i, es_s32_j,
               es_s1_i, es_s1_j,
               es_s_i, es_s_j,
               es_b_x, es_b_y,
               es_f_dn, es_f_up
              );
    CHECK(hs == hs);
    CHECK_FALSE(hs != hs);
  }

  SECTION("Construction") {
    hs_type hs_empty;
    CHECK(hs_empty.size() == 0);
    CHECK(hs_empty.total_n_bits() == 0);
    CHECK(hs_empty.dim() == 1);

    hs_type hs1(es_s32_i, es_s32_j,
                es_s1_i, es_s1_j,
                es_s_i, es_s_j,
                es_b_x, es_b_y,
                es_f_dn, es_f_up
               );
    CHECK(hs1.size() == 10);
    CHECK(hs1.total_n_bits() == 20);
    CHECK(hs1.dim() == 1048576);

    hs_type hs2(hs1);
    CHECK(hs2 == hs1);

    hs_type hs3(std::move(hs2));
    CHECK(hs3 == hs1);

    CHECK_THROWS_AS(hs_type(es_s32_i, es_s_i, es_s_i, es_f_dn),
                    hs_type::elementary_space_exists);

    CHECK(hs_type(fermion_es) == hs_type(es_f_dn, es_f_up));
    CHECK(hs_type(boson_es) == hs_type(es_b_x, es_b_y));
    CHECK(hs_type(spin_es) == hs_type(es_s_i, es_s_j));
    CHECK(hs_type(spin1_es) == hs_type(es_s1_i, es_s1_j));
    CHECK(hs_type(spin32_es) == hs_type(es_s32_i, es_s32_j));
  }

  SECTION("Assignment") {
    hs_type hs_empty;
    hs_type hs1(es_s32_i, es_s32_j, es_s1_i, es_s1_j);
    hs_type hs2(es_s_i, es_s_j, es_b_x, es_b_y, es_f_dn, es_f_up);
    hs_type hs3(es_s_i, es_b_x, es_b_y, es_f_up);
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

  SECTION("has(), index(), bit_range() and basis_state_index()") {
    hs_type hs(es_s32_i, es_s32_j,
               es_s1_j,
               es_s_i, es_s_j,
               es_b_x,
               es_f_dn, es_f_up
              );
    CHECK(hs.size() == 8);
    CHECK(hs.total_n_bits() == 14);
    CHECK(hs.dim() == 16384);
    CHECK(hs.algebra_bit_range(fermion) == std::make_pair(0, 1));
    CHECK(hs.algebra_bit_range(boson) == std::make_pair(2, 5));
    CHECK(hs.algebra_bit_range(spin) == std::make_pair(6, 13));

    CHECK(hs.has(es_f_dn));
    CHECK(hs.index(es_f_dn) == 0);
    CHECK(hs.bit_range(es_f_dn) == std::make_pair(0, 0));
    CHECK(hs.basis_state_index(es_f_dn, 0) == 0);
    CHECK(hs.basis_state_index(es_f_dn, 1) == 1);
    CHECK(hs.has(es_f_up));
    CHECK(hs.index(es_f_up) == 1);
    CHECK(hs.bit_range(es_f_up) == std::make_pair(1, 1));
    CHECK(hs.basis_state_index(es_f_up, 0) == 0);
    CHECK(hs.basis_state_index(es_f_up, 1) == 2);
    CHECK(hs.has(es_b_x));
    CHECK(hs.index(es_b_x) == 2);
    CHECK(hs.bit_range(es_b_x) == std::make_pair(2, 5));
    CHECK(hs.basis_state_index(es_b_x, 0) == 0);
    CHECK(hs.basis_state_index(es_b_x, 1) == 4);
    CHECK(hs.basis_state_index(es_b_x, 5) == 20);
    CHECK_FALSE(hs.has(es_b_y));
    CHECK_THROWS_AS(hs.index(es_b_y), hs_type::elementary_space_not_found);
    CHECK_THROWS_AS(hs.bit_range(es_b_y), hs_type::elementary_space_not_found);
    CHECK(hs.has(es_s_i));
    CHECK(hs.index(es_s_i) == 3);
    CHECK(hs.bit_range(es_s_i) == std::make_pair(6, 6));
    CHECK(hs.basis_state_index(es_s_i, 0) == 0);
    CHECK(hs.basis_state_index(es_s_i, 1) == 64);
    CHECK(hs.has(es_s_j));
    CHECK(hs.index(es_s_j) == 4);
    CHECK(hs.bit_range(es_s_j) == std::make_pair(7, 7));
    CHECK_FALSE(hs.has(es_s1_i));
    CHECK_THROWS_AS(hs.index(es_s1_i), hs_type::elementary_space_not_found);
    CHECK_THROWS_AS(hs.bit_range(es_s1_i), hs_type::elementary_space_not_found);
    CHECK(hs.has(es_s1_j));
    CHECK(hs.index(es_s1_j) == 5);
    CHECK(hs.bit_range(es_s1_j) == std::make_pair(8, 9));
    CHECK(hs.basis_state_index(es_s1_j, 0) == 0);
    CHECK(hs.basis_state_index(es_s1_j, 1) == 256);
    CHECK(hs.basis_state_index(es_s1_j, 2) == 512);
    CHECK(hs.has(es_s32_i));
    CHECK(hs.index(es_s32_i) == 6);
    CHECK(hs.bit_range(es_s32_i) == std::make_pair(10, 11));
    CHECK(hs.basis_state_index(es_s32_i, 0) == 0);
    CHECK(hs.basis_state_index(es_s32_i, 1) == 1024);
    CHECK(hs.basis_state_index(es_s32_i, 2) == 2048);
    CHECK(hs.basis_state_index(es_s32_i, 3) == 3072);
    CHECK(hs.has(es_s32_j));
    CHECK(hs.index(es_s32_j) == 7);
    CHECK(hs.bit_range(es_s32_j) == std::make_pair(12, 13));
    CHECK(hs.basis_state_index(es_s32_j, 0) == 0);
    CHECK(hs.basis_state_index(es_s32_j, 1) == 4096);
    CHECK(hs.basis_state_index(es_s32_j, 2) == 8192);
    CHECK(hs.basis_state_index(es_s32_j, 3) == 12288);
  }

  SECTION("add()") {
    hs_type hs(es_s32_i);

    auto check_hs = [&hs](es_type const& es,
                          int index,
                          int size,
                          int total_n_bits,
                          int b, int e,
                          int fermion_b, int fermion_e,
                          int boson_b, int boson_e,
                          int spin_b, int spin_e
                         ) {
      CHECK(hs.size() == size);
      CHECK(hs.total_n_bits() == total_n_bits);
      CHECK(hs.dim() == 1 << total_n_bits);
      CHECK(hs.has(es));
      CHECK(hs.index(es) == index);
      CHECK(hs.bit_range(es) == std::make_pair(b, e));

      if(fermion_b != -1) {
        auto ref_range = std::make_pair(fermion_b, fermion_e);
        CHECK(hs.algebra_bit_range(fermion) == ref_range);
      } else {
        CHECK_THROWS_AS(hs.algebra_bit_range(fermion),
                        std::runtime_error);
      }
      if(boson_b != -1) {
        auto ref_range = std::make_pair(boson_b, boson_e);
        CHECK(hs.algebra_bit_range(boson) == ref_range);
      } else {
        CHECK_THROWS_AS(hs.algebra_bit_range(boson),
                        std::runtime_error);
      }
      if(spin_b != -1) {
        auto ref_range = std::make_pair(spin_b, spin_e);
        CHECK(hs.algebra_bit_range(spin) == ref_range);
      } else {
        CHECK_THROWS_AS(hs.algebra_bit_range(spin),
                        std::runtime_error);
      }
    };

    check_hs(es_s32_i, 0, 1, 2, 0, 1, -1, -1, -1, -1, 0, 1);
    hs.add(es_s32_j);
    check_hs(es_s32_j, 1, 2, 4, 2, 3, -1, -1, -1, -1, 0, 3);
    hs.add(es_s1_j);
    check_hs(es_s1_j, 0, 3, 6, 0, 1, -1, -1, -1, -1, 0, 5);
    hs.add(es_s_i);
    check_hs(es_s_i, 0, 4, 7, 0, 0, -1, -1, -1, -1, 0, 6);
    hs.add(es_s_j);
    check_hs(es_s_j, 1, 5, 8, 1, 1, -1, -1, -1, -1, 0, 7);
    hs.add(es_b_x);
    check_hs(es_b_x, 0, 6, 12, 0, 3, -1, -1, 0, 3, 4, 11);
    hs.add(es_f_dn);
    check_hs(es_f_dn, 0, 7, 13, 0, 0, 0, 0, 1, 4, 5, 12);
    hs.add(es_f_up);
    check_hs(es_f_up, 1, 8, 14, 1, 1, 0, 1, 2, 5, 6, 13);

    CHECK_THROWS_AS(hs.add(es_s_j), hs_type::elementary_space_exists);

    check_hs(es_f_dn, 0, 8, 14, 0, 0, 0, 1, 2, 5, 6, 13);
    check_hs(es_f_up, 1, 8, 14, 1, 1, 0, 1, 2, 5, 6, 13);
    check_hs(es_b_x, 2, 8, 14, 2, 5, 0, 1, 2, 5, 6, 13);
    check_hs(es_s_i, 3, 8, 14, 6, 6, 0, 1, 2, 5, 6, 13);
    check_hs(es_s_j, 4, 8, 14, 7, 7, 0, 1, 2, 5, 6, 13);
    check_hs(es_s1_j, 5, 8, 14, 8, 9, 0, 1, 2, 5, 6, 13);
    check_hs(es_s32_i, 6, 8, 14, 10, 11, 0, 1, 2, 5, 6, 13);
    check_hs(es_s32_j, 7, 8, 14, 12, 13, 0, 1, 2, 5, 6, 13);
    CHECK(hs.algebra_bit_range(fermion) == std::make_pair(0, 1));
    CHECK(hs.algebra_bit_range(boson) == std::make_pair(2, 5));
    CHECK(hs.algebra_bit_range(spin) == std::make_pair(6, 13));
  }

  SECTION("Construction from expression") {
    auto expr = 2.0 * S_p<4>("i", 0) * S_m<4>("j", 0) +
                5.0 * n("up", 0) * n("dn", 0);

    hs_type hs1(expr);
    CHECK(hs1.size() == 4);
    CHECK(hs1.total_n_bits() == 6);
    CHECK(hs1.dim() == 64);
    CHECK(hs1.has(es_f_dn));
    CHECK(hs1.index(es_f_dn) == 0);
    CHECK(hs1.bit_range(es_f_dn) == std::make_pair(0, 0));
    CHECK(hs1.has(es_f_up));
    CHECK(hs1.has(es_f_up) == 1);
    CHECK(hs1.bit_range(es_f_up) == std::make_pair(1, 1));
    CHECK(hs1.has(es_s32_i));
    CHECK(hs1.index(es_s32_i) == 2);
    CHECK(hs1.bit_range(es_s32_i) == std::make_pair(2, 3));
    CHECK(hs1.has(es_s32_j));
    CHECK(hs1.index(es_s32_j) == 3);
    CHECK(hs1.bit_range(es_s32_j) == std::make_pair(4, 5));
    CHECK(make_hilbert_space(expr) == hs1);

    SECTION("foreach()") {
      size_t count = 0;
      foreach(hs1, [&count](int i) { count += i; });
      CHECK(count == 2016);
    }

    expr += a_dag("x", 0) + a("y", 0);

    using ex_type = es_construction_failure<std::string, int>;
    CHECK_THROWS_AS(hs_type(expr), ex_type);

    hs_type hs2(expr, boson_es_constructor(4));
    CHECK(hs2.size() == 6);
    CHECK(hs2.total_n_bits() == 14);
    CHECK(hs2.dim() == 16384);
    CHECK(hs2.has(es_f_dn));
    CHECK(hs2.index(es_f_dn) == 0);
    CHECK(hs2.bit_range(es_f_dn) == std::make_pair(0, 0));
    CHECK(hs2.has(es_f_up));
    CHECK(hs2.index(es_f_up) == 1);
    CHECK(hs2.bit_range(es_f_up) == std::make_pair(1, 1));
    CHECK(hs2.has(es_b_x));
    CHECK(hs2.index(es_b_x) == 2);
    CHECK(hs2.bit_range(es_b_x) == std::make_pair(2, 5));
    CHECK(hs2.has(es_b_y));
    CHECK(hs2.index(es_b_y) == 3);
    CHECK(hs2.bit_range(es_b_y) == std::make_pair(6, 9));
    CHECK(hs2.has(es_s32_i));
    CHECK(hs2.index(es_s32_i) == 4);
    CHECK(hs2.bit_range(es_s32_i) == std::make_pair(10, 11));
    CHECK(hs2.has(es_s32_j));
    CHECK(hs2.index(es_s32_j) == 5);
    CHECK(hs2.bit_range(es_s32_j) == std::make_pair(12, 13));
    CHECK(make_hilbert_space(expr, boson_es_constructor(4)) == hs2);
  }

  SECTION("Very big Hilbert space") {
    hs_type hs1;
    for(int i = 0; i < 32; ++i)
      hs1.add(make_space_spin(3.0/2, "s", i));

    using ex_type = hs_type::hilbert_space_too_big;
    CHECK_THROWS_AS(hs1.add(make_space_spin(3.0/2, "s", 32)), ex_type);

    auto expr = expression<double, std::string, int>(1);
    for(int i = 0; i < 32; ++i) expr *= S_p<4>("s", i);
    hs_type hs2(expr);
    CHECK(hs2.total_n_bits() == 64);

    expr *= S_p<4>("s", 32);
    CHECK_THROWS_AS(hs_type(expr), ex_type);
  }
}
