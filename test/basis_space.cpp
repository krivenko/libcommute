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

#include <vector>

using namespace libcommute;

using bs_type = basis_space<std::string, int>;

template<typename V> void check_equality(V const& v) {
  for(size_t i1 = 0; i1 < v.size(); ++i1) {
    for(size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] == *v[i2]) == (i1 == i2));
      CHECK((*v[i1] != *v[i2]) == (i1 != i2));
    }
  }
}

template<typename V> void check_less_greater(V const& v) {
  for(size_t i1 = 0; i1 < v.size(); ++i1) {
    for(size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] < *v[i2]) == (i1 < i2));
      CHECK((*v[i1] > *v[i2]) == (i1 > i2));
    }
  }
}

TEST_CASE("Basis Hilbert space", "[basis_space]") {
  using namespace static_indices;

  // Setup

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

  SECTION("fermion") {
    for(auto * bs : fermion_bs) {
      CHECK(bs->algebra_id() == FERMION_ALGEBRA_ID);
      CHECK(bs->n_bits() == 1);
    }
    check_equality(fermion_bs);
    check_less_greater(fermion_bs);
  }

  SECTION("boson") {
    for(auto * bs : boson_bs) {
      CHECK(bs->algebra_id() == BOSON_ALGEBRA_ID);
      CHECK(bs->n_bits() == 4);
    }
    check_equality(boson_bs);
    check_less_greater(boson_bs);
  }

  SECTION("spin-1/2") {
    for(auto * bs : spin_bs) {
      CHECK(bs->algebra_id() == SPIN_ALGEBRA_ID);
      CHECK(bs->n_bits() == 1);
    }
    check_equality(spin_bs);
    check_less_greater(spin_bs);
  }

  SECTION("spin-1") {
    for(auto * bs : spin1_bs) {
      CHECK(bs->algebra_id() == SPIN_ALGEBRA_ID);
      CHECK(bs->n_bits() == 2);
    }
    check_equality(spin1_bs);
    check_less_greater(spin1_bs);
  }

  SECTION("spin-3/2") {
    for(auto * bs : spin32_bs) {
      CHECK(bs->algebra_id() == SPIN_ALGEBRA_ID);
      CHECK(bs->n_bits() == 2);
    }
    check_equality(spin32_bs);
    check_less_greater(spin32_bs);
  }

  SECTION("different base spaces") {
    std::vector<bs_type*> all_bs;
    for(auto const & bs_list : {fermion_bs,
                                boson_bs,
                                spin_bs,
                                spin1_bs,
                                spin32_bs}) {
      for(auto const& bs : bs_list) all_bs.push_back(bs);
    }
    check_equality(all_bs);
    check_less_greater(all_bs);
  }
}
