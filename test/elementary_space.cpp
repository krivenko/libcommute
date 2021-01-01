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

#include "check_ordering.hpp"

#include <libcommute/loperator/elementary_space_fermion.hpp>
#include <libcommute/loperator/elementary_space_boson.hpp>
#include <libcommute/loperator/elementary_space_spin.hpp>

#include <vector>

using namespace libcommute;

using es_type = elementary_space<std::string, int>;

TEST_CASE("Elementary Hilbert space", "[elementary_space]") {
  using namespace static_indices;

  // Setup

  // Fermionic elementary spaces
  auto es_f_dn = make_space_fermion("dn", 0);
  auto es_f_up = make_space_fermion("up", 0);
  std::vector<es_type*> fermion_es = {&es_f_dn,&es_f_up};

  // Bosonic elementary spaces (4 bits)
  auto es_b_x = make_space_boson(4, "x", 0);
  auto es_b_y = make_space_boson(4, "y", 0);
  std::vector<es_type*> boson_es = {&es_b_x,&es_b_y};

  // Spin-1/2 algebra elementary spaces
  auto es_s_i = make_space_spin(0.5, "i", 0);
  auto es_s_j = make_space_spin(0.5, "j", 0);
  std::vector<es_type*> spin_es = {&es_s_i,&es_s_j};

  // Spin-1 algebra elementary spaces
  auto es_s1_i = make_space_spin(1.0, "i", 0);
  auto es_s1_j = make_space_spin(1.0, "j", 0);
  std::vector<es_type*> spin1_es = {&es_s1_i,&es_s1_j};

  // Spin-3/2 algebra elementary spaces
  auto es_s32_i = make_space_spin(3.0/2, "i", 0);
  auto es_s32_j = make_space_spin(3.0/2, "j", 0);
  std::vector<es_type*> spin32_es = {&es_s32_i,&es_s32_j};

  SECTION("fermion") {
    for(auto * es : fermion_es) {
      CHECK(es->algebra_id() == fermion);
      CHECK(es->n_bits() == 1);
    }
    check_equality(fermion_es);
    check_less_greater(fermion_es);
  }

  SECTION("boson") {
    for(auto * es : boson_es) {
      CHECK(es->algebra_id() == boson);
      CHECK(es->n_bits() == 4);
    }
    check_equality(boson_es);
    check_less_greater(boson_es);
  }

  SECTION("spin-1/2") {
    for(auto * es : spin_es) {
      CHECK(es->algebra_id() == spin);
      CHECK(es->n_bits() == 1);
    }
    check_equality(spin_es);
    check_less_greater(spin_es);
  }

  SECTION("spin-1") {
    for(auto * es : spin1_es) {
      CHECK(es->algebra_id() == spin);
      CHECK(es->n_bits() == 2);
    }
    check_equality(spin1_es);
    check_less_greater(spin1_es);
  }

  SECTION("spin-3/2") {
    for(auto * es : spin32_es) {
      CHECK(es->algebra_id() == spin);
      CHECK(es->n_bits() == 2);
    }
    check_equality(spin32_es);
    check_less_greater(spin32_es);
  }

  SECTION("different base spaces") {
    std::vector<es_type*> all_es;
    for(auto const & es_list : {fermion_es,
                                boson_es,
                                spin_es,
                                spin1_es,
                                spin32_es}) {
      for(auto const& es : es_list) all_es.push_back(es);
    }
    check_equality(all_es);
    check_less_greater(all_es);
  }
}
