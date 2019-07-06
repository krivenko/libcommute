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

#include <libcommute/generator.hpp>

#include <sstream>
#include <vector>

using namespace libcommute;

template<typename V> void check_equality(V const& v) {
  for(int i1 = 0; i1 < v.size(); ++i1) {
    for(int i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] == *v[i2]) == (i1 == i2));
      CHECK((*v[i1] != *v[i2]) == (i1 != i2));
    }
  }
}

template<typename V> void check_less_greater(V const& v) {
  for(int i1 = 0; i1 < v.size(); ++i1) {
    for(int i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] < *v[i2]) == (i1 < i2));
      CHECK((*v[i1] > *v[i2]) == (i1 > i2));
    }
  }
}

template<typename OP> void check_print(OP const& op, std::string const& ref) {
  std::stringstream ss;
  ss << op;
  CHECK(ss.str() == ref);
}

TEST_CASE("Algebra generators", "[generator]") {

  // Setup
  using gen_type = generator<std::string, int>;

  // Fermionic generators
  auto Cdag_dn = make_fermion(true, "dn", 0);
  auto Cdag_up = make_fermion(true, "up", 0);
  auto C_up = make_fermion(false, "up", 0);
  auto C_dn = make_fermion(false, "dn", 0);
  std::vector<gen_type*> fermion_ops = {&Cdag_dn,&Cdag_up,&C_up,&C_dn};

  // Bosonic generators
  auto Adag_x = make_boson(true, "x", 0);
  auto Adag_y = make_boson(true, "y", 0);
  auto A_y = make_boson(false, "y", 0);
  auto A_x = make_boson(false, "x", 0);
  std::vector<gen_type*> boson_ops = {&Adag_x,&Adag_y,&A_y,&A_x};

  // Spin algebra generators
  auto Sp_i = make_spin(spin_component::plus, "i", 0);
  auto Sm_i = make_spin(spin_component::minus, "i", 0);
  auto Sz_i = make_spin(spin_component::z, "i", 0);
  auto Sp_j = make_spin(spin_component::plus, "j", 0);
  auto Sm_j = make_spin(spin_component::minus, "j", 0);
  auto Sz_j = make_spin(spin_component::z, "j", 0);
  std::vector<gen_type*> spin_ops = {&Sp_i,&Sm_i,&Sz_i,&Sp_j,&Sm_j,&Sz_j};

  SECTION("fermion") {
    for(auto * op : fermion_ops)
      CHECK(op->algebra_id() == FERMION_ALGEBRA_ID);

    check_equality(fermion_ops);
    check_less_greater(fermion_ops);

    check_print(Cdag_dn, "C+(dn,0)");
    check_print(Cdag_up, "C+(up,0)");
    check_print(C_up, "C(up,0)");
    check_print(C_dn, "C(dn,0)");
  }

  SECTION("boson") {
    for(auto * op : boson_ops)
      CHECK(op->algebra_id() == BOSON_ALGEBRA_ID);

    check_equality(boson_ops);
    check_less_greater(boson_ops);

    check_print(Adag_x, "A+(x,0)");
    check_print(Adag_y, "A+(y,0)");
    check_print(A_y, "A(y,0)");
    check_print(A_x, "A(x,0)");
  }

  SECTION("spin") {
    for(auto * op : spin_ops)
      CHECK(op->algebra_id() == SPIN_ALGEBRA_ID);

    check_equality(spin_ops);
    check_less_greater(spin_ops);

    check_print(Sp_i, "S+(i,0)");
    check_print(Sm_i, "S-(i,0)");
    check_print(Sz_i, "Sz(i,0)");
    check_print(Sp_j, "S+(j,0)");
    check_print(Sm_j, "S-(j,0)");
    check_print(Sz_j, "Sz(j,0)");
  }

  SECTION("different_algebras") {
    std::vector<gen_type*> all_ops;
    for(auto const& op : fermion_ops) all_ops.push_back(op);
    for(auto const& op : boson_ops) all_ops.push_back(op);
    for(auto const& op : spin_ops) all_ops.push_back(op);

    check_equality(all_ops);
    check_less_greater(all_ops);
  }
}
