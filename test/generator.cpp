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
#include "utility.hpp"

#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/generator_boson.hpp>
#include <libcommute/expression/generator_spin.hpp>

#include <vector>

using namespace libcommute;

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

  // Spin-1/2 algebra generators
  auto Sp_i = make_spin(spin_component::plus, "i", 0);
  auto Sm_i = make_spin(spin_component::minus, "i", 0);
  auto Sz_i = make_spin(spin_component::z, "i", 0);
  auto Sp_j = make_spin(spin_component::plus, "j", 0);
  auto Sm_j = make_spin(spin_component::minus, "j", 0);
  auto Sz_j = make_spin(spin_component::z, "j", 0);
  std::vector<gen_type*> spin_ops = {&Sp_i,&Sm_i,&Sz_i,&Sp_j,&Sm_j,&Sz_j};

  // Spin-1 algebra generators
  auto S1p_i = make_spin(1, spin_component::plus, "i", 0);
  auto S1m_i = make_spin(1, spin_component::minus, "i", 0);
  auto S1z_i = make_spin(1, spin_component::z, "i", 0);
  auto S1p_j = make_spin(1, spin_component::plus, "j", 0);
  auto S1m_j = make_spin(1, spin_component::minus, "j", 0);
  auto S1z_j = make_spin(1, spin_component::z, "j", 0);
  std::vector<gen_type*> spin1_ops = {&S1p_i,&S1m_i,&S1z_i,
                                      &S1p_j,&S1m_j,&S1z_j};

  // Spin-3/2 algebra generators
  auto S32p_i = make_spin(3.0/2, spin_component::plus, "i", 0);
  auto S32m_i = make_spin(3.0/2, spin_component::minus, "i", 0);
  auto S32z_i = make_spin(3.0/2, spin_component::z, "i", 0);
  auto S32p_j = make_spin(3.0/2, spin_component::plus, "j", 0);
  auto S32m_j = make_spin(3.0/2, spin_component::minus, "j", 0);
  auto S32z_j = make_spin(3.0/2, spin_component::z, "j", 0);
  std::vector<gen_type*> spin32_ops = {&S32p_i,&S32m_i,&S32z_i,
                                       &S32p_j,&S32m_j,&S32z_j};

  SECTION("fermion") {
    for(auto * op : fermion_ops) {
      CHECK(op->algebra_id() == FERMION_ALGEBRA_ID);
      CHECK(op->nilpotent_power() == 2);
    }

    check_equality(fermion_ops);
    check_less_greater(fermion_ops);

    for(auto const& g : fermion_ops) {
      CHECK(is_fermion(*g));
      CHECK_FALSE(is_boson(*g));
      CHECK_FALSE(is_spin(*g));
    }

    CHECK_THAT(Cdag_dn, Prints<gen_type>("C+(dn,0)"));
    CHECK_THAT(Cdag_up, Prints<gen_type>("C+(up,0)"));
    CHECK_THAT(C_up, Prints<gen_type>("C(up,0)"));
    CHECK_THAT(C_dn, Prints<gen_type>("C(dn,0)"));
  }

  SECTION("boson") {
    for(auto * op : boson_ops) {
      CHECK(op->algebra_id() == BOSON_ALGEBRA_ID);
      CHECK(op->nilpotent_power() == -1);
    }

    check_equality(boson_ops);
    check_less_greater(boson_ops);

    for(auto const& g : boson_ops) {
      CHECK_FALSE(is_fermion(*g));
      CHECK(is_boson(*g));
      CHECK_FALSE(is_spin(*g));
    }

    CHECK_THAT(Adag_x, Prints<gen_type>("A+(x,0)"));
    CHECK_THAT(Adag_y, Prints<gen_type>("A+(y,0)"));
    CHECK_THAT(A_y, Prints<gen_type>("A(y,0)"));
    CHECK_THAT(A_x, Prints<gen_type>("A(x,0)"));
  }

  SECTION("spin-1/2") {
    for(auto * op : spin_ops) {
      CHECK(op->algebra_id() == SPIN_ALGEBRA_ID);

      auto spin_gen_p = dynamic_cast<generator_spin<std::string, int>*>(op);
      CHECK((spin_gen_p->component() == spin_component::z ||
            (op->nilpotent_power() == 2)));
    }

    check_equality(spin_ops);
    check_less_greater(spin_ops);

    for(auto const& g : spin_ops) {
      CHECK_FALSE(is_fermion(*g));
      CHECK_FALSE(is_boson(*g));
      CHECK(is_spin(*g));
    }

    CHECK_THAT(Sp_i, Prints<gen_type>("S+(i,0)"));
    CHECK_THAT(Sm_i, Prints<gen_type>("S-(i,0)"));
    CHECK_THAT(Sz_i, Prints<gen_type>("Sz(i,0)"));
    CHECK_THAT(Sp_j, Prints<gen_type>("S+(j,0)"));
    CHECK_THAT(Sm_j, Prints<gen_type>("S-(j,0)"));
    CHECK_THAT(Sz_j, Prints<gen_type>("Sz(j,0)"));
  }

  SECTION("spin-1") {
    for(auto * op : spin1_ops) {
      CHECK(op->algebra_id() == SPIN_ALGEBRA_ID);

      auto spin_gen_p = dynamic_cast<generator_spin<std::string, int>*>(op);
      CHECK((spin_gen_p->component() == spin_component::z ||
            (op->nilpotent_power() == 3)));
    }

    check_equality(spin1_ops);
    check_less_greater(spin_ops);

    for(auto const& g : spin1_ops) {
      CHECK_FALSE(is_fermion(*g));
      CHECK_FALSE(is_boson(*g));
      CHECK(is_spin(*g));
    }

    CHECK_THAT(S1p_i, Prints<gen_type>("S1+(i,0)"));
    CHECK_THAT(S1m_i, Prints<gen_type>("S1-(i,0)"));
    CHECK_THAT(S1z_i, Prints<gen_type>("S1z(i,0)"));
    CHECK_THAT(S1p_j, Prints<gen_type>("S1+(j,0)"));
    CHECK_THAT(S1m_j, Prints<gen_type>("S1-(j,0)"));
    CHECK_THAT(S1z_j, Prints<gen_type>("S1z(j,0)"));
  }

  SECTION("spin-3/2") {
    for(auto * op : spin32_ops) {
      CHECK(op->algebra_id() == SPIN_ALGEBRA_ID);

      auto spin_gen_p = dynamic_cast<generator_spin<std::string, int>*>(op);
      CHECK((spin_gen_p->component() == spin_component::z ||
            (op->nilpotent_power() == 4)));
    }

    check_equality(spin32_ops);
    check_less_greater(spin32_ops);

    for(auto const& g : spin32_ops) {
      CHECK_FALSE(is_fermion(*g));
      CHECK_FALSE(is_boson(*g));
      CHECK(is_spin(*g));
    }

    CHECK_THAT(S32p_i, Prints<gen_type>("S3/2+(i,0)"));
    CHECK_THAT(S32m_i, Prints<gen_type>("S3/2-(i,0)"));
    CHECK_THAT(S32z_i, Prints<gen_type>("S3/2z(i,0)"));
    CHECK_THAT(S32p_j, Prints<gen_type>("S3/2+(j,0)"));
    CHECK_THAT(S32m_j, Prints<gen_type>("S3/2-(j,0)"));
    CHECK_THAT(S32z_j, Prints<gen_type>("S3/2z(j,0)"));
  }

  SECTION("different_algebras") {
    std::vector<gen_type*> all_ops;
    for(auto const & op_list : {fermion_ops,
                                boson_ops,
                                spin_ops,
                                spin1_ops,
                                spin32_ops}) {
      for(auto const& op : op_list) all_ops.push_back(op);
    }

    check_equality(all_ops);
    check_less_greater(all_ops);
  }
}
