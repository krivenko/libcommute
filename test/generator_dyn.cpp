/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <catch.hpp>

#include "check_ordering.hpp"
#include "print_matcher.hpp"

#include <libcommute/expression/generator_boson.hpp>
#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/generator_spin.hpp>

#include <algorithm>
#include <initializer_list>
#include <iterator>
#include <vector>

using namespace libcommute;
using namespace dynamic_indices;

using gen_type = generator<dyn_indices>;
using linear_function_t = gen_type::linear_function_t;

// Check if a linear function of generators has only one constant term.
#define CHECK_LINEAR_FUNCTION_0(NAME, CONST_TERM)                              \
  CHECK(NAME.const_term == CONST_TERM);                                        \
  CHECK(NAME.terms.empty());

// Check if a linear function of generators has a constant term and one
// non-constant term.
#define CHECK_LINEAR_FUNCTION_1(NAME, CONST_TERM, COEFF, GEN)                  \
  CHECK(NAME.const_term == CONST_TERM);                                        \
  CHECK(NAME.terms.size() == 1);                                               \
  CHECK(*NAME.terms[0].first == GEN);                                          \
  CHECK(NAME.terms[0].second == COEFF);

void check_generator_spin_swap_with(std::vector<gen_type*> const& v,
                                    bool one_half = false) {
  linear_function_t f;
  for(std::size_t i = 0; i < v.size(); ++i) {
    for(std::size_t j = i + 1; j < v.size(); ++j) {
      double c = v[j]->swap_with(*v[i], f);
      if(one_half) {
        if(j % 3 == 1 && i == j - 1) { // S_- S_+ = 1/2 - S_z
          CHECK(c == 0);
          CHECK_LINEAR_FUNCTION_1(f, 0.5, -1, *v[i + 2]);
        } else if(j % 3 == 2 && i == j - 2) { // S_z S_+ = 1/2 S_+
          CHECK(c == 0);
          CHECK_LINEAR_FUNCTION_1(f, 0, 0.5, *v[i]);
        } else if(j % 3 == 2 && i == j - 1) { // S_z S_- = -1/2 S_-
          CHECK(c == 0);
          CHECK_LINEAR_FUNCTION_1(f, 0, -0.5, *v[i]);
        } else {
          CHECK(c == 1);
          CHECK_LINEAR_FUNCTION_0(f, 0);
        }
      } else {
        CHECK(c == 1);
        if(j % 3 == 1 && i == j - 1) { // S_- S_+ = S_+ * S_- - 2*S_z
          CHECK_LINEAR_FUNCTION_1(f, 0, -2, *v[i + 2]);
        } else if(j % 3 == 2 && i == j - 2) { // S_z S_+ = S_+ * S_z + S_+
          CHECK_LINEAR_FUNCTION_1(f, 0, 1, *v[i]);
        } else if(j % 3 == 2 && i == j - 1) { // S_z S_- = S_- * S_z - S_-
          CHECK_LINEAR_FUNCTION_1(f, 0, -1, *v[i]);
        } else {
          CHECK_LINEAR_FUNCTION_0(f, 0);
        }
      }
    }
  }
}

void check_generator_spin_simplify_prod(std::vector<gen_type*> const& v,
                                        bool one_half = false) {
  linear_function_t f;
  for(std::size_t i = 0; i < v.size(); ++i) {
    for(std::size_t j = i; j < v.size(); ++j) {
      bool c = v[i]->simplify_prod(*v[j], f);
      if(one_half) {
        if(i % 3 == 0 && j == i) { // S_+ * S_+ = 0
          CHECK(c);
          CHECK_LINEAR_FUNCTION_0(f, 0);
        } else if(i % 3 == 1 && j == i) { // S_- * S_- = 0
          CHECK(c);
          CHECK_LINEAR_FUNCTION_0(f, 0);
        } else if(i % 3 == 2 && j == i) { // S_z * S_z = 1/4
          CHECK(c);
          CHECK_LINEAR_FUNCTION_0(f, 0.25);
        } else if(i % 3 == 0 && j == i + 1) { // S_+ * S_- = 1/2 + S_z
          CHECK(c);
          CHECK_LINEAR_FUNCTION_1(f, 0.5, 1, *v[i + 2]);
        } else if(i % 3 == 0 && j == i + 2) { // S_+ * S_z = -1/2 S_+
          CHECK(c);
          CHECK_LINEAR_FUNCTION_1(f, 0, -0.5, *v[i]);
        } else if(i % 3 == 1 && j == i + 1) { // S_- * S_z = 1/2 S_-
          CHECK(c);
          CHECK_LINEAR_FUNCTION_1(f, 0, 0.5, *v[i]);
        } else {
          CHECK_FALSE(c);
        }
      } else {
        // No simplifications for higher spins
        CHECK_FALSE(c);
      }
    }
  }
}

template <typename V>
void check_conj(V const& v, std::initializer_list<int> ref) {
  linear_function_t f;
  int n = 0;
  for(auto ref_it = ref.begin(); ref_it != ref.end(); ++ref_it, ++n) {
    v[n]->conj(f);
    CHECK_LINEAR_FUNCTION_1(f, 0, 1, *v[*ref_it]);
  }
}

TEST_CASE("Algebra generators (dyn_indices)", "[generator]") {
  using namespace dynamic_indices;

  // Setup

  // Fermionic generators
  auto Cdag_dn = make_fermion(true, "dn", 0);
  auto Cdag_up = make_fermion(true, "up", 0);
  auto C_up = make_fermion(false, "up", 0);
  auto C_dn = make_fermion(false, "dn", 0);
  std::vector<gen_type*> fermion_ops = {&Cdag_dn, &Cdag_up, &C_up, &C_dn};

  // Bosonic generators
  auto Adag_x = make_boson(true, "x");
  auto Adag_y = make_boson(true, "y");
  auto A_y = make_boson(false, "y");
  auto A_x = make_boson(false, "x");
  std::vector<gen_type*> boson_ops = {&Adag_x, &Adag_y, &A_y, &A_x};

  // Spin-1/2 algebra generators
  auto Sp_i = make_spin(spin_component::plus, 1);
  auto Sm_i = make_spin(spin_component::minus, 1);
  auto Sz_i = make_spin(spin_component::z, 1);
  auto Sp_j = make_spin(spin_component::plus, 2);
  auto Sm_j = make_spin(spin_component::minus, 2);
  auto Sz_j = make_spin(spin_component::z, 2);
  std::vector<gen_type*> spin_ops = {&Sp_i, &Sm_i, &Sz_i, &Sp_j, &Sm_j, &Sz_j};

  // Spin-1 algebra generators
  auto S1p_i = make_spin(1, spin_component::plus, 1);
  auto S1m_i = make_spin(1, spin_component::minus, 1);
  auto S1z_i = make_spin(1, spin_component::z, 1);
  auto S1p_j = make_spin(1, spin_component::plus, 2);
  auto S1m_j = make_spin(1, spin_component::minus, 2);
  auto S1z_j = make_spin(1, spin_component::z, 2);
  std::vector<gen_type*> spin1_ops =
      {&S1p_i, &S1m_i, &S1z_i, &S1p_j, &S1m_j, &S1z_j};

  // Spin-3/2 algebra generators
  auto S32p_i = make_spin(3.0 / 2, spin_component::plus, 1);
  auto S32m_i = make_spin(3.0 / 2, spin_component::minus, 1);
  auto S32z_i = make_spin(3.0 / 2, spin_component::z, 1);
  auto S32p_j = make_spin(3.0 / 2, spin_component::plus, 2);
  auto S32m_j = make_spin(3.0 / 2, spin_component::minus, 2);
  auto S32z_j = make_spin(3.0 / 2, spin_component::z, 2);
  std::vector<gen_type*> spin32_ops =
      {&S32p_i, &S32m_i, &S32z_i, &S32p_j, &S32m_j, &S32z_j};

  linear_function_t lin_f;

  SECTION("fermion") {
    for(auto* op : fermion_ops) {
      CHECK(op->algebra_id() == fermion);
      CHECK_FALSE(op->reduce_power(3, lin_f));
      CHECK_FALSE(op->reduce_power(4, lin_f));
    }

    for(std::size_t i = 0; i < fermion_ops.size(); ++i) {
      auto const* fermion_gen_p =
          dynamic_cast<generator_fermion<dyn_indices>*>(fermion_ops[i]);
      CHECK(fermion_gen_p->dagger() == (i < 2));
    }

    check_equality(fermion_ops);
    check_less_greater(fermion_ops);

    for(auto const& g : fermion_ops) {
      CHECK(is_fermion(*g));
      CHECK_FALSE(is_boson(*g));
      CHECK_FALSE(is_spin(*g));
    }

    linear_function_t f;
    for(std::size_t i = 0; i < fermion_ops.size(); ++i) {
      for(std::size_t j = i + 1; j < fermion_ops.size(); ++j) {
        double c = fermion_ops[j]->swap_with(*fermion_ops[i], f);
        CHECK(c == -1);
        CHECK_LINEAR_FUNCTION_0(f, ((j == 2 && i == 1) || (j == 3 && i == 0)));
      }
    }

    check_conj(fermion_ops, {3, 2, 1, 0});

    CHECK_THAT(Cdag_dn, Prints<gen_type>("C+(dn,0)"));
    CHECK_THAT(Cdag_up, Prints<gen_type>("C+(up,0)"));
    CHECK_THAT(C_up, Prints<gen_type>("C(up,0)"));
    CHECK_THAT(C_dn, Prints<gen_type>("C(dn,0)"));
  }

  SECTION("boson") {
    for(auto* op : boson_ops) {
      CHECK(op->algebra_id() == boson);
      CHECK_FALSE(op->reduce_power(3, lin_f));
      CHECK_FALSE(op->reduce_power(4, lin_f));
    }

    for(std::size_t i = 0; i < boson_ops.size(); ++i) {
      auto const* boson_gen_p =
          dynamic_cast<generator_boson<dyn_indices>*>(boson_ops[i]);
      CHECK(boson_gen_p->dagger() == (i < 2));
    }

    check_equality(boson_ops);
    check_less_greater(boson_ops);

    for(auto const& g : boson_ops) {
      CHECK_FALSE(is_fermion(*g));
      CHECK(is_boson(*g));
      CHECK_FALSE(is_spin(*g));
    }

    linear_function_t f;
    for(std::size_t i = 0; i < boson_ops.size(); ++i) {
      for(std::size_t j = i + 1; j < boson_ops.size(); ++j) {
        double c = boson_ops[j]->swap_with(*boson_ops[i], f);
        CHECK(c == 1);
        CHECK_LINEAR_FUNCTION_0(f, ((j == 2 && i == 1) || (j == 3 && i == 0)));
      }
    }

    check_conj(boson_ops, {3, 2, 1, 0});

    CHECK_THAT(Adag_x, Prints<gen_type>("A+(x)"));
    CHECK_THAT(Adag_y, Prints<gen_type>("A+(y)"));
    CHECK_THAT(A_y, Prints<gen_type>("A(y)"));
    CHECK_THAT(A_x, Prints<gen_type>("A(x)"));
  }

  SECTION("spin-1/2") {
    for(auto* op : spin_ops) {
      CHECK(op->algebra_id() == spin);

      auto const* spin_gen_p = dynamic_cast<generator_spin<dyn_indices>*>(op);
      if(spin_gen_p->component() == spin_component::z) {
        CHECK_FALSE(op->reduce_power(3, lin_f));
        CHECK_FALSE(op->reduce_power(4, lin_f));
      } else {
        CHECK(op->reduce_power(3, lin_f));
        CHECK(lin_f.vanishing());
        CHECK(op->reduce_power(4, lin_f));
        CHECK(lin_f.vanishing());
      }
      CHECK(spin_gen_p->spin() == 0.5);
      CHECK(spin_gen_p->multiplicity() == 2);
    }

    check_equality(spin_ops);
    check_less_greater(spin_ops);

    for(auto const& g : spin_ops) {
      CHECK_FALSE(is_fermion(*g));
      CHECK_FALSE(is_boson(*g));
      CHECK(is_spin(*g));
    }

    check_generator_spin_swap_with(spin_ops, true);
    check_generator_spin_simplify_prod(spin_ops, true);

    check_conj(spin_ops, {1, 0, 2, 4, 3, 5});

    CHECK_THAT(Sp_i, Prints<gen_type>("S+(1)"));
    CHECK_THAT(Sm_i, Prints<gen_type>("S-(1)"));
    CHECK_THAT(Sz_i, Prints<gen_type>("Sz(1)"));
    CHECK_THAT(Sp_j, Prints<gen_type>("S+(2)"));
    CHECK_THAT(Sm_j, Prints<gen_type>("S-(2)"));
    CHECK_THAT(Sz_j, Prints<gen_type>("Sz(2)"));
  }

  SECTION("spin-1") {
    for(auto* op : spin1_ops) {
      CHECK(op->algebra_id() == spin);

      auto const* spin_gen_p = dynamic_cast<generator_spin<dyn_indices>*>(op);
      if(spin_gen_p->component() == spin_component::z) {
        CHECK_FALSE(op->reduce_power(3, lin_f));
        CHECK_FALSE(op->reduce_power(4, lin_f));
      } else {
        CHECK(op->reduce_power(3, lin_f));
        CHECK_LINEAR_FUNCTION_0(lin_f, 0);
        CHECK(op->reduce_power(4, lin_f));
        CHECK_LINEAR_FUNCTION_0(lin_f, 0);
      }
      CHECK(spin_gen_p->spin() == 1.0);
      CHECK(spin_gen_p->multiplicity() == 3);
    }

    check_equality(spin1_ops);
    check_less_greater(spin_ops);

    for(auto const& g : spin1_ops) {
      CHECK_FALSE(is_fermion(*g));
      CHECK_FALSE(is_boson(*g));
      CHECK(is_spin(*g));
    }

    check_generator_spin_swap_with(spin1_ops);
    check_generator_spin_simplify_prod(spin1_ops);

    check_conj(spin1_ops, {1, 0, 2, 4, 3, 5});

    CHECK_THAT(S1p_i, Prints<gen_type>("S1+(1)"));
    CHECK_THAT(S1m_i, Prints<gen_type>("S1-(1)"));
    CHECK_THAT(S1z_i, Prints<gen_type>("S1z(1)"));
    CHECK_THAT(S1p_j, Prints<gen_type>("S1+(2)"));
    CHECK_THAT(S1m_j, Prints<gen_type>("S1-(2)"));
    CHECK_THAT(S1z_j, Prints<gen_type>("S1z(2)"));
  }

  SECTION("spin-3/2") {
    for(auto* op : spin32_ops) {
      CHECK(op->algebra_id() == spin);

      auto const* spin_gen_p = dynamic_cast<generator_spin<dyn_indices>*>(op);
      if(spin_gen_p->component() == spin_component::z) {
        CHECK_FALSE(op->reduce_power(3, lin_f));
        CHECK_FALSE(op->reduce_power(4, lin_f));
      } else {
        CHECK_FALSE(op->reduce_power(3, lin_f));
        CHECK(op->reduce_power(4, lin_f));
        CHECK_LINEAR_FUNCTION_0(lin_f, 0);
      }
      CHECK(spin_gen_p->spin() == 1.5);
      CHECK(spin_gen_p->multiplicity() == 4);
    }

    check_equality(spin32_ops);
    check_less_greater(spin32_ops);

    for(auto const& g : spin32_ops) {
      CHECK_FALSE(is_fermion(*g));
      CHECK_FALSE(is_boson(*g));
      CHECK(is_spin(*g));
    }

    check_generator_spin_swap_with(spin32_ops);
    check_generator_spin_simplify_prod(spin32_ops);

    check_conj(spin32_ops, {1, 0, 2, 4, 3, 5});

    CHECK_THAT(S32p_i, Prints<gen_type>("S3/2+(1)"));
    CHECK_THAT(S32m_i, Prints<gen_type>("S3/2-(1)"));
    CHECK_THAT(S32z_i, Prints<gen_type>("S3/2z(1)"));
    CHECK_THAT(S32p_j, Prints<gen_type>("S3/2+(2)"));
    CHECK_THAT(S32m_j, Prints<gen_type>("S3/2-(2)"));
    CHECK_THAT(S32z_j, Prints<gen_type>("S3/2z(2)"));
  }

  SECTION("different_algebras") {
    std::vector<gen_type*> all_ops;
    for(auto const& op_list :
        {fermion_ops, boson_ops, spin_ops, spin1_ops, spin32_ops}) {
      std::copy(op_list.begin(), op_list.end(), std::back_inserter(all_ops));
    }

    check_equality(all_ops);
    check_less_greater(all_ops);

    // Check that generators from different algebras commute
    linear_function_t f;
    for(std::size_t i = 0; i < all_ops.size(); ++i) {
      for(std::size_t j = i + 1; j < all_ops.size(); ++j) {
        double c = swap_with(*all_ops[j], *all_ops[i], f);
        if(all_ops[j]->algebra_id() != all_ops[i]->algebra_id()) {
          CHECK(c == 1);
          CHECK_LINEAR_FUNCTION_0(f, 0);
        }
      }
    }
  }
}
