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

#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/generator_boson.hpp>
#include <libcommute/expression/generator_spin.hpp>
#include <libcommute/expression/monomial.hpp>
#include <libcommute/qoperator/basis_space_fermion.hpp>
#include <libcommute/qoperator/basis_space_boson.hpp>
#include <libcommute/qoperator/basis_space_spin.hpp>
#include <libcommute/qoperator/hilbert_space.hpp>
#include <libcommute/qoperator/monomial_action_fermion.hpp>
#include <libcommute/qoperator/monomial_action_boson.hpp>
#include <libcommute/qoperator/monomial_action_spin.hpp>

#include "./monomial_action.hpp"

#include <string>
#include <utility>

using namespace libcommute;

TEST_CASE("Different algebra tags", "[monomial_action_tags]") {

  using mon_type = monomial<std::string, int>;
  using hs_type = hilbert_space<std::string, int>;

  SECTION("Missing algebra tags") {
    using namespace static_indices;

    hs_type hs(make_space_fermion("dn", 0),
            make_space_boson(4, "x", 0),
            make_space_spin(1.0, "i", 0)
          );

    SECTION("Constant monomial") {
      mon_type const_m{};
      auto ma = monomial_action<>(std::make_pair(const_m.begin(), const_m.end()),
                                  hs);

      sv_index_type out_index;
      double coeff = 10;
      for(sv_index_type in_index = 0; in_index < (1 << 7); ++in_index) {
        bool r = ma.act(in_index, out_index, coeff);
        CHECK(r);
        CHECK(out_index == in_index);
        CHECK(coeff == 10);
      }

      mon_type m{make_boson(false, "y", 0)};

      using ex_type = unknown_generator<std::string, int>;
      CHECK_THROWS_AS(monomial_action<>(std::make_pair(m.begin(), m.end()), hs),
                      ex_type);
    }

    SECTION("Non-constant monomial") {
      mon_type m1{make_fermion(true, "dn", 0),
                  make_spin(spin_component::plus, "i", 0)};

      CHECK_NOTHROW(monomial_action<fermion, spin>(
        std::make_pair(m1.begin(), m1.end()), hs)
      );

      mon_type m2{make_fermion(true, "dn", 0),
                  make_boson(false, "x", 0),
                  make_spin(spin_component::plus, "i", 0)};
      auto m2_range = std::make_pair(m2.begin(), m2.end());

      using ex_type = unknown_generator<std::string, int>;
      CHECK_THROWS_AS((monomial_action<>(m2_range, hs)), ex_type);
      CHECK_THROWS_AS((monomial_action<fermion>(m2_range, hs)), ex_type);
      CHECK_THROWS_AS((monomial_action<boson>(m2_range, hs)), ex_type);
      CHECK_THROWS_AS((monomial_action<spin>(m2_range, hs)), ex_type);
      CHECK_THROWS_AS((monomial_action<fermion, boson>(m2_range, hs)), ex_type);
      CHECK_THROWS_AS((monomial_action<fermion, spin>(m2_range, hs)), ex_type);
      CHECK_THROWS_AS((monomial_action<boson, spin>(m2_range, hs)), ex_type);
      CHECK_NOTHROW((monomial_action<fermion, boson, spin>(m2_range, hs)));
    }
  }
}
