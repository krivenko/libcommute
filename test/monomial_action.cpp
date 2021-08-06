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

#include <catch.hpp>

#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/generator_boson.hpp>
#include <libcommute/expression/generator_spin.hpp>
#include <libcommute/expression/monomial.hpp>
#include <libcommute/loperator/elementary_space_fermion.hpp>
#include <libcommute/loperator/elementary_space_boson.hpp>
#include <libcommute/loperator/elementary_space_spin.hpp>
#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/monomial_action_fermion.hpp>
#include <libcommute/loperator/monomial_action_boson.hpp>
#include <libcommute/loperator/monomial_action_spin.hpp>

#include "./monomial_action.hpp"

#include <bitset>
#include <memory>
#include <string>
#include <utility>

using namespace libcommute;

TEST_CASE("Different algebra IDs", "[monomial_action_IDs]") {

  using mon_type = monomial<std::string, int>;
  using hs_type = hilbert_space<std::string, int>;

  SECTION("Missing algebra IDs") {
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
        out_index = in_index;
        bool r = ma.act(out_index, coeff);
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
                  make_spin(1, spin_component::plus, "i", 0)};

      CHECK_NOTHROW(monomial_action<fermion, spin>(
        std::make_pair(m1.begin(), m1.end()), hs)
      );

      mon_type m2{make_fermion(true, "dn", 0),
                  make_boson(false, "x", 0),
                  make_spin(1, spin_component::plus, "i", 0)};
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

TEST_CASE("Action of a mixed monomial", "[monomial_action]") {

  using mon_type = monomial<int>;
  using hs_type = hilbert_space<int>;
  using ma_type = monomial_action<fermion, boson, spin>;

  constexpr int total_n_bits = 11;

  using namespace static_indices;

  //
  // Hilbert space
  //

  hs_type hs(make_space_fermion(0),
             make_space_fermion(1),
             make_space_boson(2, 0),
             make_space_boson(2, 1),
             make_space_spin(0.5, 0),
             make_space_spin(1.0, 1),
             make_space_spin(1.5, 2)
        );

  //
  // Algebra generators
  //

  std::vector<std::unique_ptr<generator<int>>> gens;
  gens.emplace_back(make_fermion(true, 0).clone());
  gens.emplace_back(make_fermion(false, 0).clone());
  gens.emplace_back(make_fermion(true, 1).clone());
  gens.emplace_back(make_fermion(false, 1).clone());
  gens.emplace_back(make_boson(true, 0).clone());
  gens.emplace_back(make_boson(false, 0).clone());
  gens.emplace_back(make_boson(true, 1).clone());
  gens.emplace_back(make_boson(false, 1).clone());
  gens.emplace_back(make_spin(0.5, spin_component::plus, 0).clone());
  gens.emplace_back(make_spin(0.5, spin_component::minus, 0).clone());
  gens.emplace_back(make_spin(0.5, spin_component::z, 0).clone());
  gens.emplace_back(make_spin(1.0, spin_component::plus, 1).clone());
  gens.emplace_back(make_spin(1.0, spin_component::minus, 1).clone());
  gens.emplace_back(make_spin(1.0, spin_component::z, 1).clone());
  gens.emplace_back(make_spin(1.5, spin_component::plus, 2).clone());
  gens.emplace_back(make_spin(1.5, spin_component::minus, 2).clone());
  gens.emplace_back(make_spin(1.5, spin_component::z, 2).clone());

  //
  // Reference monomial actions
  //

  // Fermions
  auto ref_c_dag_c_action = [](generator<int> const& g,
                              sv_index_type & index,
                              double & coeff) {
    int ind = std::get<0>(g.indices());
    bool dagger = dynamic_cast<generator_fermion<int> const&>(g).dagger();
    std::bitset<total_n_bits> in_bitset(index);
    if(dagger && in_bitset.test(ind)) return false;
    if(!dagger && !in_bitset.test(ind)) return false;

    int n = 0;
    for(int i = 0; i < ind; ++i) n += in_bitset[i];
    index = dagger ? in_bitset.set(ind).to_ulong() :
                     in_bitset.reset(ind).to_ulong();
    coeff *= (n%2 == 0 ? 1 : -1);
    return true;
  };

  // Bosons
  auto ref_a_dag_a_action = [&](generator<int> const& g,
                              sv_index_type & index,
                              double & coeff) {
    int ind = std::get<0>(g.indices());
    bool dagger = dynamic_cast<generator_boson<int> const&>(g).dagger();

    static std::vector<bit_range_t> bit_ranges = {{2, 3}, {4, 5}};

    auto const& bit_range = bit_ranges[ind];
    int n_bits = bit_range.second - bit_range.first + 1;
    int n_max = (1 << n_bits) - 1;

    std::bitset<total_n_bits> in_bitset(index);

    int n = 0;
    for(int i = 0; i < n_bits; ++i) {
      n += in_bitset.test(bit_range.first + i) * (1 << i);
    }

    if(dagger) {
      ++n;
      if(n > n_max) return false;
      coeff *= std::sqrt(n);
    } else {
      --n;
      if(n < 0) return false;
      coeff *= std::sqrt(n+1);
    }

    std::bitset<total_n_bits> out_bitset(index);
    std::bitset<total_n_bits> n_bitset((unsigned int)n);
    for(int i = 0; i < n_bits; ++i) {
      out_bitset[bit_range.first + i] = n_bitset[i];
    }
    index = out_bitset.to_ulong();

    return true;
  };

  // Spins
  auto ref_spin_action = [&](generator<int> const& g,
                             sv_index_type & index,
                             double & coeff) {
    int ind = std::get<0>(g.indices());
    double s = dynamic_cast<generator_spin<int> const&>(g).spin();
    spin_component c = dynamic_cast<generator_spin<int> const&>(g).component();

    static std::vector<bit_range_t> bit_ranges = {{6, 6}, {7, 8}, {9, 10}};

    auto const& bit_range = bit_ranges[ind];
    int n_bits = bit_range.second - bit_range.first + 1;

    std::bitset<total_n_bits> in_bitset(index);

    double m = -s;
    for(int i = 0; i < n_bits; ++i) {
      m += in_bitset.test(bit_range.first + i) << i;
    }

    switch(c) {
      case plus:
        m += 1;
        if(m > s) return false;
        coeff *= std::sqrt(s*(s+1) - (m-1)*m);
        break;
      case minus:
        m -= 1;
        if(m < -s) return false;
        coeff *= std::sqrt(s*(s+1) - (m+1)*m);
        break;
      case z:
        if(m == 0) return false;
        coeff *= m;
        break;
    }

    std::bitset<total_n_bits> out_bitset(index);
    std::bitset<total_n_bits> n_bitset((unsigned int)(m + s));
    for(int i = 0; i < n_bits; ++i) {
      out_bitset[bit_range.first + i] = n_bitset[i];
    }
    index = out_bitset.to_ulong();

    return true;
  };

  // General case

  auto ref_action = [&](generator<int> const& g,
                        sv_index_type & index,
                        double & coeff) {
    switch(g.algebra_id()) {
      case fermion:
        return ref_c_dag_c_action(g, index, coeff);
      case boson:
        return ref_a_dag_a_action(g, index, coeff);
      case spin:
        return ref_spin_action(g, index, coeff);
    }
    return false;
  };

  //
  // Build list of state indices
  //

  std::vector<sv_index_type> in_index_list;
  for(sv_index_type f1 : {0, 1}) {
    for(sv_index_type f2 : {0, 1}) {
      for(sv_index_type b1 : {0, 1, 2, 3}) {
        for(sv_index_type b2 : {0, 1, 2, 3}) {
          for(sv_index_type s12 : {0, 1}) {
            for(sv_index_type s1 : {0, 1, 2}) {
              for(sv_index_type s32 : {0, 1, 2, 3}) {
                sv_index_type index = 0;
                index += f1;
                index += f2 << 1;
                index += b1 << 2;
                index += b2 << 4;
                index += s12 << 6;
                index += s1 << 7;
                index += s32 << 9;
                in_index_list.push_back(index);
              }
            }
          }
        }
      }
    }
  }

  SECTION("No operators") {
    mon_type mon;
    ma_type ma(std::make_pair(mon.begin(), mon.end()), hs);
    for(auto in_index : in_index_list) {
      double coeff = 2;
      sv_index_type out_index = in_index;
      bool nonzero = ma.act(out_index, coeff);
      CHECK(nonzero);
      CHECK(out_index == in_index);
      CHECK(coeff == 2);
    }
  }
  SECTION("1 operator") {
    for(unsigned int i = 0; i < gens.size(); ++i) {
      mon_type mon{gens[i]->clone()};
      check_monomial_action<ma_type>(mon, hs, ref_action, in_index_list);
    }
  }
  SECTION("2 operators") {
    for(unsigned int i = 0; i < gens.size(); ++i) {
      for(unsigned int j = 0; j < gens.size(); ++j) {
        if(!(*gens[i] < *gens[j])) continue;
        mon_type mon{gens[i]->clone(), gens[j]->clone()};
        check_monomial_action<ma_type>(mon, hs, ref_action, in_index_list);
      }
    }
  }
  SECTION("3 operators") {
    for(unsigned int i = 0; i < gens.size(); ++i) {
      for(unsigned int j = 0; j < gens.size(); ++j) {
        for(unsigned int k = 0; k < gens.size(); ++k) {
          if(!(*gens[i] < *gens[j]) || !(*gens[j] < *gens[k])) continue;
          mon_type mon{gens[i]->clone(), gens[j]->clone(), gens[k]->clone()};
          check_monomial_action<ma_type>(mon, hs, ref_action, in_index_list);
        }
      }
    }
  }
  SECTION("4 operators") {
    for(unsigned int i = 0; i < gens.size(); ++i) {
      for(unsigned int j = 0; j < gens.size(); ++j) {
        for(unsigned int k = 0; k < gens.size(); ++k) {
          for(unsigned int l = 0; l < gens.size(); ++l) {
            if(!(*gens[i] < *gens[j]) ||
               !(*gens[j] < *gens[k]) ||
               !(*gens[k] < *gens[l])) continue;
            mon_type mon{gens[i]->clone(),
                         gens[j]->clone(),
                         gens[k]->clone(),
                         gens[l]->clone()};
            check_monomial_action<ma_type>(mon, hs, ref_action, in_index_list);
          }
        }
      }
    }
  }
}
