/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <catch.hpp>

#include <libcommute/expression/generator_fermion.hpp>
#include <libcommute/expression/monomial.hpp>
#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/monomial_action.hpp>
#include <libcommute/loperator/monomial_action_fermion.hpp>

#include "./monomial_action.hpp"

#include <bitset>
#include <utility>
#include <vector>

using namespace libcommute;

TEST_CASE("Action of a fermionic monomial on an index",
          "[monomial_action_fermion]") {

  using namespace static_indices;

  using mon_type = monomial<int>;
  using hs_type = hilbert_space<int>;
  using pad_es_type = elementary_space_padding<int>;
  using ma_type = monomial_action<fermion>;

  constexpr int n_ops = 4;
  constexpr int n_pad_spaces = 2;
  constexpr int n_pad_bits = 2 * n_pad_spaces;

  std::vector<generator_fermion<int>> gens;
  for(int i = 0; i < n_ops; ++i) {
    gens.emplace_back(true, i);
    gens.emplace_back(false, i);
  }

  auto ref_c_dag_c_action =
      [](generator<int> const& g, sv_index_type& index, double& coeff) {
        int ind = std::get<0>(g.indices());
        bool dagger = dynamic_cast<generator_fermion<int> const&>(g).dagger();
        std::bitset<n_ops + n_pad_bits> in_bitset(index);
        if(dagger && in_bitset.test(ind + n_pad_bits)) return false;
        if(!dagger && !in_bitset.test(ind + n_pad_bits)) return false;
        // Count particles
        int n = 0;
        for(int i = 0; i < ind; ++i)
          n += in_bitset[i + n_pad_bits];
        index = dagger ? in_bitset.set(ind + n_pad_bits).to_ulong() :
                         in_bitset.reset(ind + n_pad_bits).to_ulong();
        coeff *= (n % 2 == 0 ? 1 : -1);
        return true;
      };

  hs_type hs;
  for(int i = 0; i < n_pad_spaces; ++i)
    hs.add(pad_es_type(i));
  for(int i = 0; i < n_ops; ++i)
    hs.add(make_space_fermion(i));

  std::vector<sv_index_type> in_index_list(1 << n_ops);
  for(unsigned int i = 0; i < in_index_list.size(); i++)
    in_index_list[i] = (i << n_pad_bits) + (1 << n_pad_bits) - 1;

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
    // NOLINTNEXTLINE(modernize-loop-convert)
    for(unsigned int i = 0; i < gens.size(); ++i) {
      mon_type mon(gens[i]);
      check_monomial_action<ma_type>(mon,
                                     hs,
                                     ref_c_dag_c_action,
                                     in_index_list);
    }
  }
  SECTION("2 operators") {
    for(unsigned int i = 0; i < gens.size(); ++i) {
      for(unsigned int j = 0; j < gens.size(); ++j) {
        if(!(gens[i] < gens[j])) continue;
        mon_type mon(gens[i], gens[j]);
        check_monomial_action<ma_type>(mon,
                                       hs,
                                       ref_c_dag_c_action,
                                       in_index_list);
      }
    }
  }
  SECTION("3 operators") {
    for(unsigned int i = 0; i < gens.size(); ++i) {
      for(unsigned int j = 0; j < gens.size(); ++j) {
        for(unsigned int k = 0; k < gens.size(); ++k) {
          if(!(gens[i] < gens[j]) || !(gens[j] < gens[k])) continue;
          mon_type mon(gens[i], gens[j], gens[k]);
          check_monomial_action<ma_type>(mon,
                                         hs,
                                         ref_c_dag_c_action,
                                         in_index_list);
        }
      }
    }
  }
  SECTION("4 operators") {
    for(unsigned int i = 0; i < gens.size(); ++i) {
      for(unsigned int j = 0; j < gens.size(); ++j) {
        for(unsigned int k = 0; k < gens.size(); ++k) {
          for(unsigned int l = 0; l < gens.size(); ++l) {
            if(!(gens[i] < gens[j]) || !(gens[j] < gens[k]) ||
               !(gens[k] < gens[l]))
              continue;
            mon_type mon(gens[i], gens[j], gens[k], gens[l]);
            check_monomial_action<ma_type>(mon,
                                           hs,
                                           ref_c_dag_c_action,
                                           in_index_list);
          }
        }
      }
    }
  }
}
