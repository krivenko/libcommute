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

#include <libcommute/expression/generator_boson.hpp>
#include <libcommute/expression/monomial.hpp>
#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/monomial_action.hpp>
#include <libcommute/loperator/monomial_action_boson.hpp>

#include "./monomial_action.hpp"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <utility>
#include <vector>

using namespace libcommute;

TEST_CASE("Action of a bosonic monomial on an index",
          "[monomial_action_fermion]") {

  using namespace static_indices;

  using mon_type = monomial<int>;
  using hs_type = hilbert_space<int>;
  using pad_es_type = elementary_space_padding<int>;
  using ma_type = monomial_action<boson>;

  constexpr int n_ops = 3;
  constexpr int n_op_bits = 6;
  constexpr int n_pad_spaces = 2;
  constexpr int n_pad_bits = 2*n_pad_spaces;
  constexpr int total_n_bits = n_op_bits + n_pad_bits;

  std::vector<generator_boson<int>> gens;
  for(int i = 0; i < n_ops; ++i) {
    gens.emplace_back(true, i);
    gens.emplace_back(false, i);
  }

  hs_type hs;
  for(int i = 0; i < n_pad_spaces; ++i) hs.add(pad_es_type(i));
  std::vector<bit_range_t> bit_ranges = {{4, 4}, {5, 7}, {8, 9}};
  for(int i = 0; i < n_ops; ++i) {
    int n_bits = bit_ranges[i].second - bit_ranges[i].first + 1;
    hs.add(make_space_boson(n_bits, i));
  }

  auto ref_a_dag_a_action = [&](generator<int> const& g,
                                sv_index_type & index,
                                double & coeff) {
    int ind = std::get<0>(g.indices());
    bool dagger = dynamic_cast<generator_boson<int> const&>(g).dagger();

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

  std::vector<sv_index_type> in_index_list(1 << n_op_bits);
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
    for(unsigned int i = 0; i < gens.size(); ++i) {
      mon_type mon(gens[i]);
      check_monomial_action<ma_type>(mon,
                                     hs,
                                     ref_a_dag_a_action,
                                     in_index_list);
    }
  }
  SECTION("2 operators") {
    for(unsigned int i = 0; i < gens.size(); ++i) {
      for(unsigned int j = 0; j < gens.size(); ++j) {
        if(gens[i] > gens[j]) continue;
        mon_type mon(gens[i], gens[j]);
        check_monomial_action<ma_type>(mon,
                                       hs,
                                       ref_a_dag_a_action,
                                       in_index_list);
      }
    }
  }
  SECTION("3 operators") {
    for(unsigned int i = 0; i < gens.size(); ++i) {
      for(unsigned int j = 0; j < gens.size(); ++j) {
        for(unsigned int k = 0; k < gens.size(); ++k) {
          if(gens[i] > gens[j] || gens[j] > gens[k]) continue;
          mon_type mon(gens[i], gens[j], gens[k]);
          check_monomial_action<ma_type>(mon,
                                         hs,
                                         ref_a_dag_a_action,
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
            if(gens[i] > gens[j] ||
               gens[j] > gens[k] ||
               gens[k] > gens[l]) continue;
            mon_type mon(gens[i], gens[j], gens[k], gens[l]);
            check_monomial_action<ma_type>(mon,
                                           hs,
                                           ref_a_dag_a_action,
                                           in_index_list);
          }
        }
      }
    }
  }
}
