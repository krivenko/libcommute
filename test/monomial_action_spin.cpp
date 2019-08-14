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

#include <libcommute/expression/generator_spin.hpp>
#include <libcommute/expression/monomial.hpp>
#include <libcommute/qoperator/hilbert_space.hpp>
#include <libcommute/qoperator/monomial_action.hpp>
#include <libcommute/qoperator/monomial_action_spin.hpp>

#include "./monomial_action.hpp"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <utility>
#include <vector>

using namespace libcommute;

TEST_CASE("Action of a spin monomial on an index",
          "[monomial_action_spin]") {

  using namespace static_indices;

  using mon_type = monomial<int>;
  using hs_type = hilbert_space<int>;
  using pad_bs_type = basis_space_padding<int>;
  using ma_type = monomial_action<spin>;

  constexpr int n_op_bits = 5;
  constexpr int n_pad_spaces = 2;
  constexpr int n_pad_bits = 2*n_pad_spaces;
  constexpr int total_n_bits = n_op_bits + n_pad_bits;

  std::vector<generator_spin<int>> gens;
  gens.emplace_back(0.5, spin_component::plus, 0);
  gens.emplace_back(0.5, spin_component::minus, 0);
  gens.emplace_back(0.5, spin_component::z, 0);
  gens.emplace_back(1.0, spin_component::plus, 1);
  gens.emplace_back(1.0, spin_component::minus, 1);
  gens.emplace_back(1.0, spin_component::z, 1);
  gens.emplace_back(1.5, spin_component::plus, 2);
  gens.emplace_back(1.5, spin_component::minus, 2);
  gens.emplace_back(1.5, spin_component::z, 2);

  hs_type hs;
  for(int i = 0; i < n_pad_spaces; ++i) hs.add(pad_bs_type(i));
  std::vector<bit_range_t> bit_ranges = {{4, 4}, {5, 6}, {7, 8}};
  hs.add(make_space_spin(0.5, 0));
  hs.add(make_space_spin(1.0, 1));
  hs.add(make_space_spin(1.5, 2));

  auto ref_spin_action = [&](generator<int> const& g,
                             sv_index_type in_index,
                             sv_index_type & out_index,
                             double & coeff) {
    int ind = std::get<0>(g.indices());
    double s = dynamic_cast<generator_spin<int> const&>(g).spin();
    spin_component c = dynamic_cast<generator_spin<int> const&>(g).component();

    auto const& bit_range = bit_ranges[ind];
    int n_bits = bit_range.second - bit_range.first + 1;

    std::bitset<total_n_bits> in_bitset(in_index);

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

    std::bitset<total_n_bits> out_bitset(in_index);
    std::bitset<total_n_bits> n_bitset((unsigned int)(m + s));
    for(int i = 0; i < n_bits; ++i) {
      out_bitset[bit_range.first + i] = n_bitset[i];
    }
    out_index = out_bitset.to_ulong();

    return true;
  };

  std::vector<sv_index_type> in_index_list;
  for(sv_index_type n12 : {0, 1}) {
    for(sv_index_type n1 : {0, 1, 2}) {
      for(sv_index_type n32 : {0, 1, 2, 3}) {
        sv_index_type index = (1 << n_pad_bits) - 1;
        index += n12 << n_pad_bits;
        index += n1 << (n_pad_bits + 1);
        index += n32 << (n_pad_bits + 3);
        in_index_list.push_back(index);
      }
    }
  }

  SECTION("No operators") {
    mon_type mon;
    ma_type ma(std::make_pair(mon.begin(), mon.end()), hs);
    for(auto in_index : in_index_list) {
      double coeff = 2;
      sv_index_type out_index;
      bool nonzero = ma.act(in_index, out_index, coeff);
      CHECK(nonzero);
      CHECK(out_index == in_index);
      CHECK(coeff == 2);
    }
  }
  SECTION("1 operator") {
    for(int i = 0; i < gens.size(); ++i) {
      mon_type mon(gens[i]);
      check_monomial_action<ma_type>(mon,
                                     hs,
                                     ref_spin_action,
                                     in_index_list);
    }
  }
  SECTION("2 operators") {
    for(int i = 0; i < gens.size(); ++i) {
      for(int j = 0; j < gens.size(); ++j) {
        if(gens[i] > gens[j]) continue;
        mon_type mon(gens[i], gens[j]);
        check_monomial_action<ma_type>(mon,
                                       hs,
                                       ref_spin_action,
                                       in_index_list);
      }
    }
  }
  SECTION("3 operators") {
    for(int i = 0; i < gens.size(); ++i) {
      for(int j = 0; j < gens.size(); ++j) {
        for(int k = 0; k < gens.size(); ++k) {
          if(gens[i] > gens[j] || gens[j] > gens[k]) continue;
          mon_type mon(gens[i], gens[j], gens[k]);
          check_monomial_action<ma_type>(mon,
                                         hs,
                                         ref_spin_action,
                                         in_index_list);
        }
      }
    }
  }
  SECTION("4 operators") {
    for(int i = 0; i < gens.size(); ++i) {
      for(int j = 0; j < gens.size(); ++j) {
        for(int k = 0; k < gens.size(); ++k) {
          for(int l = 0; l < gens.size(); ++l) {
            if(gens[i] > gens[j] ||
               gens[j] > gens[k] ||
               gens[k] > gens[l]) continue;
            mon_type mon(gens[i], gens[j], gens[k], gens[l]);
            check_monomial_action<ma_type>(mon,
                                           hs,
                                           ref_spin_action,
                                           in_index_list);
          }
        }
      }
    }
  }
}
