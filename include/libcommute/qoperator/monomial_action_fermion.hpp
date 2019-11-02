/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_QOPERATOR_MONOMIAL_ACTION_FERMION_HPP_
#define LIBCOMMUTE_QOPERATOR_MONOMIAL_ACTION_FERMION_HPP_

#include "basis_space_fermion.hpp"
#include "hilbert_space.hpp"
#include "monomial_action.hpp"
#include "state_vector.hpp"
#include "../expression/generator_fermion.hpp"

#include <algorithm>
#include <vector>

//
// Action of a monomial comprised of fermionic algebra generators
//

namespace libcommute {

template<> class monomial_action<fermion> {

  // Is this monomial a constant?
  bool is_const = false;
  // Bit masks used to change bits
  sv_index_type annihilation_mask = 0;
  sv_index_type creation_mask = 0;
  // Bit masks for particle counting
  sv_index_type annihilation_count_mask;
  sv_index_type creation_count_mask;

public:

  template<typename... IndexTypes>
  monomial_action(detail::monomial_range_t<IndexTypes...> const& m_range,
                  hilbert_space<IndexTypes...> const& hs) {

    if(m_range.second == m_range.first) {
      is_const = true;
      return;
    }

    std::vector<int> creation_set_bits;
    std::vector<int> annihilation_set_bits;

    for(auto it = m_range.first; it != m_range.second; ++it) {
      if(!is_fermion(*it))
        throw unknown_generator<IndexTypes...>(*it);

      basis_space_fermion<IndexTypes...> bs(it->indices());
      if(!hs.has(bs))
        throw unknown_generator<IndexTypes...>(*it);

      auto br = hs.bit_range(bs);
      // All fermionic basis spaces are 2-dimensional
      assert(br.first == br.second);

      bool dagger =
        dynamic_cast<generator_fermion<IndexTypes...> const&>(*it).dagger();

      (dagger ? creation_set_bits : annihilation_set_bits)
        .emplace_back(br.first);
      (dagger ? creation_mask : annihilation_mask)
        |= (sv_index_type(1) << br.first);
    }

    auto const& range = hs.algebra_bit_range(fermion::algebra_id());

    creation_count_mask = compute_count_mask(creation_set_bits, range);
    annihilation_count_mask = compute_count_mask(annihilation_set_bits, range);
  }

  template<typename ScalarType>
  inline bool act(sv_index_type & index,
                  ScalarType & coeff) const {

    if(is_const) return true;

    // Fermions
    if ((index & annihilation_mask) != annihilation_mask)
      return false; // Zero after acting with the annihilation operators

    sv_index_type inter_index = index & ~annihilation_mask;

    if (((inter_index ^ creation_mask) & creation_mask) != creation_mask)
      return false; // Zero after acting with the creation operators

    index = ~(~inter_index & ~creation_mask);
    bool minus = parity_number_of_bits((inter_index & annihilation_count_mask) ^
                                       (index & creation_count_mask)
                                      );
    if(minus) coeff *= -1;
    return true;
  }

private:

  // Compute parity of the number of set bits in i
  inline static bool parity_number_of_bits(sv_index_type i) {
    i ^= i >> 32;
    i ^= i >> 16;
    i ^= i >> 8;
    i ^= i >> 4;
    i ^= i >> 2;
    i ^= i >> 1;
    return i & 0x01;
  }

  //
  inline static sv_index_type
  compute_count_mask(std::vector<int> const &d, bit_range_t const& bit_range) {
    sv_index_type mask = 0;
    bool is_on = (d.size() % 2 == 1);
    for(int i = bit_range.first; i <= bit_range.second; ++i) {
      if(std::find(d.begin(), d.end(), i) != d.end())
        is_on = !is_on;
      else if(is_on)
        mask |= (sv_index_type(1) << i);
    }
    return mask;
  }
};

} // namespace libcommute

#endif
