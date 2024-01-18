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
#ifndef LIBCOMMUTE_LOPERATOR_MONOMIAL_ACTION_SPIN_HPP_
#define LIBCOMMUTE_LOPERATOR_MONOMIAL_ACTION_SPIN_HPP_

#include "../expression/generator_spin.hpp"
#include "elementary_space_spin.hpp"
#include "hilbert_space.hpp"
#include "monomial_action.hpp"
#include "state_vector.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

//
// Action of a monomial comprised of spin algebra generators
//

#ifndef LIBCOMMUTE_SPIN_MAX_NUM_PRECOMPUTED_SQRT
// Maximum allowed size of `sqr_roots` in monomial_action<spin>
#define LIBCOMMUTE_SPIN_MAX_NUM_PRECOMPUTED_SQRT 128
#endif

namespace libcommute {

template <> class monomial_action<spin> {

  //
  // Calculations in this class are perfomed using the shifted magnetic quantum
  // number n = m + s, n = 0,...,2*s
  //

  // Update of a single spin mode
  struct single_spin_update_t {
    // Spin times 2
    sv_index_type s2;

    // The bit range corresponding to one spin mode is selected as
    // (in_index >> shift) & mask
    int shift;
    sv_index_type mask;

    // Generator: +, - or z
    spin_component c;

    // Power of generator
    sv_index_type power;
  };

  // List of single-boson updates
  std::vector<single_spin_update_t> updates_;

  // Precomputed square roots of integers
  std::vector<double> sqr_roots_;

  inline double sqr_root(sv_index_type n) const {
    if(n < LIBCOMMUTE_SPIN_MAX_NUM_PRECOMPUTED_SQRT)
      return sqr_roots_[n];
    else
      return std::sqrt(n);
  }

public:
  template <typename... IndexTypes>
  monomial_action(detail::monomial_range_t<IndexTypes...> const& m_range,
                  hilbert_space<IndexTypes...> const& hs) {
    sv_index_type sqr_roots_size = 0;

    auto it = m_range.first;
    auto next_it = it;
    ++next_it;
    auto end_it = m_range.second;

    sv_index_type power = 1;
    for(; it != end_it; ++it, ++next_it) {
      if(!is_spin(*it)) throw unknown_generator<IndexTypes...>(*it);

      if(next_it == end_it || *next_it != *it) {
        auto const& g = dynamic_cast<generator_spin<IndexTypes...> const&>(*it);
        double s = g.spin();

        elementary_space_spin<IndexTypes...> es(s, g.indices());
        if(!hs.has(es)) throw unknown_generator<IndexTypes...>(g);

        bit_range_t const& bit_range = hs.bit_range(es);
        int shift = bit_range.first;
        int n_bits = bit_range.second - bit_range.first + 1;

        sv_index_type ss =
            g.multiplicity() % 2 == 0 ? (s + 0.5) * (s + 0.5) : s * (s + 1);
        sqr_roots_size = std::max(sqr_roots_size, ss + 1);

        updates_.emplace_back(
            single_spin_update_t{sv_index_type(2 * s),
                                 shift,
                                 (sv_index_type(1) << n_bits) - 1,
                                 g.component(),
                                 power});

        power = 1;
      } else
        ++power;
    }

    if(!updates_.empty()) {
      sqr_roots_size =
          std::min(sqr_roots_size,
                   sv_index_type(LIBCOMMUTE_SPIN_MAX_NUM_PRECOMPUTED_SQRT));
      sqr_roots_.resize(sqr_roots_size);
      for(sv_index_type n = 0; n < sqr_roots_size; ++n)
        sqr_roots_[n] = std::sqrt(double(n));
    }
  }

  template <typename ScalarType>
  inline bool act(sv_index_type& index, ScalarType& coeff) const {

    for(std::size_t b = updates_.size(); b-- != 0;) {
      auto const& update = updates_[b];
      sv_index_type n = (index >> update.shift) & update.mask;
      switch(update.c) {
      case plus: {
        if(n + update.power > update.s2) return false;
        for(sv_index_type d = 0; d < update.power; ++d)
          mul_assign(coeff,
                     scalar_traits<ScalarType>::make_const(
                         sqr_root((update.s2 - (n + d)) * (n + d + 1))));
        index += update.power << update.shift;
      } break;
      case minus: {
        if(n < update.power) return false;
        for(sv_index_type d = 0; d < update.power; ++d)
          mul_assign(coeff,
                     scalar_traits<ScalarType>::make_const(
                         sqr_root((update.s2 - (n - d) + 1) * (n - d))));
        index -= update.power << update.shift;
      } break;
      case z: {
        if((update.s2 % 2 == 0) && n == update.s2 / 2) return false;
        mul_assign(
            coeff,
            scalar_traits<ScalarType>::make_const(
                std::pow(double(n) - double(update.s2) / 2, update.power)));
      } break;
      }
    }
    return true;
  }
};

} // namespace libcommute

#endif
