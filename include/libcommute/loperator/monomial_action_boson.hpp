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
#ifndef LIBCOMMUTE_LOPERATOR_MONOMIAL_ACTION_BOSON_HPP_
#define LIBCOMMUTE_LOPERATOR_MONOMIAL_ACTION_BOSON_HPP_

#include "../expression/generator_boson.hpp"
#include "elementary_space_boson.hpp"
#include "hilbert_space.hpp"
#include "monomial_action.hpp"
#include "state_vector.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

//
// Action of a monomial comprised of bosonic algebra generators
//

#ifndef LIBCOMMUTE_BOSON_MAX_NUM_PRECOMPUTED_SQRT
// Maximum allowed size of `sqr_roots` in monomial_action<boson>
#define LIBCOMMUTE_BOSON_MAX_NUM_PRECOMPUTED_SQRT 128
#endif

namespace libcommute {

template <> class monomial_action<boson> {

  // Update of a single bosonic mode
  struct single_boson_update_t {

    // The bit range corresponding to one bosonic mode is selected as
    // (in_index >> shift) & mask
    int shift;
    sv_index_type mask;

    // Change of the occupation number within this bosonic mode
    std::int64_t n_change;

    // Maximal occupation number allowed within this bosonic mode
    sv_index_type n_max;

    // n_change * (1 << shift)
    std::int64_t state_change;
  };

  // List of single-boson updates
  std::vector<single_boson_update_t> updates_;

  // Precomputed square roots of occupation numbers
  std::vector<double> sqr_roots_;

  // Will action of this monomial always result in a zero?
  bool vanishing_ = false;

  inline double sqr_root(sv_index_type n) const {
    if(n < LIBCOMMUTE_BOSON_MAX_NUM_PRECOMPUTED_SQRT)
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
      if(!is_boson(*it)) throw unknown_generator<IndexTypes...>(*it);

      if(next_it == end_it || *next_it != *it) {
        elementary_space_boson<IndexTypes...> es(0, it->indices());
        if(!hs.has(es)) throw unknown_generator<IndexTypes...>(*it);

        bit_range_t const& bit_range = hs.bit_range(es);
        int shift = bit_range.first;
        int n_bits = bit_range.second - bit_range.first + 1;
        sv_index_type n_max = (~sv_index_type(0)) >> (sv_index_width - n_bits);

        sqr_roots_size = std::max(sqr_roots_size, n_max + 1);

        // Application of this monomial would inevitably give 0. Skipping.
        if(power > n_max) {
          vanishing_ = true;
          return;
        }

        bool dagger =
            dynamic_cast<generator_boson<IndexTypes...> const&>(*it).dagger();

        auto n_change = static_cast<std::int64_t>(dagger ? power : -power);
        std::int64_t state_change = n_change * (std::int64_t(1) << shift);

        updates_.emplace_back(single_boson_update_t{
            shift,
            (~sv_index_type(0)) >> (sv_index_width - n_bits),
            n_change,
            n_max,
            state_change});

        power = 1;
      } else
        ++power;
    }

    if(!updates_.empty()) {
      sqr_roots_size =
          std::min(sqr_roots_size,
                   sv_index_type(LIBCOMMUTE_BOSON_MAX_NUM_PRECOMPUTED_SQRT));
      sqr_roots_.resize(sqr_roots_size);
      for(sv_index_type n = 0; n < sqr_roots_size; ++n)
        sqr_roots_[n] = std::sqrt(double(n));
    }
  }

  template <typename ScalarType>
  inline bool act(sv_index_type& index, ScalarType& coeff) const {

    if(vanishing_) return false;

    for(std::size_t b = updates_.size(); b-- != 0;) {
      auto const& update = updates_[b];

      auto n_part =
          static_cast<std::int64_t>((index >> update.shift) & update.mask);
      std::int64_t new_n_part = n_part + update.n_change;
      if(new_n_part < 0 || new_n_part > std::int64_t(update.n_max))
        return false;

      if(update.n_change > 0) {
        for(int d = 1; d <= update.n_change; ++d)
          mul_assign(
              coeff,
              scalar_traits<ScalarType>::make_const(sqr_root(n_part + d)));
      } else {
        for(int d = 0; d <= -update.n_change - 1; ++d)
          mul_assign(
              coeff,
              scalar_traits<ScalarType>::make_const(sqr_root(n_part - d)));
      }
      index += update.state_change;
    }
    return true;
  }
};

} // namespace libcommute

#endif
