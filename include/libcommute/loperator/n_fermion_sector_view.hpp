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
#ifndef LIBCOMMUTE_LOPERATOR_N_FERMION_SECTOR_VIEW_HPP_
#define LIBCOMMUTE_LOPERATOR_N_FERMION_SECTOR_VIEW_HPP_

#include "../algebra_ids.hpp"

#include "state_vector.hpp"

namespace libcommute {

namespace detail {

// Binomial coefficient C(n, k)
inline sv_index_type binomial(unsigned int n, unsigned int k) {
  if(k > n) return 0;
  if(k > n / 2) k = n - k;
  sv_index_type C = 1;
  for(unsigned int i = 0; i < k; ++i)
    C = (C * (n - i)) / (i + 1);
  return C;
}

} // namespace detail

// Size of the fermionic sector with N particles
template <typename HSType>
inline sv_index_type n_fermion_sector_size(HSType const& hs, unsigned int N) {
  auto total_n_bits = hs.total_n_bits();
  if(total_n_bits == 0) return 0;
  unsigned int M = 0;
  if(hs.has_algebra(fermion)) {
    auto fermion_bit_range = hs.algebra_bit_range(fermion);
    M = fermion_bit_range.second - fermion_bit_range.first + 1;
  }
  return detail::binomial(M, N) * (sv_index_type(1) << (total_n_bits - M));
}

} // namespace libcommute

#endif
