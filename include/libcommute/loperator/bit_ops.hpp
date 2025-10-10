/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2025 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_LOPERATOR_BIT_OPS_HPP_
#define LIBCOMMUTE_LOPERATOR_BIT_OPS_HPP_

#include <cstdint>

#if (defined(__GNUC__) || defined(__clang__)) && defined(__BMI2__)
#include <immintrin.h>
#endif

namespace libcommute {

namespace detail {

//
// Advanced bit-wise operations.
// Implemented using compiler intrinsics wherever possible.
//

// 2^n
inline constexpr std::uint64_t pow2(unsigned int n) {
  return std::uint64_t(1) << n;
}

// The lowest power of 2 that is not smaller than i
inline std::uint64_t next_pow2(std::uint64_t i) {
  if(i < 2) return 1;
#if defined(__GNUC__) || defined(__clang__)
  return std::uint64_t(1) << (64 - __builtin_clzll(i - 1));
#else
  // https://graphics.stanford.edu/%7Eseander/bithacks.html#RoundUpPowerOf2
  --i;
  i |= i >> 1;
  i |= i >> 2;
  i |= i >> 4;
  i |= i >> 8;
  i |= i >> 16;
  i |= i >> 32;
  ++i;
  return i;
#endif
}

// Count set bits
inline unsigned int popcount(std::uint64_t i) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_popcountll(i);
#else
  i ^= i >> 32;
  i ^= i >> 16;
  i ^= i >> 8;
  i ^= i >> 4;
  i ^= i >> 2;
  i ^= i >> 1;
  return i;
#endif
}

// Compute parity of the number of set bits
inline static bool parity_popcount(std::uint64_t i) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_parityll(i);
#else
  i ^= i >> 32;
  i ^= i >> 16;
  i ^= i >> 8;
  i ^= i >> 4;
  i ^= i >> 2;
  i ^= i >> 1;
  return i & 0x1;
#endif
}

// Count trailing zeroes
inline unsigned int count_trailing_zeros(std::uint64_t i) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_ctzll(i);
#else
  // https://blog.stephencleary.com/2010/10/
  //   implementing-gccs-builtin-functions.html
  unsigned int n = 1;
  if((i & 0xffffffff) == 0) {
    n = n + 32;
    i = i >> 32;
  }
  if((i & 0x0000ffff) == 0) {
    n = n + 16;
    i = i >> 16;
  }
  if((i & 0x000000ff) == 0) {
    n = n + 8;
    i = i >> 8;
  }
  if((i & 0x0000000f) == 0) {
    n = n + 4;
    i = i >> 4;
  }
  if((i & 0x00000003) == 0) {
    n = n + 2;
    i = i >> 2;
  }
  return n - (i & 1);
#endif
}

// Parallel bits deposit
inline std::uint64_t deposit_bits(std::uint64_t src, std::uint64_t mask) {
#if (defined(__GNUC__) || defined(__clang__)) && defined(__BMI2__)
  return _pdep_u64(src, mask);
#else
  // https://www.chessprogramming.org/BMI2#Serial_Implementation
  std::uint64_t res = 0;
  for(std::uint64_t bb = 1; mask; bb += bb) {
    // cppcheck-suppress oppositeExpression
    if(src & bb) res |= mask & -mask;
    mask &= mask - 1;
  }
  return res;
#endif
}

// Parallel bits extract
inline std::uint64_t extract_bits(std::uint64_t val, std::uint64_t mask) {
#if (defined(__GNUC__) || defined(__clang__)) && defined(__BMI2__)
  return _pext_u64(val, mask);
#else
  // From https://www.chessprogramming.org/BMI2#Serial_Implementation_2
  std::uint64_t res = 0;
  for(std::uint64_t bb = 1; mask; bb += bb) {
    if(val & mask & -mask) res |= bb;
    mask &= mask - 1;
  }
  return res;
#endif
}

} // namespace detail

} // namespace libcommute

#endif
