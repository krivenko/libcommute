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

#include <catch.hpp>

#include <libcommute/loperator/bit_ops.hpp>

using namespace libcommute::detail;

TEST_CASE("Advanced bit-wise operations", "[bit_ops]") {

  // cppcheck-suppress-begin knownConditionTrueFalse
  SECTION("2^n") {
    CHECK(pow2(0) == 1);
    CHECK(pow2(1) == 2);
    CHECK(pow2(4) == 16);
  }

  SECTION("next_pow2") {
    CHECK(next_pow2(0) == 1);
    CHECK(next_pow2(1) == 1);
    CHECK(next_pow2(2) == 2);
    CHECK(next_pow2(3) == 4);
    CHECK(next_pow2(4) == 4);
    CHECK(next_pow2(5) == 8);
    CHECK(next_pow2(6) == 8);
    CHECK(next_pow2(7) == 8);
    CHECK(next_pow2(8) == 8);
  }

  SECTION("popcount") {
    for(unsigned int i = 0; i < 63; ++i)
      CHECK(popcount(std::uint64_t(1) << i) == 1);

    std::uint64_t one = 1;
    CHECK(popcount((one << 32) + (one << 8)) == 2);
    CHECK(popcount((one << 32) + (one << 4) + (one << 2)) == 3);
    CHECK(popcount((one << 32) + (one << 8) + (one << 4) + (one << 2)) == 4);
  }

  SECTION("count_trailing_zeros") {
    for(unsigned int i = 0; i < 63; ++i)
      CHECK(count_trailing_zeros(std::uint64_t(1) << i) == i);
  }

  SECTION("deposit_bits") {
    CHECK(deposit_bits(0x3F, 0x07020408) == 0x07020408);
    CHECK(deposit_bits(0x38, 0x07020408) == 0x07000000);
    CHECK(deposit_bits(0x0F, 0x07020408) == 0x01020408);
  }

  SECTION("extract_bits") {
    CHECK(extract_bits(0x00FF00FF, 0x07020408) == 5 /* 0b101 */);
    CHECK(extract_bits(0xFF00FF00, 0x07020408) == 58 /* 0b111010 */);
    CHECK(extract_bits(0xFF00000000000000, 0x0100000000000000) == 1);
  }
  // cppcheck-suppress-end knownConditionTrueFalse
}
