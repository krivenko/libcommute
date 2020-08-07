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

#include "catch2/catch.hpp"

#include <libcommute/expression/expression.hpp>

#include <type_traits>

using namespace libcommute;

//
// Mock types used to test return types of arithmetic operations
//

struct ST1 {};
struct ST2 {};
struct ST1p2 {}; // Result type of ST1{} + ST2{}
struct ST1m2 {}; // Result type of ST1{} - ST2{}
struct ST1t2 {}; // Result type of ST1{} * ST2{}
struct STum {};  // Result type of -ST1{}

ST1 operator+(ST1 const&, ST1 const&) { return ST1{}; }
ST1p2 operator+(ST1 const&, ST2 const&) { return ST1p2{}; }

ST1 operator-(ST1 const&, ST1 const&) { return ST1{}; }
ST1m2 operator-(ST1 const&, ST2 const&) { return ST1m2{}; }

ST1 operator*(ST1 const&, ST1 const&) { return ST1{}; }
ST1t2 operator*(ST1 const&, ST2 const&) { return ST1t2{}; }

STum operator-(ST1 const&) { return STum{}; }

namespace libcommute {

template<> struct scalar_traits<ST1> {
  static ST1 make_const(double x) { return ST1{}; }
  static bool is_zero(ST1 const& x) { return true; }
};

template<> struct scalar_traits<ST2> {
  static ST2 make_const(double x) { return ST2{}; }
  static bool is_zero(ST2 const& x) { return true; }
};

template<> struct scalar_traits<ST1p2> {
  static ST1p2 make_const(double x) { return ST1p2{}; }
  static bool is_zero(ST1p2 const& x) { return true; }
};

template<> struct scalar_traits<ST1m2> {
  static ST1m2 make_const(double x) { return ST1m2{}; }
  static bool is_zero(ST1m2 const& x) { return true; }
};

template<> struct scalar_traits<ST1t2> {
  static ST1t2 make_const(double x) { return ST1t2{}; }
  static bool is_zero(ST1t2 const& x) { return true; }
};

}

TEST_CASE("Result types of arithmetic operations with mixed scalar types",
          "[mixed_types]") {
  expression<ST1> expr1;
  expression<ST2> expr2;

  SECTION("Addition") {
    CHECK(std::is_same<decltype(expr1 + expr1), expression<ST1>>::value);
    CHECK(std::is_same<decltype(expr1 + expr2), expression<ST1p2>>::value);
    CHECK(std::is_same<decltype(expr1 += expr1), expression<ST1> &>::value);
    CHECK(std::is_same<decltype(expr1 += expr2), expression<ST1> &>::value);
  }

  SECTION("Addition of a constant") {
    CHECK(std::is_same<decltype(expr1 + ST1{}), expression<ST1>>::value);
    CHECK(std::is_same<decltype(ST1{} + expr1), expression<ST1>>::value);
    CHECK(std::is_same<decltype(ST1{} + expr2), expression<ST1p2>>::value);
    CHECK(std::is_same<decltype(expr1 + ST2{}), expression<ST1p2>>::value);
    CHECK(std::is_same<decltype(expr1 += ST1{}), expression<ST1> &>::value);
    CHECK(std::is_same<decltype(expr1 += ST2{}), expression<ST1> &>::value);
  }

  SECTION("Subtraction") {
    CHECK(std::is_same<decltype(expr1 - expr1), expression<ST1>>::value);
    CHECK(std::is_same<decltype(expr1 - expr2), expression<ST1m2>>::value);
    CHECK(std::is_same<decltype(expr1 -= expr1), expression<ST1> &>::value);
    CHECK(std::is_same<decltype(expr1 -= expr2), expression<ST1> &>::value);
  }

  SECTION("Subtraction of a constant") {
    CHECK(std::is_same<decltype(expr1 - ST1{}), expression<ST1>>::value);
    CHECK(std::is_same<decltype(ST1{} - expr1), expression<ST1>>::value);
    CHECK(std::is_same<decltype(ST1{} - expr2), expression<ST1m2>>::value);
    CHECK(std::is_same<decltype(expr1 - ST2{}), expression<ST1m2>>::value);
    CHECK(std::is_same<decltype(expr1 -= ST1{}), expression<ST1> &>::value);
    CHECK(std::is_same<decltype(expr1 -= ST2{}), expression<ST1> &>::value);
  }

  SECTION("Multiplication") {
    CHECK(std::is_same<decltype(expr1 * expr1), expression<ST1>>::value);
    CHECK(std::is_same<decltype(expr1 * expr2), expression<ST1t2>>::value);
    CHECK(std::is_same<decltype(expr1 *= expr1), expression<ST1> &>::value);
    CHECK(std::is_same<decltype(expr1 *= expr2), expression<ST1> &>::value);
  }

  SECTION("Multiplication by a constant") {
    CHECK(std::is_same<decltype(expr1 * ST1{}), expression<ST1>>::value);
    CHECK(std::is_same<decltype(ST1{} * expr1), expression<ST1>>::value);
    CHECK(std::is_same<decltype(ST1{} * expr2), expression<ST1t2>>::value);
    CHECK(std::is_same<decltype(expr1 * ST2{}), expression<ST1t2>>::value);
    CHECK(std::is_same<decltype(expr1 *= ST1{}), expression<ST1> &>::value);
    CHECK(std::is_same<decltype(expr1 *= ST2{}), expression<ST1> &>::value);
  }

  SECTION("Unary minus") {
    CHECK(std::is_same<decltype(-expr1), expression<STum>>::value);
  }
}
