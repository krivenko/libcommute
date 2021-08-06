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

#include <catch.hpp>

#include "my_complex.hpp"

#include <libcommute/scalar_traits.hpp>

#include <complex>
#include <type_traits>

using namespace libcommute;

TEST_CASE("Traits of scalar types", "[scalar_traits]") {

  SECTION("is_complex") {
    CHECK_FALSE(is_complex<double>::value);
    CHECK(is_complex<std::complex<double>>::value);
  }

  SECTION("int") {
    CHECK(scalar_traits<int>::is_zero(0));
    CHECK_FALSE(scalar_traits<int>::is_zero(4));
    CHECK(scalar_traits<int>::make_const(0) == 0);
    CHECK(scalar_traits<int>::make_const(1) == 1);
    CHECK(scalar_traits<int>::conj(4) == 4);
  }

  SECTION("long") {
    CHECK(scalar_traits<long>::is_zero(0));
    CHECK_FALSE(scalar_traits<long>::is_zero(4));
    CHECK(scalar_traits<long>::make_const(0) == 0);
    CHECK(scalar_traits<long>::make_const(1) == 1);
    CHECK(scalar_traits<long>::conj(4) == 4);
  }

  SECTION("float") {
    CHECK(scalar_traits<float>::is_zero(.0));
    CHECK(scalar_traits<float>::is_zero(1e-50));
    CHECK_FALSE(scalar_traits<float>::is_zero(4.0));
    CHECK(scalar_traits<float>::make_const(0) == .0);
    CHECK(scalar_traits<float>::make_const(1) == 1.0);
    CHECK(scalar_traits<float>::conj(4.0) == 4.0);
  }

  SECTION("double") {
    CHECK(scalar_traits<double>::is_zero(.0));
    CHECK(scalar_traits<double>::is_zero(1e-50));
    CHECK_FALSE(scalar_traits<double>::is_zero(4.0));
    CHECK(scalar_traits<double>::make_const(0) == .0);
    CHECK(scalar_traits<double>::make_const(1) == 1.0);
    CHECK(scalar_traits<double>::conj(4.0) == 4.0);
  }

  SECTION("std::complex<float>") {
    using cmplx = std::complex<float>;

    CHECK(scalar_traits<cmplx>::is_zero(.0));
    CHECK(scalar_traits<cmplx>::is_zero(cmplx(1e-50, 1e-50)));
    CHECK_FALSE(scalar_traits<cmplx>::is_zero(cmplx(.0, 4.0)));
    CHECK_FALSE(scalar_traits<cmplx>::is_zero(cmplx(4.0, .0)));
    CHECK_FALSE(scalar_traits<cmplx>::is_zero(cmplx(4.0, 4.0)));
    CHECK(scalar_traits<cmplx>::make_const(0) == cmplx(0));
    CHECK(scalar_traits<cmplx>::make_const(1) == cmplx(1.0));
    CHECK(scalar_traits<cmplx>::conj(cmplx(1.0, 2.0)) == cmplx(1.0,-2.0));
  }

  SECTION("std::complex<double>") {
    using cmplx = std::complex<double>;

    CHECK(scalar_traits<cmplx>::is_zero(.0));
    CHECK(scalar_traits<cmplx>::is_zero(cmplx(1e-50, 1e-50)));
    CHECK_FALSE(scalar_traits<cmplx>::is_zero(cmplx(.0, 4.0)));
    CHECK_FALSE(scalar_traits<cmplx>::is_zero(cmplx(4.0, .0)));
    CHECK_FALSE(scalar_traits<cmplx>::is_zero(cmplx(4.0, 4.0)));
    CHECK(scalar_traits<cmplx>::make_const(0) == cmplx(0));
    CHECK(scalar_traits<cmplx>::make_const(1) == cmplx(1.0));
    CHECK(scalar_traits<cmplx>::conj(cmplx(1.0, 2.0)) == cmplx(1.0,-2.0));
  }

  SECTION("my_complex") {
    CHECK(scalar_traits<my_complex>::is_zero(my_complex{0, 0}));
    CHECK_FALSE(scalar_traits<my_complex>::is_zero(my_complex{0, 4}));
    CHECK_FALSE(scalar_traits<my_complex>::is_zero(my_complex{4, 0}));
    CHECK_FALSE(scalar_traits<my_complex>::is_zero(my_complex{4, 4}));
    CHECK(scalar_traits<my_complex>::make_const(0) == my_complex{0, 0});
    CHECK(scalar_traits<my_complex>::make_const(1) == my_complex{1, 0});
    CHECK(scalar_traits<my_complex>::conj(my_complex{1, 2}) ==
          my_complex{1, -2});
  }
}

TEST_CASE("Result types of arithmetic operations", "[arithmetic_result_type]") {
  using cmplx = std::complex<double>;

  SECTION("minus") {
    CHECK(std::is_same<minus_type<long>, long>::value);
    CHECK(std::is_same<minus_type<double>, double>::value);
    CHECK(std::is_same<minus_type<cmplx>, cmplx>::value);
  }

  SECTION("sum") {
    CHECK(std::is_same<sum_type<long, long>, long>::value);
    CHECK(std::is_same<sum_type<long, double>, double>::value);
    CHECK(std::is_same<sum_type<double, long>, double>::value);
    CHECK(std::is_same<sum_type<double, double>, double>::value);

    CHECK(std::is_same<sum_type<double, cmplx>, cmplx>::value);
    CHECK(std::is_same<sum_type<cmplx, double>, cmplx>::value);
    CHECK(std::is_same<sum_type<cmplx, cmplx>, cmplx>::value);
  }

  SECTION("difference") {
    CHECK(std::is_same<diff_type<long, long>, long>::value);
    CHECK(std::is_same<diff_type<long, double>, double>::value);
    CHECK(std::is_same<diff_type<double, long>, double>::value);
    CHECK(std::is_same<diff_type<double, double>, double>::value);

    CHECK(std::is_same<diff_type<double, cmplx>, cmplx>::value);
    CHECK(std::is_same<diff_type<cmplx, double>, cmplx>::value);
    CHECK(std::is_same<diff_type<cmplx, cmplx>, cmplx>::value);
  }

  SECTION("multiplication") {
    CHECK(std::is_same<mul_type<long, long>, long>::value);
    CHECK(std::is_same<mul_type<long, double>, double>::value);
    CHECK(std::is_same<mul_type<double, long>, double>::value);
    CHECK(std::is_same<mul_type<double, double>, double>::value);

    CHECK(std::is_same<mul_type<double, cmplx>, cmplx>::value);
    CHECK(std::is_same<mul_type<cmplx, double>, cmplx>::value);
    CHECK(std::is_same<mul_type<cmplx, cmplx>, cmplx>::value);
  }

  SECTION("my_complex") {
    CHECK(std::is_same<minus_type<my_complex>, my_complex>::value);

    CHECK(std::is_same<sum_type<my_complex, int>, my_complex>::value);
    CHECK(std::is_same<sum_type<int, my_complex>, my_complex>::value);
    CHECK(std::is_same<sum_type<my_complex, my_complex>,my_complex>::value);

    CHECK(std::is_same<diff_type<my_complex, int>, my_complex>::value);
    CHECK(std::is_same<diff_type<int, my_complex>, my_complex>::value);
    CHECK(std::is_same<diff_type<my_complex, my_complex>,my_complex>::value);

    CHECK(std::is_same<mul_type<my_complex, int>, my_complex>::value);
    CHECK(std::is_same<mul_type<int, my_complex>, my_complex>::value);
    CHECK(std::is_same<mul_type<my_complex, my_complex>, my_complex>::value);
  }
}

TEST_CASE("Detect availability of compound assignments",
          "[has_compound_assignment]") {
  SECTION("has_add_assign") {
    CHECK(has_add_assign<double, double>::value);
    CHECK(has_add_assign<std::complex<double>, double>::value);
    CHECK_FALSE(has_add_assign<double, std::complex<double>>::value);
    CHECK_FALSE(has_add_assign<my_complex, my_complex>::value);
    CHECK_FALSE(has_add_assign<my_complex, double>::value);
  }

  SECTION("has_sub_assign") {
    CHECK(has_sub_assign<double, double>::value);
    CHECK(has_sub_assign<std::complex<double>, double>::value);
    CHECK_FALSE(has_sub_assign<double, std::complex<double>>::value);
    CHECK_FALSE(has_sub_assign<my_complex, my_complex>::value);
    CHECK_FALSE(has_sub_assign<my_complex, double>::value);
  }

  SECTION("has_mul_assign") {
    CHECK(has_mul_assign<double, double>::value);
    CHECK(has_mul_assign<std::complex<double>, double>::value);
    CHECK_FALSE(has_mul_assign<double, std::complex<double>>::value);
    CHECK_FALSE(has_mul_assign<my_complex, my_complex>::value);
    CHECK_FALSE(has_mul_assign<my_complex, double>::value);
  }
}

//
// Mock types used to test functions *_assign()
//

struct ST1 { int a; };
struct ST2 { int a; };

struct ST3 {
  int a;
  ST3 & operator=(ST3 const&) = default;

  ST3 & operator+=(ST1 const& x) { a += x.a; return *this; }
  ST3 operator+(ST2 const& x) { return {a + 2*x.a}; }

  ST3 & operator-=(ST1 const& x) { a -= x.a; return *this; }
  ST3 operator-(ST2 const& x) { return {a - 2*x.a}; }

  ST3 & operator*=(ST1 const& x) { a *= x.a; return *this; }
  ST3 operator*(ST2 const& x) { return {a * 2*x.a}; }
};

TEST_CASE("Functions *_assign()", "[op_assign]") {
  ST1 x = {1};
  ST2 y = {1};

  SECTION("add_assign") {
    ST3 z1 = {2}, z2 = {2};
    CHECK(add_assign(z1, x).a == 3);
    CHECK(z1.a == 3);
    CHECK(add_assign(z2, y).a == 4);
    CHECK(z2.a == 4);
  }

  SECTION("sub_assign") {
    ST3 z1 = {2}, z2 = {2};
    CHECK(sub_assign(z1, x).a == 1);
    CHECK(z1.a == 1);
    CHECK(sub_assign(z2, y).a == 0);
    CHECK(z2.a == 0);
  }

  SECTION("mul_assign") {
    ST3 z1 = {2}, z2 = {2};
    CHECK(mul_assign(z1, x).a == 2);
    CHECK(z1.a == 2);
    CHECK(mul_assign(z2, y).a == 4);
    CHECK(z2.a == 4);
  }
}
