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

#include "int_complex.hpp"

#include <libcommute/expression/scalar_traits.hpp>

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
    CHECK(scalar_traits<int>::real(4) == 4);
    CHECK(scalar_traits<int>::imag(4) == 0);
    CHECK(scalar_traits<int>::conj(4) == 4);
  }

  SECTION("long") {
    CHECK(scalar_traits<long>::is_zero(0));
    CHECK_FALSE(scalar_traits<long>::is_zero(4));
    CHECK(scalar_traits<long>::make_const(0) == 0);
    CHECK(scalar_traits<long>::make_const(1) == 1);
    CHECK(scalar_traits<long>::real(4) == 4);
    CHECK(scalar_traits<long>::imag(4) == 0);
    CHECK(scalar_traits<long>::conj(4) == 4);
  }

  SECTION("float") {
    CHECK(scalar_traits<float>::is_zero(.0));
    CHECK(scalar_traits<float>::is_zero(1e-50));
    CHECK_FALSE(scalar_traits<float>::is_zero(4.0));
    CHECK(scalar_traits<float>::make_const(0) == .0);
    CHECK(scalar_traits<float>::make_const(1) == 1.0);
    CHECK(scalar_traits<float>::real(4.0) == 4.0);
    CHECK(scalar_traits<float>::imag(4.0) == 0);
    CHECK(scalar_traits<float>::conj(4.0) == 4.0);
  }

  SECTION("double") {
    CHECK(scalar_traits<double>::is_zero(.0));
    CHECK(scalar_traits<double>::is_zero(1e-50));
    CHECK_FALSE(scalar_traits<double>::is_zero(4.0));
    CHECK(scalar_traits<double>::make_const(0) == .0);
    CHECK(scalar_traits<double>::make_const(1) == 1.0);
    CHECK(scalar_traits<double>::real(4.0) == 4.0);
    CHECK(scalar_traits<double>::imag(4.0) == 0);
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
    CHECK(scalar_traits<cmplx>::real(cmplx(1.0, 2.0)) == cmplx(1.0));
    CHECK(scalar_traits<cmplx>::imag(cmplx(1.0, 2.0)) == cmplx(2.0));
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
    CHECK(scalar_traits<cmplx>::real(cmplx(1.0, 2.0)) == cmplx(1.0));
    CHECK(scalar_traits<cmplx>::imag(cmplx(1.0, 2.0)) == cmplx(2.0));
    CHECK(scalar_traits<cmplx>::conj(cmplx(1.0, 2.0)) == cmplx(1.0,-2.0));
  }

  SECTION("int_complex") {
    CHECK(scalar_traits<int_complex>::is_zero(int_complex{0, 0}));
    CHECK_FALSE(scalar_traits<int_complex>::is_zero(int_complex{0, 4}));
    CHECK_FALSE(scalar_traits<int_complex>::is_zero(int_complex{4, 0}));
    CHECK_FALSE(scalar_traits<int_complex>::is_zero(int_complex{4, 4}));
    CHECK(scalar_traits<int_complex>::make_const(0) == int_complex{0, 0});
    CHECK(scalar_traits<int_complex>::make_const(1) == int_complex{1, 0});
    CHECK(scalar_traits<int_complex>::real(int_complex{1, 2}) ==
          int_complex{1, 0});
    CHECK(scalar_traits<int_complex>::imag(int_complex{1, 2}) ==
          int_complex{2, 0});
    CHECK(scalar_traits<int_complex>::conj(int_complex{1, 2}) ==
          int_complex{1, -2});
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

  SECTION("int_complex") {
    CHECK(std::is_same<minus_type<int_complex>, int_complex>::value);

    CHECK(std::is_same<sum_type<int_complex, int>, int_complex>::value);
    CHECK(std::is_same<sum_type<int, int_complex>, int_complex>::value);
    CHECK(std::is_same<sum_type<int_complex, int_complex>,int_complex>::value);

    CHECK(std::is_same<diff_type<int_complex, int>, int_complex>::value);
    CHECK(std::is_same<diff_type<int, int_complex>, int_complex>::value);
    CHECK(std::is_same<diff_type<int_complex, int_complex>,int_complex>::value);

    CHECK(std::is_same<mul_type<int_complex, int>, int_complex>::value);
    CHECK(std::is_same<mul_type<int, int_complex>, int_complex>::value);
    CHECK(std::is_same<mul_type<int_complex, int_complex>, int_complex>::value);
  }
}
