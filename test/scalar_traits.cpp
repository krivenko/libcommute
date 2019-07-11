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

#include <libcommute/expression/scalar_traits.hpp>

#include <complex>
#include <type_traits>

using namespace libcommute;

// Non-default-constructible custom scalar type
struct my_complex {
  int re; int im;

  my_complex() = delete;

  friend bool operator==(my_complex const& c1, my_complex const& c2) {
    return c1.re == c2.re && c1.im == c2.im;
  }

  // Arithmetics
  my_complex operator-() {return {-re, -im};}

  friend my_complex operator+(my_complex c, int n) {return {c.re + n, c.im};}
  friend my_complex operator+(int n, my_complex c) {return c + n;}
  friend my_complex operator+(my_complex c1, my_complex c2) {
    return {c1.re + c2.re, c1.im + c2.im};
  }

  friend my_complex operator-(my_complex c, int n) {return {c.re - n, c.im};}
  friend my_complex operator-(int n, my_complex c) {return {-c.re + n, -c.im};}
  friend my_complex operator-(my_complex c1, my_complex c2) {
    return {c1.re - c2.re, c1.im - c2.im};
  }

  friend my_complex operator*(my_complex c, int n) {return {c.re*n, c.im*n};}
  friend my_complex operator*(int n, my_complex c) {return c*n;}
  friend my_complex operator*(my_complex c1, my_complex c2) {
    return {c1.re*c2.re - c1.im*c2.im, c1.re*c2.im + c1.im*c2.re};
  }
};

namespace libcommute {

// Specialize scalar_traits for my_complex
template<> struct scalar_traits<my_complex> {
  // Zero value test
  static bool is_zero(my_complex const& x) { return x.re == 0 && x.im == 0;  }
  // Real part of x
  static my_complex real(my_complex const& x) { return {x.re, 0}; }
  // Imaginary part of x
  static my_complex imag(my_complex const& x) { return {x.im, 0}; }
  // Complex conjugate of x
  static my_complex conj(my_complex const& x) { return {x.re, -x.im}; }
};

}

TEST_CASE("Traits of scalar types", "[scalar_traits]") {

  SECTION("is_complex") {
    CHECK_FALSE(is_complex<double>::value);
    CHECK(is_complex<std::complex<double>>::value);
  }

  SECTION("int") {
    CHECK(scalar_traits<int>::is_zero(0));
    CHECK_FALSE(scalar_traits<int>::is_zero(4));
    CHECK(scalar_traits<int>::real(4) == 4);
    CHECK(scalar_traits<int>::imag(4) == 0);
    CHECK(scalar_traits<int>::conj(4) == 4);
  }

  SECTION("long") {
    CHECK(scalar_traits<long>::is_zero(0));
    CHECK_FALSE(scalar_traits<long>::is_zero(4));
    CHECK(scalar_traits<long>::real(4) == 4);
    CHECK(scalar_traits<long>::imag(4) == 0);
    CHECK(scalar_traits<long>::conj(4) == 4);
  }

  SECTION("float") {
    CHECK(scalar_traits<float>::is_zero(.0));
    CHECK(scalar_traits<float>::is_zero(1e-50));
    CHECK_FALSE(scalar_traits<float>::is_zero(4.0));
    CHECK(scalar_traits<float>::real(4.0) == 4.0);
    CHECK(scalar_traits<float>::imag(4.0) == 0);
    CHECK(scalar_traits<float>::conj(4.0) == 4.0);
  }

  SECTION("double") {
    CHECK(scalar_traits<double>::is_zero(.0));
    CHECK(scalar_traits<double>::is_zero(1e-50));
    CHECK_FALSE(scalar_traits<double>::is_zero(4.0));
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
    CHECK(scalar_traits<cmplx>::real(cmplx(1.0, 2.0)) == cmplx(1.0));
    CHECK(scalar_traits<cmplx>::imag(cmplx(1.0, 2.0)) == cmplx(2.0));
    CHECK(scalar_traits<cmplx>::conj(cmplx(1.0, 2.0)) == cmplx(1.0,-2.0));
  }

  SECTION("my_complex") {
    CHECK(scalar_traits<my_complex>::is_zero(my_complex{0, 0}));
    CHECK_FALSE(scalar_traits<my_complex>::is_zero(my_complex{0, 4}));
    CHECK_FALSE(scalar_traits<my_complex>::is_zero(my_complex{4, 0}));
    CHECK_FALSE(scalar_traits<my_complex>::is_zero(my_complex{4, 4}));
    CHECK(scalar_traits<my_complex>::real(my_complex{1, 2}) ==
          my_complex{1, 0});
    CHECK(scalar_traits<my_complex>::imag(my_complex{1, 2}) ==
          my_complex{2, 0});
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
    CHECK(std::is_same<sum_type<my_complex, my_complex>, my_complex>::value);

    CHECK(std::is_same<diff_type<my_complex, int>, my_complex>::value);
    CHECK(std::is_same<diff_type<int, my_complex>, my_complex>::value);
    CHECK(std::is_same<diff_type<my_complex, my_complex>, my_complex>::value);

    CHECK(std::is_same<mul_type<my_complex, int>, my_complex>::value);
    CHECK(std::is_same<mul_type<int, my_complex>, my_complex>::value);
    CHECK(std::is_same<mul_type<my_complex, my_complex>, my_complex>::value);
  }
}
