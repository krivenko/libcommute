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

#ifndef LIBCOMMUTE_TEST_INT_COMPLEX_HPP_
#define LIBCOMMUTE_TEST_INT_COMPLEX_HPP_

#include <libcommute/expression/scalar_traits.hpp>

#include <iostream>

//
// Non-default-constructible mock scalar type
//
struct int_complex {
  int re; int im;

  int_complex() = delete;

  int_complex(int re, int im) : re(re), im(im) {}
  int_complex(int re) : re(re), im(0) {}

  friend bool operator==(int_complex const& c1, int_complex const& c2) {
    return c1.re == c2.re && c1.im == c2.im;
  }

  // Arithmetics
  int_complex operator-() const {return {-re, -im};}

  int_complex & operator+=(int_complex c) {
    re += c.re;
    im += c.im;
    return *this;
  }
  int_complex & operator-=(int_complex c) {
    re -= c.re;
    im -= c.im;
    return *this;
  }
  int_complex & operator*=(int_complex c2) {
    auto c1 = *this;
    re = c1.re*c2.re - c1.im*c2.im;
    im = c1.re*c2.im + c1.im*c2.re;
    return *this;
  }
  int_complex & operator+=(int n) { return *this += {n, 0}; }
  int_complex & operator-=(int n) { return *this -= {n, 0}; }
  int_complex & operator*=(int n) { return *this *= {n, 0}; }

  friend int_complex operator+(int_complex c, int n) {return {c.re + n, c.im};}
  friend int_complex operator+(int n, int_complex c) {return c + n;}
  friend int_complex operator+(int_complex c1, int_complex c2) {
    return {c1.re + c2.re, c1.im + c2.im};
  }

  friend int_complex operator-(int_complex c, int n) {return {c.re - n, c.im};}
  friend int_complex operator-(int n, int_complex c) {return {-c.re+n, -c.im};}
  friend int_complex operator-(int_complex c1, int_complex c2) {
    return {c1.re - c2.re, c1.im - c2.im};
  }

  friend int_complex operator*(int_complex c, int n) {return {c.re*n, c.im*n};}
  friend int_complex operator*(int n, int_complex c) {return c*n;}
  friend int_complex operator*(int_complex c1, int_complex c2) {
    return {c1.re*c2.re - c1.im*c2.im, c1.re*c2.im + c1.im*c2.re};
  }
};

std::ostream & operator<<(std::ostream & os, int_complex const& c) {
  return os << "{" << c.re << "," << c.im << "}";
}

namespace libcommute {

// Specialize scalar_traits for int_complex
template<> struct scalar_traits<int_complex> {
  // Zero value test
  static bool is_zero(int_complex const& x) { return x.re == 0 && x.im == 0; }
  // Zero value
  static int_complex zero() { return {0, 0}; }
  // Unitary value
  static int_complex one() { return {1, 0}; }
  // Real part of x
  static int_complex real(int_complex const& x) { return {x.re, 0}; }
  // Imaginary part of x
  static int_complex imag(int_complex const& x) { return {x.im, 0}; }
  // Complex conjugate of x
  static int_complex conj(int_complex const& x) { return {x.re, -x.im}; }
};

} // namespace libcommute

#endif