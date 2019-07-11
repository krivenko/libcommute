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

#ifndef LIBCOMMUTE_TEST_MY_COMPLEX_HPP_
#define LIBCOMMUTE_TEST_MY_COMPLEX_HPP_

#include <libcommute/expression/scalar_traits.hpp>

//
// Non-default-constructible mock scalar type
//
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
  static bool is_zero(my_complex const& x) { return x.re == 0 && x.im == 0; }
  // Unitary value
  static constexpr my_complex one() { return {1, 0}; }
  // Real part of x
  static my_complex real(my_complex const& x) { return {x.re, 0}; }
  // Imaginary part of x
  static my_complex imag(my_complex const& x) { return {x.im, 0}; }
  // Complex conjugate of x
  static my_complex conj(my_complex const& x) { return {x.re, -x.im}; }
};

} // namespace libcommute

#endif
