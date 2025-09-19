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
#ifndef LIBCOMMUTE_TEST_INT_COMPLEX_HPP_
#define LIBCOMMUTE_TEST_INT_COMPLEX_HPP_

#include <libcommute/scalar_traits.hpp>

#include <iostream>

//
// Non-default-constructible mock scalar type
//
struct my_complex {
  float re;
  float im;

  my_complex() = delete;

  my_complex(float re, float im) : re(re), im(im) {}
  // cppcheck-suppress noExplicitConstructor
  my_complex(float re) : re(re), im(0) {}

  friend bool operator==(my_complex const& c1, my_complex const& c2) {
    return c1.re == c2.re && c1.im == c2.im;
  }

  // Arithmetics
  my_complex operator-() const { return {-re, -im}; }

  friend my_complex operator+(my_complex c, float x) {
    return {c.re + x, c.im};
  }
  friend my_complex operator+(float x, my_complex c) { return c + x; }
  friend my_complex operator+(my_complex c1, my_complex c2) {
    return {c1.re + c2.re, c1.im + c2.im};
  }

  friend my_complex operator-(my_complex c, float x) {
    return {c.re - x, c.im};
  }
  friend my_complex operator-(float x, my_complex c) {
    return {-c.re + x, -c.im};
  }
  friend my_complex operator-(my_complex c1, my_complex c2) {
    return {c1.re - c2.re, c1.im - c2.im};
  }

  friend my_complex operator*(my_complex c, float x) {
    return {c.re * x, c.im * x};
  }
  friend my_complex operator*(float x, my_complex c) { return c * x; }
  friend my_complex operator*(my_complex c1, my_complex c2) {
    return {c1.re * c2.re - c1.im * c2.im, c1.re * c2.im + c1.im * c2.re};
  }

  float operator()(int m) const { return re * static_cast<float>(m); }
  my_complex operator()(int m1, int m2) const {
    return {static_cast<float>(m1 + m2) * re, static_cast<float>(m1 - m2) * im};
  }
};

inline std::ostream& operator<<(std::ostream& os, my_complex const& c) {
  return os << "{" << (c.re == -0 ? 0 : c.re) << "," << (c.im == -0 ? 0 : c.im)
            << "}";
}

namespace libcommute {

// Specialize scalar_traits for my_complex
template <> struct scalar_traits<my_complex> {
  // Zero value test
  static bool is_zero(my_complex const& x) { return x.re == 0 && x.im == 0; }
  // Make a constant from a double value
  static my_complex make_const(var_number const& x) {
    return {static_cast<float>(double(x)), 0};
  }
  // Real part of x
  static my_complex real(my_complex const& x) { return {x.re, 0}; }
  // Imaginary part of x
  static my_complex imag(my_complex const& x) { return {x.im, 0}; }
  // Complex conjugate of x
  static my_complex conj(my_complex const& x) { return {x.re, -x.im}; }
};

} // namespace libcommute

#endif
