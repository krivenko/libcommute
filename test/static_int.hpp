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

#ifndef LIBCOMMUTE_TEST_STATIC_INT_HPP_
#define LIBCOMMUTE_TEST_STATIC_INT_HPP_

#include <libcommute/expression/scalar_traits.hpp>

#include <iostream>

//
// Integer type with its value known at compile type
// Implements a variety of arithmetic operations
//

template<int N> struct static_int {
  static_int() = default;

  static_int<-N> operator-() const {return {};}

  static_int<N> operator*=(static_int<1>) {return *this;}
  static_int<N> operator+=(static_int<0>) {return *this;}
  static_int<N> operator-=(static_int<0>) {return *this;}
};

template<int N1, int N2>
bool operator==(static_int<N1> const&, static_int<N2> const&) {
  return N1 == N2;
}

template<int N1, int N2>
static_int<N1+N2>
operator+(static_int<N1> const&, static_int<N2> const&) {return {};}

template<int N1, int N2>
static_int<N1-N2>
operator-(static_int<N1> const&, static_int<N2> const&) {return {};}

template<int N1, int N2>
static_int<N1*N2>
operator*(static_int<N1> const&, static_int<N2> const&) {return {};}

template<int N>
std::ostream & operator<<(std::ostream & os, static_int<N>) {
  return os << "static_int<" << N << ">" << std::endl;
}

namespace libcommute {

// Specialize scalar_traits for static_int<N>
template<int N> struct scalar_traits<static_int<N>> {
  // Zero value test
  static bool is_zero(static_int<N> const& x) { return N == 0; }
  // Unitary value
  static constexpr static_int<1> one() { return {}; }
  // Real part of x
  static static_int<N> real(static_int<1> const&) { return {}; }
  // Imaginary part of x
  static static_int<0> imag(static_int<1> const&) { return {}; }
  // Complex conjugate of x
  static static_int<N> conj(static_int<1> const&) { return {}; }
};

} // namespace libcommute

#endif
