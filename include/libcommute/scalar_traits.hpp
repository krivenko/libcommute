/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_SCALAR_TRAITS_HPP_
#define LIBCOMMUTE_SCALAR_TRAITS_HPP_

#include "metafunctions.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>
#include <utility>

#ifndef LIBCOMMUTE_FLOATING_POINT_TOL_EPS
// Tolerance for detection of zero floating point values expressed in units of
// the machine epsilon for the respective floating point type.
#define LIBCOMMUTE_FLOATING_POINT_TOL_EPS 100
#endif

namespace libcommute {

// Metafunction to detect complex types
template <typename T> struct is_complex : std::false_type {};
template <typename T> struct is_complex<std::complex<T>> : std::true_type {};

// Enable template instantiation if Trait<T>::value is true
template <template <typename> class Trait, typename T>
using with_trait = typename std::enable_if<Trait<T>::value>::type;

// Traits of types that can be used as the ScalarType parameter of `expression`.
// User-defined scalar types need to spcialize this structure.
template <typename S, typename = void> struct scalar_traits {};

//
// Integral types
//
template <typename S> struct scalar_traits<S, with_trait<std::is_integral, S>> {
  // Zero value test
  static bool is_zero(S const& x) { return x == 0; }
  // Make a constant from a double value
  static S make_const(double x) {
    assert(std::nearbyint(x) == x);
    return S(std::nearbyint(x));
  }
  // Complex conjugate of x
  static S conj(S const& x) { return x; }
};

//
// Floating point types
//
template <typename S>
struct scalar_traits<S, with_trait<std::is_floating_point, S>> {
  // Zero value test
  static bool is_zero(S const& x) {
    return std::abs(x) < LIBCOMMUTE_FLOATING_POINT_TOL_EPS *
                             std::numeric_limits<S>::epsilon();
  }
  // Make a constant from a double value
  static S make_const(double x) { return x; }
  // Complex conjugate of x
  static S conj(S const& x) { return x; }
};

//
// Complex types
//
template <typename S> struct scalar_traits<S, with_trait<is_complex, S>> {
  // Zero value test
  static bool is_zero(S const& x) {
    using real_t = typename S::value_type;
    return scalar_traits<real_t>::is_zero(x.real()) &&
           scalar_traits<real_t>::is_zero(x.imag());
  }
  // Make a constant from a double value
  static S make_const(double x) { return S(x); }
  // Complex conjugate of x
  static S conj(S const& x) { return std::conj(x); }
};

//
// Result types of arithmetic operations
//

// Type of unary minus result
template <typename S> using minus_type = decltype(-std::declval<S>());

// Type of sum of two objects with types S1 and S2
template <typename S1, typename S2>
using sum_type = decltype(std::declval<S1>() + std::declval<S2>());

// Type of difference of two objects with types S1 and S2
template <typename S1, typename S2>
using diff_type = decltype(std::declval<S1>() - std::declval<S2>());

// Type of product of two objects with types S1 and S2
template <typename S1, typename S2>
using mul_type = decltype(std::declval<S1>() * std::declval<S2>());

//
// Given a pair of types S1, S2 and a binary operation OP, define a trait
// structure that detect whether 'S1 OP S2' is a valid expression.
//
#define DEFINE_HAS_OP(NAME, OP)                                                \
  template <typename S1, typename S2, typename = void>                         \
  struct has_##NAME : std::false_type {};                                      \
  template <typename S1, typename S2>                                          \
  struct has_##NAME<                                                           \
      S1,                                                                      \
      S2,                                                                      \
      void_t<decltype(std::declval<S1&>() OP std::declval<S2>())>>             \
    : std::true_type {};

DEFINE_HAS_OP(add_assign, +=)
DEFINE_HAS_OP(sub_assign, -=)
DEFINE_HAS_OP(mul_assign, *=)
#undef DEFINE_HAS_OP

//
// Define a function NAME_assign(S1 & a, S2 const& b) that calls the compound
// assignment operator 'a COMPOUND_OP b' if it is available, and 'a = a OP b'
// otherwise.
//
#define DEFINE_OP_ASSIGN_FUNC(NAME, COMPOUND_OP, OP)                           \
  template <typename S1, typename S2>                                          \
  inline S1& NAME##_assign_impl(S1& a, S2 const& b, std::true_type) {          \
    return a COMPOUND_OP b;                                                    \
  }                                                                            \
                                                                               \
  template <typename S1, typename S2>                                          \
  inline S1& NAME##_assign_impl(S1& a, S2 const& b, std::false_type) {         \
    a = a OP b;                                                                \
    return a;                                                                  \
  }                                                                            \
                                                                               \
  template <typename S1, typename S2>                                          \
  inline S1& NAME##_assign(S1& a, S2 const& b) {                               \
    return NAME##_assign_impl(a, b, has_##NAME##_assign<S1, S2>());            \
  }

DEFINE_OP_ASSIGN_FUNC(add, +=, +)
DEFINE_OP_ASSIGN_FUNC(sub, -=, -)
DEFINE_OP_ASSIGN_FUNC(mul, *=, *)
#undef DEFINE_OP_ASSIGN_FUNC

} // namespace libcommute

#endif
