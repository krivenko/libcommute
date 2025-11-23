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
#ifndef LIBCOMMUTE_SCALARS_GMPXX_HPP_
#define LIBCOMMUTE_SCALARS_GMPXX_HPP_

#include "../scalar_traits.hpp"

#include <gmpxx.h>

#include <stdexcept>
#include <string>
#include <type_traits>

namespace libcommute {

//
// Explicitly set result types of arithmetic operations with GMP types
//
// Normally, expressions like 'a + b' involving mp*_class objects return special
// expression template types [1]. The following specializations ensure that the
// results of the arithmetic operations are cast back to the type of
// the arguments. This way GMP's expression templates do not proliferate into
// expression's ScalarType.
//
// [1] https://gmplib.org/manual/C_002b_002b-Interface-Internals
//

template <> struct uminus_res<mpz_class> {
  using type = mpz_class;
};
template <> struct sum_res<mpz_class, mpz_class> {
  using type = mpz_class;
};
template <> struct diff_res<mpz_class, mpz_class> {
  using type = mpz_class;
};
template <> struct mul_res<mpz_class, mpz_class> {
  using type = mpz_class;
};

template <> struct uminus_res<mpq_class> {
  using type = mpq_class;
};
template <> struct sum_res<mpq_class, mpq_class> {
  using type = mpq_class;
};
template <> struct diff_res<mpq_class, mpq_class> {
  using type = mpq_class;
};
template <> struct mul_res<mpq_class, mpq_class> {
  using type = mpq_class;
};

template <> struct uminus_res<mpf_class> {
  using type = mpf_class;
};
template <> struct sum_res<mpf_class, mpf_class> {
  using type = mpf_class;
};
template <> struct diff_res<mpf_class, mpf_class> {
  using type = mpf_class;
};
template <> struct mul_res<mpf_class, mpf_class> {
  using type = mpf_class;
};

//
// mpz_class
//
template <> struct scalar_traits<mpz_class> {

  // Zero value test
  static bool is_zero(mpz_class const& x) { return x == 0; }
  // Make a constant from a double value
  static mpz_class make_const(double x) {
    if(std::nearbyint(x) != x) {
      throw std::runtime_error("Cannot convert " + std::to_string(x) +
                               " to mpz_class");
    }
    return std::nearbyint(x);
  }
  // Make a constant from a variant number
  static mpz_class make_const(var_number const& vn) {
    if(vn.number_type == var_number::integer)
      return int(vn);
    else {
      std::stringstream ss;
      ss << vn;
      throw std::runtime_error("Cannot convert the variant number " +
                               ss.str() + " to mpz_class");
    }
  }
  // Complex conjugate of x
  static mpz_class conj(mpz_class const& x) { return x; }
};

//
// mpq_class
//
template <> struct scalar_traits<mpq_class> {

  // Zero value test
  static bool is_zero(mpq_class const& x) { return x == 0; }
  // Make a constant from a double value
  static mpq_class make_const(double x) {
    if(std::nearbyint(x) != x) {
      throw std::runtime_error("Cannot convert " + std::to_string(x) +
                               " to mpq_class");
    }
    return std::nearbyint(x);
  }
  // Make a constant from a variant number
  static mpq_class make_const(var_number const& vn) {
    if(vn.number_type == var_number::integer)
      return int(vn);
    else if(vn.number_type == var_number::rational)
      return {vn.numerator(), vn.denominator()};
    else {
      std::stringstream ss;
      ss << vn;
      throw std::runtime_error("Cannot convert the variant number " +
                               ss.str() + " to mpq_class");
    }
  }
  // Complex conjugate of x
  static mpq_class conj(mpq_class const& x) { return x; }
};

//
// mpf_class
//
template <> struct scalar_traits<mpf_class> {

  // Zero value test
  static bool is_zero(mpf_class const& x) { return x == 0; }
  // Make a constant from a double value
  static mpf_class make_const(double x) { return x; }
  // Make a constant from a variant number
  static mpf_class make_const(var_number const& vn) { return double(vn); }
  // Complex conjugate of x
  static mpf_class conj(mpf_class const& x) { return x; }
};

} // namespace libcommute

#endif
