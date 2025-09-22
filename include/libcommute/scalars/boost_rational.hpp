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
#ifndef LIBCOMMUTE_SCALARS_BOOST_RATIONAL_HPP_
#define LIBCOMMUTE_SCALARS_BOOST_RATIONAL_HPP_

#include "../scalar_traits.hpp"

#include <boost/rational.hpp>

#include <stdexcept>
#include <string>

namespace libcommute {

//
// Scalar traits of boost::rational
//
template <typename IntType> struct scalar_traits<boost::rational<IntType>> {

  using rational_t = boost::rational<IntType>;

  // Zero value test
  static bool is_zero(rational_t const& x) { return x == 0; }
  // Make a constant from a double value
  static rational_t make_const(double x) {
    if(std::nearbyint(x) != x) {
      throw std::runtime_error("Cannot convert " + std::to_string(x) +
                               " to boost::rational");
    }
    return rational_t(IntType(std::nearbyint(x)));
  }
  // Make a constant from a variadic number
  static rational_t make_const(var_number const& vn) {
    if(vn.number_type == var_number::integer)
      return rational_t(int(vn));
    else if(vn.number_type == var_number::rational)
      return rational_t(vn.numerator(), vn.denominator());
    else {
      std::stringstream ss;
      ss << vn;
      throw std::runtime_error("Cannot convert the variadic number " +
                               ss.str() + " to boost::rational");
    }
  }
  // Complex conjugate of x
  static rational_t conj(rational_t const& x) { return x; }
};

} // namespace libcommute

#endif
