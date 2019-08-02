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
#ifndef LIBCOMMUTE_STATE_VECTOR_HPP_
#define LIBCOMMUTE_STATE_VECTOR_HPP_

#include "../expressions/scalar_traits.hpp"

#include <vector>

// qoperator is a quantum-mechanical operator that acts on a state vector.
// The following C++20 concept describes a valid state vector type.
//template<typename T>
//concept StateVector = requires(T v) {
//  { get_size(v) } -> std::convertible_to<std::uint64_t>;
//  { get_element(v, std::uint64_t{}) };
//  { zeros_like(v) } -> T;
//  { set_zeros(v) } -> T;
//};

namespace libcommute {

//
// Functions implementing StateVector interface for std::vector<T>
//

template<typename T>
std::uint64_t get_size(std::vector<T> const& v) {
  return v.size();
}

template<typename T>
T const& get_element(std::vector<T> const& v, std::uint64_t n) {
  return v[n];
}

template<typename T>
T & get_element(std::vector<T> & v, std::uint64_t n) {
  return v[n];
}

template<typename T>
std::vector<T> zeros_like(std::vector<T> const& v) {
  return std::vector<T>(v.size(), scalar_traits<T>::make_const(0));
}

template<typename T>
void set_zeros(std::vector<T> & v) {
  std::fill(v.begin(), v.end(), scalar_traits<T>::make_const(0));
}

} // namespace libcommute
