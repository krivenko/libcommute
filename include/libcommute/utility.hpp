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
#ifndef LIBCOMMUTE_UTILITY_HPP_
#define LIBCOMMUTE_UTILITY_HPP_

#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <tuple>
#include <vector>

//
// Various utility functions and metafunctions
//

namespace libcommute {

// std::make_unique<T>() from C++14
#if __cplusplus < 201402L
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
#define LIBCOMMUTE_NO_STD_MAKE_UNIQUE
#endif

// std::remove_cvref from C++20
template<typename T> struct remove_cvref {
  using type = typename std::remove_cv<
    typename std::remove_reference<T>::type
  >::type;
};
template<typename T> using remove_cvref_t = typename remove_cvref<T>::type;

// Identity metafunction for all T except for the C-string types
// that is mapped to std::string
template<typename T> struct c_str_to_string { using type = T; };
template<size_t N> struct c_str_to_string<const char (&)[N]> {
  using type = std::string;
};
template<> struct c_str_to_string<const char *> {
  using type = std::string;
};
template<typename T>
using c_str_to_string_t = typename c_str_to_string<T>::type;

namespace detail {

template<size_t N> struct print_tuple_impl {
  template<typename... T>
  static void apply(std::ostream & os, std::tuple<T...> const& t) {
    os << std::get<sizeof...(T) - 1 - N>(t) << ",";
    print_tuple_impl<N - 1>::apply(os, t);
  }
};
template<> struct print_tuple_impl<0> {
  template<typename... T>
  static void apply(std::ostream & os, std::tuple<T...> const& t) {
    os << std::get<sizeof...(T) - 1>(t);
  }
};

} // namespace detail

// Print an std::tuple to an output stream as a comma-separated list
template<typename... T>
void print_tuple(std::ostream & os, std::tuple<T...> const& t) {
  detail::print_tuple_impl<sizeof...(T) - 1>::apply(os, t);
}

// Linear function of objects
template<typename T> struct linear_function {
  // Constant term
  double const_term;
  // Generators and their respective coefficients
  std::vector<std::pair<T, double>> terms;
};

} // namespace libcommute

#endif
