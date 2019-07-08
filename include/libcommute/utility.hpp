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

//
// Various utility functions and metafunctions
//

namespace libcommute {

// std::make_unique<T>() from C++14
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// Identity metafunction for all T except for the C-string type
// that is mapped to std::string
template<typename T> struct c_str_to_string_t { using type = T; };
template<size_t N> struct c_str_to_string_t<const char (&)[N]> {
  using type = std::string;
};

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

} // namespace libcommute

#endif
