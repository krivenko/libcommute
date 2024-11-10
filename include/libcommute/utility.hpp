/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2024 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_UTILITY_HPP_
#define LIBCOMMUTE_UTILITY_HPP_

#include "metafunctions.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

//
// Various utility functions and metafunctions
//

namespace libcommute {

//
// c_str_to_string_t<T> maps T to itself, unless T is a C-string that is
// mapped to std::string.
//

template <typename T> struct c_str_to_string {
  using type = remove_cvref_t<T>;
};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
template <std::size_t N> struct c_str_to_string<char const (&)[N]> {
  using type = std::string;
};
template <> struct c_str_to_string<char const*> {
  using type = std::string;
};
template <> struct c_str_to_string<char const*&> {
  using type = std::string;
};
template <> struct c_str_to_string<char const* const&> {
  using type = std::string;
};
template <> struct c_str_to_string<char const*&&> {
  using type = std::string;
};
template <typename T>
using c_str_to_string_t = typename c_str_to_string<T>::type;

//
// Check that all Types... are derived from Base
//

namespace detail {

template <typename Base, typename T, typename... Tail>
struct all_derived_from_impl {
  // NOLINTNEXTLINE(modernize-type-traits)
  using T_ = typename std::remove_reference<T>::type;
  // NOLINTNEXTLINE(modernize-type-traits)
  static constexpr bool value = std::is_base_of<Base, T_>::value &&
                                all_derived_from_impl<Base, Tail...>::value;
};
template <typename Base, typename T> struct all_derived_from_impl<Base, T> {
  // NOLINTNEXTLINE(modernize-type-traits)
  using T_ = typename std::remove_reference<T>::type;
  // NOLINTNEXTLINE(modernize-type-traits)
  static constexpr bool value = std::is_base_of<Base, T_>::value;
};

} // namespace detail

template <typename Base, typename... Types>
struct all_derived_from : detail::all_derived_from_impl<Base, Types...> {};
template <typename Base> struct all_derived_from<Base> : std::false_type {};

//
// Check that T is not in TypeList
//

namespace detail {

template <typename T, typename TypeListHead, typename... TypeListTail>
struct not_in_type_list_impl {
  static constexpr bool value =
      // NOLINTNEXTLINE(modernize-type-traits)
      (!std::is_same<T, TypeListHead>::value) &&
      not_in_type_list_impl<T, TypeListTail...>::value;
};
template <typename T, typename TypeListHead>
struct not_in_type_list_impl<T, TypeListHead> {
  // NOLINTNEXTLINE(modernize-type-traits)
  static constexpr bool value = !std::is_same<T, TypeListHead>::value;
};

} // namespace detail

template <typename T, typename... TypeList>
struct not_in_type_list : detail::not_in_type_list_impl<T, TypeList...> {};
template <typename T> struct not_in_type_list<T> : std::true_type {};

//
// Print an std::tuple to an output stream as a comma-separated list
//

namespace detail {

template <std::size_t N> struct print_tuple_impl {
  template <typename... T>
  static void apply(std::ostream& os, std::tuple<T...> const& t) {
    os << std::get<sizeof...(T) - 1 - N>(t) << ",";
    print_tuple_impl<N - 1>::apply(os, t);
  }
};
template <> struct print_tuple_impl<0> {
  template <typename... T>
  static void apply(std::ostream& os, std::tuple<T...> const& t) {
    os << std::get<sizeof...(T) - 1>(t);
  }
};

} // namespace detail

template <typename... T>
// cppcheck-suppress constParameterReference
inline void print_tuple(std::ostream& os, std::tuple<T...> const& t) {
  detail::print_tuple_impl<sizeof...(T) - 1>::apply(os, t);
}
// cppcheck-suppress constParameterReference
inline void print_tuple(std::ostream& os, std::tuple<> const& t) {}

//
// Trivial copyable and non-copyable objects
//

struct copyable {};
struct noncopyable {
  noncopyable() = default;
  ~noncopyable() = default;
  noncopyable(noncopyable const&) = delete;
  noncopyable(noncopyable&&) noexcept = default;
  noncopyable& operator=(noncopyable const&) = delete;
  noncopyable& operator=(noncopyable&&) noexcept = default;
};

//
// Linear function of basis objects
//

template <typename T>
// NOLINTNEXTLINE(modernize-type-traits)
struct linear_function : std::conditional<std::is_copy_constructible<T>::value,
                                          copyable,
                                          noncopyable>::type {

  using basis_type = T;

  linear_function() = default;
  explicit linear_function(double const_term) : const_term(const_term) {}
  template <typename... Args>
  linear_function(double const_term, Args&&... args) : const_term(const_term) {
    static_assert(sizeof...(Args) % 2 == 0,
                  "This constructor requires an odd number of arguments");
    construct_impl(std::forward<Args>(args)...);
  }
  linear_function(double const_term, std::vector<std::pair<T, double>> terms)
    : const_term(const_term), terms(std::move(terms)) {}

  linear_function(linear_function const&) = default;
  // NOLINTNEXTLINE(performance-noexcept-move-constructor)
  linear_function(linear_function&&) = default;
  linear_function& operator=(linear_function const&) = default;
  // NOLINTNEXTLINE(performance-noexcept-move-constructor)
  linear_function& operator=(linear_function&&) = default;
  ~linear_function() = default;

  // Reset contents
  template <typename... Args> void set(double const_term, Args&&... args) {
    this->const_term = const_term;
    terms.clear();
    construct_impl(std::forward<Args>(args)...);
  }
  void set(double const_term) {
    this->const_term = const_term;
    terms.clear();
  }

  // Is this linear function identically zero?
  bool vanishing() const { return const_term == 0 && terms.empty(); }

  // Constant term
  double const_term = 0;
  // Basis objects and their respective coefficients
  std::vector<std::pair<T, double>> terms;

private:
  template <typename T1, typename T2, typename... Tail>
  void construct_impl(T1&& arg1, T2&& arg2, Tail&&... tail) {
    terms.emplace_back(std::forward<T1>(arg1), std::forward<T2>(arg2));
    construct_impl(std::forward<Tail>(tail)...);
  }
  template <typename T1, typename T2>
  void construct_impl(T1&& arg1, T2&& arg2) {
    terms.emplace_back(std::forward<T1>(arg1), std::forward<T2>(arg2));
  }
};

} // namespace libcommute

#endif
