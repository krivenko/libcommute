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
#include <type_traits>
#include <vector>

//
// Various utility functions and metafunctions
//

namespace libcommute {

//
// std::make_unique<T>() from C++14
//

#if __cplusplus < 201402L
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
#define LIBCOMMUTE_NO_STD_MAKE_UNIQUE
#endif

//
// std::remove_cvref and std::remove_cvref_t from C++20
//

template<typename T> struct remove_cvref {
  using type = typename std::remove_cv<
    typename std::remove_reference<T>::type
  >::type;
};
template<typename T> using remove_cvref_t = typename remove_cvref<T>::type;

//
// std::invoke_result and std::invoke_result_t from C++17
//

#if __cplusplus < 201703L

// Adopted from https://en.cppreference.com/w/cpp/types/result_of
namespace detail {
template<typename T>
struct is_reference_wrapper : std::false_type {};
template<typename U>
struct is_reference_wrapper<std::reference_wrapper<U>> : std::true_type {};

template<typename T>
struct invoke_impl {
  template<class F, class... Args>
  static auto call(F&& f, Args&&... args)
    -> decltype(std::forward<F>(f)(std::forward<Args>(args)...));
};

template<typename B, typename MT>
struct invoke_impl<MT B::*> {
  template<typename T, typename Td = typename std::decay<T>::type,
    typename = typename std::enable_if<std::is_base_of<B, Td>::value>::type
  >
  static auto get(T&& t) -> T&&;

  template<typename T, typename Td = typename std::decay<T>::type,
    typename = typename std::enable_if<is_reference_wrapper<Td>::value>::type
  >
  static auto get(T&& t) -> decltype(t.get());

  template<typename T, typename Td = typename std::decay<T>::type,
    typename = typename std::enable_if<!std::is_base_of<B, Td>::value>::type,
    typename = typename std::enable_if<!is_reference_wrapper<Td>::value>::type
  >
  static auto get(T&& t) -> decltype(*std::forward<T>(t));

  template<typename T, typename... Args, typename MT1,
    typename = typename std::enable_if<std::is_function<MT1>::value>::type
  >
  static auto call(MT1 B::*pmf, T&& t, Args&&... args)
    -> decltype((invoke_impl::get(std::forward<T>(t)).*pmf)(
      std::forward<Args>(args)...));

  template<typename T>
  static auto call(MT B::*pmd, T&& t)
    -> decltype(invoke_impl::get(std::forward<T>(t)).*pmd);
};

template<typename F,
         typename... Args,
         typename Fd = typename std::decay<F>::type>
auto INVOKE(F&& f, Args&&... args)
  -> decltype(invoke_impl<Fd>::call(std::forward<F>(f),
                                    std::forward<Args>(args)...));

template<typename AlwaysVoid, typename, typename...> struct invoke_result {};
template<typename F, typename... Args>
struct invoke_result<decltype(void(detail::INVOKE(std::declval<F>(),
                                                  std::declval<Args>()...))),
                 F, Args...> {
  using type = decltype(detail::INVOKE(std::declval<F>(),
                                       std::declval<Args>()...));
};
} // namespace libcommute::detail

template<typename F, typename... ArgTypes>
struct invoke_result : detail::invoke_result<void, F, ArgTypes...> {};

template<typename F, typename... ArgTypes>
using invoke_result_t = typename invoke_result<F, ArgTypes...>::type;

#define LIBCOMMUTE_NO_STD_INVOKE_RESULT
#endif

//
// Identity metafunction for all T except for the C-string types
// that is mapped to std::string
//

template<typename T> struct c_str_to_string { using type = T; };
template<size_t N> struct c_str_to_string<const char (&)[N]> {
  using type = std::string;
};
template<> struct c_str_to_string<const char *> {
  using type = std::string;
};
template<typename T>
using c_str_to_string_t = typename c_str_to_string<T>::type;

//
// Print an std::tuple to an output stream as a comma-separated list
//

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

} // namespace libcommute::detail

template<typename... T>
void print_tuple(std::ostream & os, std::tuple<T...> const& t) {
  detail::print_tuple_impl<sizeof...(T) - 1>::apply(os, t);
}
void print_tuple(std::ostream & os, std::tuple<> const& t) {}

//
// Multilinear function of basis objects
//

template<typename T> struct linear_function {
  // Constant term
  double const_term;
  // Basis objects and their respective coefficients
  std::vector<std::pair<T, double>> terms;
};

} // namespace libcommute

#endif
