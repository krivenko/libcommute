/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include "catch2/catch.hpp"

#include <libcommute/metafunctions.hpp>

#include <functional>
#include <string>
#include <type_traits>
#include <utility>

using namespace libcommute;

TEST_CASE("make_unique<T>()", "[make_unique]") {
  using T = std::pair<int, std::string>;
  auto p = make_unique<T>(3, "test");

  CHECK(std::is_same<decltype(p), std::unique_ptr<T>>::value);
  CHECK(p->first == 3);
  CHECK(p->second == "test");
}

TEST_CASE("remove_cvref<T> metafunction", "[remove_cvref]") {
  CHECK(std::is_same<typename remove_cvref<int>::type, int>::value);
  CHECK(std::is_same<typename remove_cvref<const int>::type, int>::value);
  CHECK(std::is_same<typename remove_cvref<int&>::type, int>::value);
  CHECK(std::is_same<typename remove_cvref<int&&>::type, int>::value);
  CHECK(std::is_same<typename remove_cvref<int const&>::type, int>::value);
  CHECK(std::is_same<typename remove_cvref<int volatile>::type, int>::value);
  CHECK(std::is_same<typename remove_cvref<int volatile&>::type, int>::value);
  CHECK(std::is_same<typename remove_cvref<int volatile&&>::type, int>::value);
  CHECK(std::is_same<remove_cvref_t<int>, int>::value);
  CHECK(std::is_same<remove_cvref_t<const int>, int>::value);
  CHECK(std::is_same<remove_cvref_t<const int>, int>::value);
  CHECK(std::is_same<remove_cvref_t<int&&>, int>::value);
  CHECK(std::is_same<remove_cvref_t<int const&>, int>::value);
  CHECK(std::is_same<remove_cvref_t<int volatile>, int>::value);
  CHECK(std::is_same<remove_cvref_t<int volatile&>, int>::value);
  CHECK(std::is_same<remove_cvref_t<int volatile&&>, int>::value);
}

// Free function
double f(int x) { return 0; }

// Lambda-functions
auto lambda1 = [](int x) -> double { return 0; };
auto lambda2 = [](int x, std::string const& s) -> std::string { return ""; };

// Functor
struct MyFunctor {
  double operator()(int x) { return 0; }
  std::string operator()(int x, std::string const& s) { return ""; }
};

// Static method
class MyClass {
public:
  static double method(int x) { return 0; }
};

// std::function
auto std_f1 = std::function<double(int)>(lambda1);
auto std_f2 = std::function<std::string(int x, std::string const& s)>(lambda2);

template<typename RefT, typename T, typename... Args>
void check_invoke_result() {
  CHECK(std::is_same<typename invoke_result<T, Args...>::type, RefT>::value);
  CHECK(std::is_same<invoke_result_t<T, Args...>, RefT>::value);
}

TEST_CASE("invoke_result<F> metafunction", "[invoke_result]") {
  check_invoke_result<double, decltype(f), int>();

  check_invoke_result<double, decltype(lambda1), int>();
  check_invoke_result<std::string, decltype(lambda2), int, std::string>();

  check_invoke_result<double, MyFunctor, int>();
  check_invoke_result<std::string, MyFunctor, int, std::string>();

  check_invoke_result<double, decltype(MyClass::method), int>();

  check_invoke_result<double, decltype(std_f1), int>();
  check_invoke_result<std::string, decltype(std_f2), int, std::string>();
}
