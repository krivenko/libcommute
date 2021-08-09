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

#include <catch.hpp>

#include <libcommute/utility.hpp>

#include <complex>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>

using namespace libcommute;

TEST_CASE("c_str_to_string_t<T> metafunction", "[c_str_to_string_t]") {
  CHECK(std::is_same<c_str_to_string_t<int>, int>::value);
  CHECK(std::is_same<c_str_to_string_t<double>, double>::value);
  CHECK(std::is_same<c_str_to_string_t<std::string>, std::string>::value);
  CHECK(std::is_same<c_str_to_string_t<void>, void>::value);
  CHECK(std::is_same<c_str_to_string_t<decltype("Hello, world!")>,
                     std::string>::value);

  char const* hw = "Hello, world!";
  // cppcheck-suppress constVariable
  auto& hw_ref = hw;
  auto const& hw_cref = hw;
  auto&& hw_rref = hw;
  CHECK(std::is_same<c_str_to_string_t<decltype(hw_ref)>, std::string>::value);
  CHECK(std::is_same<c_str_to_string_t<decltype(hw_cref)>, std::string>::value);
  CHECK(std::is_same<c_str_to_string_t<decltype(hw_rref)>, std::string>::value);
}

TEST_CASE("all_derived_from<Base, Types...> metafunction",
          "[all_derived_from]") {
  class my_class : public std::string {};

  CHECK_FALSE(all_derived_from<std::string>::value);
  CHECK_FALSE(all_derived_from<std::string, int>::value);
  CHECK(all_derived_from<std::string, std::string>::value);
  CHECK(all_derived_from<std::string, my_class>::value);
  CHECK_FALSE(all_derived_from<std::string, int, double>::value);
  CHECK_FALSE(all_derived_from<std::string, int, std::string>::value);
  CHECK_FALSE(all_derived_from<std::string, int, my_class>::value);
  CHECK_FALSE(all_derived_from<std::string, std::string, double>::value);
  CHECK(all_derived_from<std::string, std::string, std::string>::value);
  CHECK(all_derived_from<std::string, std::string, my_class>::value);
  CHECK_FALSE(all_derived_from<std::string, my_class, double>::value);
  CHECK(all_derived_from<std::string, my_class, std::string>::value);
  CHECK(all_derived_from<std::string, my_class, my_class>::value);
}

TEST_CASE("not_in_type_list<T, TypeList...> metafunction",
          "[not_in_type_list]") {
  CHECK(not_in_type_list<int>::value);
  CHECK_FALSE(not_in_type_list<int, int>::value);
  CHECK(not_in_type_list<int, double>::value);
  CHECK_FALSE(not_in_type_list<int, int, double>::value);
  CHECK_FALSE(not_in_type_list<int, int, int>::value);
  CHECK_FALSE(not_in_type_list<int, double, int>::value);
  CHECK(not_in_type_list<int, double, double>::value);
}

TEST_CASE("print_tuple()", "[print_tuple]") {
  std::stringstream ss;

  std::tuple<> t0;
  print_tuple(ss, t0);
  CHECK(ss.str() == "");

  ss.str().clear();
  std::tuple<int, std::string, double> t(5, "Hello, World!", 1.2);
  print_tuple(ss, t);
  CHECK(ss.str() == "5,Hello, World!,1.2");
}

TEST_CASE("linear_function<T>", "[linear_function]") {

  linear_function<std::string> f0;
  CHECK(f0.terms.empty());
  CHECK(f0.vanishing());

  linear_function<std::string> f1(3.0);
  CHECK(f1.const_term == 3.0);
  CHECK(f1.terms.empty());
  CHECK_FALSE(f1.vanishing());

  std::vector<std::pair<std::string, double>> terms = {
      std::make_pair("obj1", 2.0),
      std::make_pair("obj2", 3.0)};
  linear_function<std::string> f2(4.0, terms);
  CHECK(f2.const_term == 4.0);
  CHECK(f2.terms.size() == 2);
  CHECK(f2.terms[0].first == std::string("obj1"));
  CHECK(f2.terms[0].second == 2.0);
  CHECK(f2.terms[1].first == std::string("obj2"));
  CHECK(f2.terms[1].second == 3.0);
  CHECK_FALSE(f2.vanishing());

  linear_function<std::string> f(4.0, "obj1", 2.0, "obj2", 3.0);
  CHECK(f.const_term == 4.0);
  CHECK(f.terms.size() == 2);
  CHECK(f.terms[0].first == std::string("obj1"));
  CHECK(f.terms[0].second == 2.0);
  CHECK(f.terms[1].first == std::string("obj2"));
  CHECK(f.terms[1].second == 3.0);
  CHECK_FALSE(f.vanishing());

  f.set(5.0, "obj3", 4.0);
  CHECK(f.const_term == 5.0);
  CHECK(f.terms.size() == 1);
  CHECK(f.terms[0].first == std::string("obj3"));
  CHECK(f.terms[0].second == 4.0);
  CHECK_FALSE(f.vanishing());

  f.set(6.0);
  CHECK(f.const_term == 6.0);
  CHECK(f.terms.empty());
  CHECK_FALSE(f.vanishing());
}
