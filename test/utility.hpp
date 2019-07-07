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

#ifndef LIBCOMMUTE_TEST_UTILITY_HPP_
#define LIBCOMMUTE_TEST_UTILITY_HPP_

#include "catch2/catch.hpp"

#include <sstream>
#include <string>

//
// Utility functions used in multiple tests
//

// Checks result of calling operator<<() on an object of type `T`
template<typename T> class PrintMatcher : public Catch::MatcherBase<T> {
  std::string const& ref_;

public:
  PrintMatcher(std::string const& ref) : ref_(ref) {}

  // Performs the test for this matcher
  bool match(T const& x) const override {
    std::stringstream ss;
    ss << x;
    return ss.str() == ref_;
  }

  virtual std::string describe() const override {
    return "prints " + ref_ ;
  }
};
template<typename T>
inline PrintMatcher<T> Prints(std::string const& ref) { return {ref}; }

#endif
