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
#ifndef LIBCOMMUTE_TEST_CHECK_ORDERING_HPP_
#define LIBCOMMUTE_TEST_CHECK_ORDERING_HPP_

#include <catch.hpp>

#include <type_traits>
#include <vector>

//
// Check that elements of `v` are pairwise distinct
//

// T is a pointer
template<typename T>
void check_equality_impl(std::vector<T> const& v, std::true_type) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] == *v[i2]) == (i1 == i2));
      CHECK((*v[i1] != *v[i2]) == (i1 != i2));
    }
  }
}

// T is not a pointer
template<typename T>
void check_equality_impl(std::vector<T> const& v, std::false_type) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((v[i1] == v[i2]) == (i1 == i2));
      CHECK((v[i1] != v[i2]) == (i1 != i2));
    }
  }
}

template<typename T>
void check_equality(std::vector<T> const& v) {
  check_equality_impl(v,
    std::integral_constant<bool, std::is_pointer<T>::value>()
  );
}

//
// Check that elements of `v` are ordered
//

template<typename T>
void check_less_greater_impl(std::vector<T> const& v, std::true_type) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] < *v[i2]) == (i1 < i2));
      CHECK((*v[i1] > *v[i2]) == (i1 > i2));
    }
  }
}

template<typename T>
void check_less_greater_impl(std::vector<T> const& v, std::false_type) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((v[i1] < v[i2]) == (i1 < i2));
      CHECK((v[i1] > v[i2]) == (i1 > i2));
    }
  }
}

template<typename T>
void check_less_greater(std::vector<T> const& v) {
  check_less_greater_impl(v,
    std::integral_constant<bool, std::is_pointer<T>::value>()
  );
}

#endif
