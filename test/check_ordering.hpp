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
#ifndef LIBCOMMUTE_TEST_CHECK_ORDERING_HPP_
#define LIBCOMMUTE_TEST_CHECK_ORDERING_HPP_

#include <catch.hpp>

#include <memory>
#include <vector>

//
// Check that elements of `v` are pairwise distinct
//

template <typename T>
void check_equality(std::vector<std::shared_ptr<T>> const& v) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] == *v[i2]) == (i1 == i2));
      CHECK((*v[i1] != *v[i2]) == (i1 != i2));
    }
  }
}

template <typename T> void check_equality(std::vector<T*> const& v) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] == *v[i2]) == (i1 == i2));
      CHECK((*v[i1] != *v[i2]) == (i1 != i2));
    }
  }
}

template <typename T> void check_equality(std::vector<T> const& v) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((v[i1] == v[i2]) == (i1 == i2));
      CHECK((v[i1] != v[i2]) == (i1 != i2));
    }
  }
}

//
// Check that elements of `v` are ordered
//

template <typename T>
void check_less_greater(std::vector<std::shared_ptr<T>> const& v) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] < *v[i2]) == (i1 < i2));
      CHECK((*v[i1] > *v[i2]) == (i1 > i2));
    }
  }
}

template <typename T> void check_less_greater(std::vector<T*> const& v) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((*v[i1] < *v[i2]) == (i1 < i2));
      CHECK((*v[i1] > *v[i2]) == (i1 > i2));
    }
  }
}

template <typename T> void check_less_greater(std::vector<T> const& v) {
  for(std::size_t i1 = 0; i1 < v.size(); ++i1) {
    for(std::size_t i2 = 0; i2 < v.size(); ++i2) {
      CHECK((v[i1] < v[i2]) == (i1 < i2));
      CHECK((v[i1] > v[i2]) == (i1 > i2));
    }
  }
}

#endif
