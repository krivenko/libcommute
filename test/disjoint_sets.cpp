/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <catch.hpp>

#include <libcommute/loperator/disjoint_sets.hpp>

using namespace libcommute;

#include <iostream>

TEST_CASE("Disjoint sets data structure", "[disjoint_sets]") {

  disjoint_sets ds(10);

  SECTION("Basic tests") {
    CHECK(ds.size() == 10);
    CHECK(ds.n_sets() == 10);
    ds.make_set();
    CHECK(ds.size() == 11);
    CHECK(ds.n_sets() == 11);
    ds.make_sets(4);
    CHECK(ds.size() == 15);
    CHECK(ds.n_sets() == 15);

    for(unsigned int i = 0; i < ds.n_sets(); ++i)
      CHECK(ds.find_root(i) == i);

    CHECK_FALSE(ds.in_same_set(2, 3));
  }

  SECTION("set_union() and find_root()") {
    ds.set_union(2, 3);
    CHECK(ds.n_sets() == 9);
    CHECK(ds.in_same_set(2, 3));
    CHECK(ds.find_root(3) == 2);
    ds.set_union(3, 2);
    CHECK(ds.n_sets() == 9);
    CHECK(ds.in_same_set(2, 3));
    CHECK(ds.find_root(3) == 2);
  }

  SECTION("More tests") {
    CHECK(ds.make_set() == 10);
    CHECK(ds.n_sets() == 11);
    CHECK(ds.root_union(8, 7) == 8);
    CHECK(ds.set_union(5, 6) == 5);
    CHECK(ds.set_union(8, 5) == 8);
    CHECK(ds.n_sets() == 8);
    CHECK(ds.find_root(6) == 8);
    ds.set_union(2, 6);
    CHECK(ds.find_root(2) == 8);
    std::size_t root1 = ds.find_root(6);
    std::size_t root2 = ds.find_root(2);
    CHECK(ds.root_union(root1, root2) == 8);
    CHECK(ds.set_union(5, 6) == 8);
  }

  SECTION("normalize_sets()") {
    ds.set_union(4, 7);
    ds.set_union(7, 1);
    CHECK(ds.find_root(7) == 4);
    ds.compress_sets();
    ds.normalize_sets();
    CHECK(ds.find_root(7) == 1);
  }
}
