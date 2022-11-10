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

#include <libcommute/loperator/sparse_state_vector.hpp>

#include <type_traits>

using namespace libcommute;

TEST_CASE("Implementation of the StateVector concept based on a sparse storage",
          "[sparse_state_vector]") {

  sparse_state_vector<double> v(100);
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;

  CHECK(std::is_same<element_type_t<decltype(v)>, double>::value);
  CHECK(v.size() == 100);
  CHECK(v.n_nonzeros() == 3);
  foreach(v, [](int i, double a) { CHECK(i + 1 == a); });

  CHECK(get_element(v, 1) == 2);
  update_add_element(v, 1, 4);
  CHECK(get_element(v, 1) == 6);
  CHECK(v.n_nonzeros() == 3);
  update_add_element(v, 50, 10);
  CHECK(get_element(v, 50) == 10);
  CHECK(v.n_nonzeros() == 4);
  update_add_element(v, 50, -10);
  CHECK(get_element(v, 50) == 0);
  CHECK(v.n_nonzeros() == 3);

  auto v2 = zeros_like(v);
  CHECK(v2.n_nonzeros() == 0);
  set_zeros(v);
  CHECK(v.n_nonzeros() == 0);

  sparse_state_vector<double> v3(100);
  v3[0] = 4;
  v3[1] = 5;
  v3[2] = 6;
  v3[10] = 0;
  v3[11] = 7;
  REQUIRE(v3.n_nonzeros() == 5);
  v3.prune();
  CHECK(v3.n_nonzeros() == 4);
  v3.prune([](double a) { return int(a) % 2 == 0; });
  CHECK(v3.n_nonzeros() == 2);
  CHECK(get_element(v3, 1) == 5);
  CHECK(get_element(v3, 11) == 7);
}
