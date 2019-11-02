/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include "catch2/catch.hpp"

#include <libcommute/qoperator/state_vector.hpp>

#include <type_traits>

using namespace libcommute;

TEST_CASE("Implementation of StateVector interface", "[state_vector]") {

  SECTION("real") {
    std::vector<double> v{1, 2, 3};

    CHECK(std::is_same<element_type_t<decltype(v)>, double>::value);
    CHECK(get_size(v) == 3);
    CHECK(get_element(v, 1) == 2);
    update_add_element(v, 1, 4);
    CHECK(get_element(v, 1) == 6);
    CHECK(zeros_like(v) == std::vector<double>{0, 0, 0});
    set_zeros(v);
    CHECK(v == std::vector<double>{0, 0, 0});

    SECTION("foreach()") {
      std::vector<double> v{1, 2, 0, 4};
      foreach(v, [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("complex") {
    std::vector<std::complex<double>> v{1, 2, 3};

    CHECK(std::is_same<element_type_t<decltype(v)>,
          std::complex<double>>::value);
    CHECK(get_size(v) == 3);
    CHECK(get_element(v, 1) == 2.0);
    update_add_element(v, 1, 4);
    CHECK(get_element(v, 1) == 6.0);
    CHECK(zeros_like(v) == std::vector<std::complex<double>>{0, 0, 0});
    set_zeros(v);
    CHECK(v == std::vector<std::complex<double>>{0, 0, 0});

    SECTION("foreach()") {
      std::vector<std::complex<double>> v{1, 2, 0, 4};
      foreach(v, [](int i, std::complex<double> a) {
        CHECK(double(i + 1) == a);
      });
    }
  }
}
