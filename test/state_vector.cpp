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

#include "catch2/catch.hpp"

#include <libcommute/qoperator/state_vector.hpp>

using namespace libcommute;

TEST_CASE("Implementation of StateVector interface", "[state_vector]") {

  SECTION("real") {
    std::vector<double> v{1, 2, 3};

    CHECK(get_size(v) == 3);
    CHECK(get_element(v, 1) == 2);
    get_element(v, 1) = 4;
    CHECK(get_element(v, 1) == 4);
    CHECK(zeros_like(v) == std::vector<double>{0, 0, 0});
    set_seros(v);
    CHECK(v == std::vector<double>{0, 0, 0});
  }

  SECTION("complex") {
    std::vector<std::complex<double>> v{1, 2, 3};

    CHECK(get_size(v) == 3);
    CHECK(get_element(v, 1) == 2);
    get_element(v, 1) = 4;
    CHECK(get_element(v, 1) == 4);
    CHECK(zeros_like(v) == std::vector<std::complex<double>>{0, 0, 0});
    set_seros(v);
    CHECK(v == std::vector<std::complex<double>>{0, 0, 0});
  }
}
