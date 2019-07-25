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

#include <libcommute/qoperator/qoperator.hpp>
#include <libcommute/expression/factories.hpp>

using namespace libcommute;

TEST_CASE("Quantum-mechanical operator", "[qoperator]") {
  using namespace static_indices::real;

  auto expr = c_dag(1, "up") * c(2, "dn");
  auto qop = make_qoperator(expr);

}
