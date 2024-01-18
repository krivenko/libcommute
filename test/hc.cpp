/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <catch.hpp>

#include "my_complex.hpp"

#include <libcommute/expression/factories.hpp>
#include <libcommute/expression/hc.hpp>

using namespace libcommute;
using namespace static_indices;

TEST_CASE("Hermitian conjugate object (double)", "[hc_double]") {
  auto expr = 2 * c_dag("up", 1) * c("up", 2);
  CHECK((expr + hc) == (expr + conj(expr)));
  CHECK((expr - hc) == (expr - conj(expr)));
}

TEST_CASE("Hermitian conjugate object (complex)", "[hc_complex]") {
  auto expr = std::complex<double>(2, 2) * c_dag("up", 1) * c("up", 2);
  CHECK((expr + hc) == (expr + conj(expr)));
  CHECK((expr - hc) == (expr - conj(expr)));
}

TEST_CASE("Hermitian conjugate object (my_complex)", "[hc_my_complex]") {
  auto expr =
      my_complex(2, 2) * c_dag<my_complex>("up", 1) * c<my_complex>("up", 2);
  CHECK((expr + hc) == (expr + conj(expr)));
  CHECK((expr - hc) == (expr - conj(expr)));
}
