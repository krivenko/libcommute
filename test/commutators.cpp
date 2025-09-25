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

#include <catch.hpp>

#include "commutators.hpp"
#include "my_complex.hpp"

TEST_CASE("Commutation relations (int)", "[commutators_int]") {
  make_commutators_test_case_fermions<int>();
  make_commutators_test_case_bosons<int>();
}

TEST_CASE("Commutation relations (double)", "[commutators_double]") {
  make_commutators_test_case_fermions<double>();
  make_commutators_test_case_bosons<double>();
  make_commutators_test_case_spins<double>();
}

TEST_CASE("Commutation relations (complex)", "[commutators_complex]") {
  make_commutators_test_case_fermions<std::complex<double>>();
  make_commutators_test_case_bosons<std::complex<double>>();
  make_commutators_test_case_spins<std::complex<double>>();
}

TEST_CASE("Commutation relations (my_complex)", "[commutators_my_complex]") {
  make_commutators_test_case_fermions<my_complex>();
  make_commutators_test_case_bosons<my_complex>();
  make_commutators_test_case_spins<my_complex>();
}
