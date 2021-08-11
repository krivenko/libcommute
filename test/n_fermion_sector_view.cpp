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

#include <catch.hpp>

#include <libcommute/loperator/elementary_space_boson.hpp>
#include <libcommute/loperator/elementary_space_fermion.hpp>
#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/n_fermion_sector_view.hpp>

using namespace libcommute;

TEST_CASE("Binomial coefficient", "[binomial]") {
  using detail::binomial;

  REQUIRE(binomial(0, 0) == 1);
  REQUIRE(binomial(0, 1) == 0);
  REQUIRE(binomial(0, 2) == 0);
  REQUIRE(binomial(1, 0) == 1);
  REQUIRE(binomial(1, 1) == 1);
  REQUIRE(binomial(1, 2) == 0);

  REQUIRE(binomial(10, 4) == 210);
  REQUIRE(binomial(10, 5) == 252);
  REQUIRE(binomial(10, 6) == 210);
}

TEST_CASE("View of a state vector projected on a single N-fermion sector",
          "[n_fermion_sector_view]") {

  using namespace static_indices;

  using hs_type = hilbert_space<int>;

  SECTION("n_fermion_sector_size") {
    hs_type hs;

    // Empty Hilbert space
    CHECK(n_fermion_sector_size(hs, 0) == 0);
    CHECK(n_fermion_sector_size(hs, 1) == 0);

    // Purely fermionic Hilbert spaces
    hs.add(make_space_fermion(0));
    CHECK(n_fermion_sector_size(hs, 0) == 1);
    CHECK(n_fermion_sector_size(hs, 1) == 1);
    hs.add(make_space_fermion(1));
    CHECK(n_fermion_sector_size(hs, 0) == 1);
    CHECK(n_fermion_sector_size(hs, 1) == 2);
    CHECK(n_fermion_sector_size(hs, 2) == 1);
    hs.add(make_space_fermion(2));
    CHECK(n_fermion_sector_size(hs, 0) == 1);
    CHECK(n_fermion_sector_size(hs, 1) == 3);
    CHECK(n_fermion_sector_size(hs, 2) == 3);
    CHECK(n_fermion_sector_size(hs, 3) == 1);

    // Fermions and bosons
    hs.add(make_space_boson(4, 3));
    CHECK(n_fermion_sector_size(hs, 0) == 16);
    CHECK(n_fermion_sector_size(hs, 1) == 48);
    CHECK(n_fermion_sector_size(hs, 2) == 48);
    CHECK(n_fermion_sector_size(hs, 3) == 16);

    // Purely bosonic Hilbert space
    hs_type hs_b(make_space_boson(2, 0), make_space_boson(3, 1));
    CHECK(n_fermion_sector_size(hs_b, 0) == 32);
    CHECK(n_fermion_sector_size(hs_b, 1) == 0);
  }
}
