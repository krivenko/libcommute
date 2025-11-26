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

#include "check_ordering.hpp"
#include "print_matcher.hpp"

#include <libcommute/expression/dyn_indices.hpp>

#include <string>
#include <vector>

using namespace libcommute::dynamic_indices;

TEST_CASE("Dynamic indices", "[dyn_indices]") {

  // Setup
  using dyn_ind_t = dyn_indices_generic<std::string, int, double>;

  std::vector<dyn_ind_t> all_ind = {dyn_ind_t(),

                                    dyn_ind_t("xxx"),
                                    dyn_ind_t("yyy"),
                                    dyn_ind_t(0),
                                    dyn_ind_t(1),
                                    dyn_ind_t(0.1),
                                    dyn_ind_t(0.2),

                                    dyn_ind_t("xxx", "xxx"),
                                    dyn_ind_t("xxx", "yyy"),
                                    dyn_ind_t("xxx", 0),
                                    dyn_ind_t("xxx", 1),
                                    dyn_ind_t("xxx", 0.1),
                                    dyn_ind_t("xxx", 0.2),
                                    dyn_ind_t("yyy", "xxx"),
                                    dyn_ind_t("yyy", "yyy"),
                                    dyn_ind_t("yyy", 0),
                                    dyn_ind_t("yyy", 1),
                                    dyn_ind_t("yyy", 0.1),
                                    dyn_ind_t("yyy", 0.2),
                                    dyn_ind_t(0, "xxx"),
                                    dyn_ind_t(0, "yyy"),
                                    dyn_ind_t(0, 0),
                                    dyn_ind_t(0, 1),
                                    dyn_ind_t(0, 0.1),
                                    dyn_ind_t(0, 0.2),
                                    dyn_ind_t(1, "xxx"),
                                    dyn_ind_t(1, "yyy"),
                                    dyn_ind_t(1, 0),
                                    dyn_ind_t(1, 1),
                                    dyn_ind_t(1, 0.1),
                                    dyn_ind_t(1, 0.2),
                                    dyn_ind_t(0.1, "xxx"),
                                    dyn_ind_t(0.1, "yyy"),
                                    dyn_ind_t(0.1, 0),
                                    dyn_ind_t(0.1, 1),
                                    dyn_ind_t(0.1, 0.1),
                                    dyn_ind_t(0.1, 0.2),
                                    dyn_ind_t(0.2, "xxx"),
                                    dyn_ind_t(0.2, "yyy"),
                                    dyn_ind_t(0.2, 0),
                                    dyn_ind_t(0.2, 1),
                                    dyn_ind_t(0.2, 0.1),
                                    dyn_ind_t(0.2, 0.2)};

  check_equality(all_ind);
  check_less_greater(all_ind);

  SECTION("Printing") {
    CHECK_THAT(dyn_ind_t(), Prints<dyn_ind_t>(""));
    CHECK_THAT(dyn_ind_t("xxx"), Prints<dyn_ind_t>("xxx"));
    CHECK_THAT(dyn_ind_t(0), Prints<dyn_ind_t>("0"));
    CHECK_THAT(dyn_ind_t(0.1), Prints<dyn_ind_t>("0.100000"));

    CHECK_THAT(dyn_ind_t("xxx", "yyy"), Prints<dyn_ind_t>("xxx,yyy"));
    CHECK_THAT(dyn_ind_t("xxx", 1), Prints<dyn_ind_t>("xxx,1"));
    CHECK_THAT(dyn_ind_t("xxx", 0.2), Prints<dyn_ind_t>("xxx,0.200000"));
    CHECK_THAT(dyn_ind_t(0, "yyy"), Prints<dyn_ind_t>("0,yyy"));
    CHECK_THAT(dyn_ind_t(0, 1), Prints<dyn_ind_t>("0,1"));
    CHECK_THAT(dyn_ind_t(0, 0.2), Prints<dyn_ind_t>("0,0.200000"));
    CHECK_THAT(dyn_ind_t(0.1, "yyy"), Prints<dyn_ind_t>("0.100000,yyy"));
    CHECK_THAT(dyn_ind_t(0.1, 1), Prints<dyn_ind_t>("0.100000,1"));
    CHECK_THAT(dyn_ind_t(0.1, 0.2), Prints<dyn_ind_t>("0.100000,0.200000"));
  }
}
