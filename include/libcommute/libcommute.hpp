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
#ifndef LIBCOMMUTE_LIBCOMMUTE_HPP_
#define LIBCOMMUTE_LIBCOMMUTE_HPP_

//
// Main header of libcommute
//

#include "version.hpp"

#include "expression/generator_fermion.hpp"
#include "expression/generator_boson.hpp"
#include "expression/generator_spin.hpp"
#include "expression/expression.hpp"
#include "expression/factories.hpp"
#include "qoperator/basis_space_fermion.hpp"
#include "qoperator/basis_space_boson.hpp"
#include "qoperator/basis_space_spin.hpp"
#include "qoperator/qoperator.hpp"
#include "qoperator/space_partition.hpp"

// C++17-only headers
#if __cplusplus >= 201703L
#include "expression/dyn_indices.hpp"
#include "expression/factories_dyn.hpp"
#endif

#endif
