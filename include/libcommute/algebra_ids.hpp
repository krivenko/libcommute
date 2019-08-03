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
#ifndef LIBCOMMUTE_ALGEBRA_IDS_HPP_
#define LIBCOMMUTE_ALGEBRA_IDS_HPP_

#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

//
// Predefined algebra IDs
//

// ID of the fermionic algebra
static constexpr int FERMION_ALGEBRA_ID = -3;
// ID of the bosonic algebra
static constexpr int BOSON_ALGEBRA_ID = -2;
// ID of the spin/angular momentum algebra
static constexpr int SPIN_ALGEBRA_ID = -1;

} // namespace libcommute

#endif
