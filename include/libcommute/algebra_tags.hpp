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
#ifndef LIBCOMMUTE_ALGEBRA_TAGS_HPP_
#define LIBCOMMUTE_ALGEBRA_TAGS_HPP_

#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

//
// Predefined algebra tag types and respective algebra IDs
//

struct fermion { static constexpr int algebra_id() { return -3; } };
struct boson { static constexpr int algebra_id() { return -2; } };
struct spin { static constexpr int algebra_id() { return -1; } };

} // namespace libcommute

#endif
