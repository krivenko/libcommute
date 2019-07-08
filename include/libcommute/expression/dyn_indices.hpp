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
#ifndef LIBCOMMUTE_DYN_INDICES_HPP_
#define LIBCOMMUTE_DYN_INDICES_HPP_

#include <string>
#include <variant>

// TODO: #error if included in C++ < 17
// TODO: document

namespace libcommute {

template<typename... IndexTypes> class dyn_indices_generic {
  // TODO: std::variant<IndexTypes...>
};

using dyn_indices = dyn_indices_generic<int, std::string>;

} // namespace libcommute

#endif
