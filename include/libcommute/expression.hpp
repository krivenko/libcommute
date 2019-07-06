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
#ifndef LIBCOMMUTE_EXPRESSION_HPP_
#define LIBCOMMUTE_EXPRESSION_HPP_

#include "generator.hpp"

#include <complex>
#include <map>
#include <tuple>

//
// Polynomial expression involving quantum-mechanical operators
//

namespace libcommute {

template <typename ScalarType, typename... IndexTypes>
class expression {

public:

  using scalar_type = ScalarType;
  using index_types = std::tuple<IndexTypes...>;

  // Value semantics
  expression() = default;
  expression(expression const&) = default;
  expression(expression&&) noexcept = default;
  expression& operator=(expression const&) = default;
  expression& operator=(expression&&) noexcept = default;

  // TODO

};

// Aliases for specific scalar types
template<typename... IndexTypes>
using expression_real = expression<double, IndexTypes...>;
template<typename... IndexTypes>
using expression_complex = expression<std::complex<double>, IndexTypes...>;

} // namespace libcommute

#endif
