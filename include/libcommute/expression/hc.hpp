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
#ifndef LIBCOMMUTE_EXPRESSION_HC_HPP_
#define LIBCOMMUTE_EXPRESSION_HC_HPP_

#include "expression.hpp"

namespace libcommute {

//
// A placeholder object that adds/subtracts the Hermitian conjugate
// to/from an expression.
//
// Inspired by https://github.com/dafer45/TBTK
//
static constexpr struct {
} hc;

template <typename ScalarType, typename... IndexTypes>
inline expression<ScalarType, IndexTypes...>
operator+(expression<ScalarType, IndexTypes...> const& expr,
          decltype(hc) const&) {
  return expr + conj(expr);
}

template <typename ScalarType, typename... IndexTypes>
inline expression<ScalarType, IndexTypes...>
operator-(expression<ScalarType, IndexTypes...> const& expr,
          decltype(hc) const&) {
  return expr - conj(expr);
}

} // namespace libcommute

#endif
