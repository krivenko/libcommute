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
#ifndef LIBCOMMUTE_QOPERATOR_HPP_
#define LIBCOMMUTE_QOPERATOR_HPP_

#include <libcommute/expression/expression.hpp>
#include "hilbert_space.hpp"
#include "monomial_action.hpp"
#include "../utility.hpp"

#include <utility>
#include <vector>

namespace libcommute {

// Quantum-mechanical operator acting on a state
template<typename ScalarType, typename... IndexTypes>
class qoperator {

public:

  using expression_type = expression<ScalarType, IndexTypes...>;
  using hilbert_space_type = hilbert_space<IndexTypes...>;

  qoperator() = delete;
  qoperator(expression_type const& expr) {
    // TODO
  }
  qoperator(expression_type const& expr, hilbert_space_type const& hs) {
    // TODO
  }

  // Value semantics
  qoperator(qoperator const&) = default;
  qoperator(qoperator&&) noexcept = default;
  qoperator& operator=(qoperator const&) = default;
  qoperator& operator=(qoperator&&) noexcept = default;

  // Apply operator to state and return the resulting state.
  template<typename StateVector>
  inline StateVector operator()(StateVector const& sv) {
    StateVector res = zeros_like(sv);
    apply_impl(sv, res);
    return res;
  }

  // Apply operator to state `src` and return the resulting state via `dst`.
  template<typename StateVector>
  inline void operator()(StateVector const& src, StateVector & dst) {
    assert(src.size() == dst.size());
    set_zeros(dst);
    apply_impl(src, dst);
  }

  // Apply operator to state and return the resulting state.
  template<typename StateVector>
  inline StateVector operator*(StateVector const& sv) {
    return operator()(sv);
  }

  // Apply operator to state and return the resulting state.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `args` as arguments to produce
  // the actual coefficient values.
  template<typename StateVector, typename... CoeffArgs>
  inline StateVector apply_at(StateVector const& sv, CoeffArgs&&... args) {
    StateVector res = zeros_like(sv);
    apply_at_impl(sv, res, std::forward<CoeffArgs>(args)...);
    return res;
  }

  // Apply operator to state `src` and return the resulting state via `dst`.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `coeff_args` as arguments to produce
  // the actual coefficient values.
  template<typename StateVector, typename... CoeffArgs>
  inline void apply_inplace_at(StateVector const& src,
                               StateVector & dst,
                               CoeffArgs&&... args) {
    set_zeros(dst);
    apply_at_impl(src, dst, std::forward<CoeffArgs>(args)...);
  }

private:

  // Implemetation details of 'apply' methods.
  template<typename StateVector>
  inline void apply_impl(StateVector const& src, StateVector & dst) {
    // TODO
  }

  template<typename StateVector, typename... CoeffArgs>
  inline void apply_at_impl(StateVector const& src,
                            StateVector & dst,
                            CoeffArgs&&... args) {
    // TODO
  }

  std::vector<monomial_action> m_actions_;
};

// Factory function for qoperator
template<typename ScalarType, typename... IndexTypes, typename... Args>
qoperator<ScalarType, IndexTypes...>
make_qoperator(expression<ScalarType, IndexTypes...> const& e, Args&&... args) {
  return qoperator<ScalarType, IndexTypes...>(e, std::forward<Args>(args)...);
}

} // namespace libcommute

#endif
