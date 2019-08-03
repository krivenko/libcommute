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

// Quantum-mechanical operator acting in a Hilbert space
// TODO: take a list of algebra tag types instead of MonomialAction
template<typename MonomialAction> class qoperator {

public:

  qoperator() = delete;

  template<typename ScalarType, typename... IndexTypes>
  qoperator(expression<ScalarType, IndexTypes...> const& expr,
            hilbert_space<IndexTypes...> const& hs) {
    for(auto const& m : expr) {
      m_actions_.emplace_back(m.monomial, m.coeff, hs);
    }
  }

  // Value semantics
  qoperator(qoperator const&) = default;
  qoperator(qoperator&&) noexcept = default;
  qoperator& operator=(qoperator const&) = default;
  qoperator& operator=(qoperator&&) noexcept = default;

  // Act on state and return the resulting state.
  template<typename StateVector>
  inline StateVector operator()(StateVector const& sv) {
    StateVector res = zeros_like(sv);
    act_impl(sv, res);
    return res;
  }

  // Act on state `src` and return the resulting state via `dst`.
  template<typename StateVector>
  inline void operator()(StateVector const& src, StateVector & dst) {
    assert(src.size() == dst.size());
    set_zeros(dst);
    act_impl(src, dst);
  }

  // Act on state and return the resulting state.
  template<typename StateVector>
  inline StateVector operator*(StateVector const& sv) {
    return operator()(sv);
  }

  // Act on state and return the resulting state.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `args` as arguments to produce
  // the actual coefficient values.
  template<typename StateVector, typename... CoeffArgs>
  inline StateVector act_at(StateVector const& sv, CoeffArgs&&... args) {
    StateVector res = zeros_like(sv);
    act_at_impl(sv, res, std::forward<CoeffArgs>(args)...);
    return res;
  }

  // Act on state `src` and return the resulting state via `dst`.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `coeff_args` as arguments to produce
  // the actual coefficient values.
  template<typename StateVector, typename... CoeffArgs>
  inline void act_inplace_at(StateVector const& src,
                             StateVector & dst,
                             CoeffArgs&&... args) {
    set_zeros(dst);
    act_at_impl(src, dst, std::forward<CoeffArgs>(args)...);
  }

private:

  // Implemetation details of 'act_*' methods.
  template<typename StateVector>
  inline void act_impl(StateVector const& src, StateVector & dst) {
    for(auto const& ma : m_actions_) ma.act(src, dst);
  }

  template<typename StateVector, typename... CoeffArgs>
  inline void act_at_impl(StateVector const& src,
                          StateVector & dst,
                          CoeffArgs&&... args) {
    for(auto const& ma : m_actions_)
      ma.act(src, dst, std::forward<CoeffArgs>(args)...);
  }

  std::vector<MonomialAction> m_actions_;
};

// Factory function for qoperator
template<typename ScalarType, typename... IndexTypes>
qoperator<default_monomial_action>
make_qoperator(expression<ScalarType, IndexTypes...> const& expr,
               hilbert_space<IndexTypes...> const& hs) {
  return qoperator<default_monomial_action>(expr, hs);
}

} // namespace libcommute

#endif
