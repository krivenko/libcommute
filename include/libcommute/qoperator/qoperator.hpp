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

#include <utility>

namespace libcommute {

// qoperator is a quantum-mechanical operator that acts on a state vector.
// The following C++20 concept describes a valid state vector type.
//template<typename T>
//concept StateVector = requires(T v) {
//    { v.size() } -> std::convertible_to<std::size_t>;
//    { v[std::size_t{}] };
//};

// Quantum-mechanical operator acting on a state
template<typename ScalarType, typename... IndexTypes>
class qoperator {

public:

  using expression_type = expression<ScalarType, IndexTypes...>;

  qoperator() = delete;
  qoperator(expression_type const& expr) {
    // TODO
  }
  // TODO: fine-tuned constructor(s)

  // Value semantics
  qoperator(qoperator const&) = default;
  qoperator(qoperator&&) noexcept = default;
  qoperator& operator=(qoperator const&) = default;
  qoperator& operator=(qoperator&&) noexcept = default;

  // Apply operator to state and return the resulting state.
  template<typename StateVector>
  StateVector operator()(StateVector const& sv) {
    // TODO
  }

  // Apply operator to state `in` and return the resulting state via `out`.
  template<typename StateVector>
  void operator()(StateVector const& in, StateVector & out) {
    // TODO
  }

  // Apply operator to state and return the resulting state.
  template<typename StateVector>
  inline StateVector operator*(StateVector const& sv) {
    return operator()(sv);
  }

  // Apply operator to state and return the resulting state.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `coeff_args` as arguments to produce
  // the actual coefficient values.
  template<typename StateVector, typename... CoeffArgs>
  StateVector apply_at(StateVector const& sv, CoeffArgs&&... coeff_args) {
    // TODO
  }

  // Apply operator to state `in` and return the resulting state via `out`.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `coeff_args` as arguments to produce
  // the actual coefficient values.
  template<typename StateVector, typename... CoeffArgs>
  void apply_inplace_at(StateVector const& in,
                        StateVector & out,
                        CoeffArgs&&... coeff_args) {
    // TODO
  }
};

// Factory function for qoperator
template<typename ScalarType, typename... IndexTypes, typename... Args>
qoperator<ScalarType, IndexTypes...>
make_qoperator(expression<ScalarType, IndexTypes...> const& e, Args&&... args) {
  return qoperator<ScalarType, IndexTypes...>(e, std::forward<Args>(args)...);
}

} // namespace libcommute

#endif
