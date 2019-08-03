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
#ifndef LIBCOMMUTE_MONOMIAL_ACTION_HPP_
#define LIBCOMMUTE_MONOMIAL_ACTION_HPP_

#include "hilbert_space.hpp"
#include "../expression/monomial.hpp"

namespace libcommute {

// Action of a monomial on a quantum state vector
// TODO

class default_monomial_action {
public:

  template<typename ScalarType, typename... IndexTypes>
  default_monomial_action(monomial<IndexTypes...> const& m,
                          ScalarType&& coeff,
                          hilbert_space<IndexTypes...> const& hs) {
    // TODO
  }

  template<typename StateVector>
  inline void act(StateVector const& src, StateVector & dst) {
    // TODO
  }

  template<typename StateVector, typename... CoeffArgs>
  inline void act_at(StateVector const& src,
                     StateVector & dst,
                     CoeffArgs&&... args) {
    // TODO
  }

};

} // namespace libcommute

#endif
