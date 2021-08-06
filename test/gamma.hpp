/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#ifndef LIBCOMMUTE_TEST_GAMMA_HPP_
#define LIBCOMMUTE_TEST_GAMMA_HPP_

#include <libcommute/algebra_ids.hpp>
#include <libcommute/expression/generator.hpp>

#include <cassert>

//
// Implement algebra of 4-dimensional \Gamma-matrices
//

namespace libcommute {

// Use the lowest algebra ID available to user-defined algebras
static constexpr int gamma = LIBCOMMUTE_MIN_USER_DEFINED_ALGEBRA_ID;

class generator_gamma : public generator<int> {

  using base = generator<int>;
  using linear_function_t = typename base::linear_function_t;

public:

  int algebra_id() const override { return libcommute::gamma; }

  // Value semantics
  generator_gamma(int index) : base(index) {}
  generator_gamma(generator_gamma const&) = default;
  generator_gamma(generator_gamma&&) noexcept = default;
  generator_gamma& operator=(generator_gamma const&) = default;
  generator_gamma& operator=(generator_gamma&&) noexcept = default;
  ~generator_gamma() final = default;

  // Make a smart pointer that manages a copy of this generator
  std::unique_ptr<base> clone() const final {
    return make_unique<generator_gamma>(*this);
  }

  // c = -1, f(g) = 2\eta(g1, g2)
  double
  swap_with(base const& g2, linear_function_t & f) const final {
    assert(*this > g2);
    bool diag = base::equal(g2);
    f.set(diag * (std::get<0>(indices_) == 0 ? 2 : -2));
    return -1;
  }

  // (\Gamma^0)^2 = I_4
  // (\Gamma^k)^2 = -I_4 for k=1,2,3
  bool
  simplify_prod(base const& g2, linear_function_t & f) const final {
    if(*this == g2) {
      f.set(std::get<0>(indices_) == 0 ? 1 : -1);
      return true;
    } else
      return false;
  }

  // Gamma^0 is Hermitian and Gamma^k are anti-Hermitian
  void conj(linear_function_t & f) const final {
    f.set(0, clone(), std::get<0>(indices_) == 0 ? 1 : -1);
  }
};

} // namespace libcommute

#endif
