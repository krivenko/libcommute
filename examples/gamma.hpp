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

//
// Custom algebra of Dirac gamma matrices: algebra generators
//

#include <libcommute/expression/generator.hpp>

// For LIBCOMMUTE_MIN_USER_DEFINED_ALGEBRA_ID
#include <libcommute/algebra_ids.hpp>

#include <complex>
#include <iostream>

namespace libcommute {

//
// First, we define an ID for our new algebra
//

// Use the lowest algebra ID available to user-defined algebras
static constexpr int gamma = LIBCOMMUTE_MIN_USER_DEFINED_ALGEBRA_ID;

//
// Our generator type: a gamma matrix with one integer index
//
class generator_gamma : public generator<int> {

  using base = generator<int>;
  using linear_function_t = typename base::linear_function_t;

public:

  // Algebra ID of this generator
  int algebra_id() const final { return libcommute::gamma; }

  // Constructor: Just pass the index to the base class
  generator_gamma(int index) : base(index) {}

  // Virtual copy-constructor.
  // Make a smart pointer that manages a copy of this generator
  std::unique_ptr<base> clone() const final {
    // With C++14 or newer, libcommute::make_unique() will be resolved to
    // std::make_unique(). Otherwise, a custom implementation will be used.
    return make_unique<generator_gamma>(*this);
  }

  // This function will be called for g1 = *this and g2 such that g1 > g2.
  // We must transform the product g1 * g2 and put it into the
  // canonical order,
  //
  // g1 * g2 -> -g2 * g1 + 2\eta(g1, g2)
  double
  swap_with(base const& g2, linear_function_t & f) const final {

    // Do g1 and g2 have the same indices?
    bool diag = base::equal(g2);

    // Minkowski metric
    int index = std::get<0>(base::indices());
    double eta = diag ? (index == 0 ? 1 : -1) : 0;

    // Set f(g) to be the constant 2\eta(g1, g2)
    f.set(2 * eta);

    // Return coefficient in front of g2 * g1 in the transformed expression
    return -1;
  }

  // This function will be called for g1 = *this and g2 that are already
  // canonically ordered, g1 <= g2.
  //
  // It tries to simplify squares of gamma matrices,
  // (\gamma^0)^2 = I_4
  // (\gamma^k)^2 = -I_4 for k = 1,2,3
  bool
  simplify_prod(base const& g2, linear_function_t & f) const final {
    if(*this == g2) {
      // Replace the square with a constant
      int index = std::get<0>(base::indices());
      f.set(index == 0 ? 1 : -1);
      return true;
    } else
      // Not a square, cannot simplify
      return false;
  }

  // Hermitian conjugate of this generator as a linear function of generators.
  // \gamma^0 is Hermitian and \gamma^k are anti-Hermitian
  void conj(linear_function_t & f) const final {
    int index = std::get<0>(base::indices());
    if(index == 0) {
      // f(g) = 0 + 1*(*this)
      f.set(0, clone(), 1);
    } else {
      // f(g) = 0 + (-1)*(*this)
      f.set(0, clone(), -1);
    }
  }

  // Stream output function
  std::ostream & print(std::ostream & os) const final {
    int index = std::get<0>(base::indices());
    return os << "\\gamma^" << index;
  }
};

} // namespace libcommute

#include <libcommute/expression/expression.hpp>

namespace libcommute {

//
// Define a factory function to simplify construction of expressions
//

expression<std::complex<double>, int> make_gamma(int index) {
  using expr_t = expression<std::complex<double>, int>;
  // Return an expression made of one monomial with coefficient = 1.
  return expr_t(1.0, expr_t::monomial_t(generator_gamma(index)));
}

} // namespace libcommute
