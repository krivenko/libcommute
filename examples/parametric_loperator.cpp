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

//
// Hamiltonian of the displaced quantum harmonic oscillator with
// the displacement being a polynomial of an external parameter.
//
// Action of the Hamiltonian on the vacuum state is derived for a few values
// of the parameter using 'parametric_loperator'.
//

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <vector>

#include <libcommute/libcommute.hpp>

using namespace libcommute;

// Polynomial of one real variable - a (very) basic implementation.
// For production code, it is strongly advised to use optimized libraries
// such as 'polynomial' from Boost Math Toolkit.
class polynomial {

  // Polynomial coefficients {C_0, C_1, C_2, ..., C_n}.
  // P_n(x) = C_0 + C_1 x + C_2 x^2 + ... + C_n * x^n.
  std::vector<double> coefficients;

public:
  // Construct zero polynomial of degree n
  explicit polynomial(std::size_t n) : coefficients(n + 1, 0) {}

  // Construct from a list of coefficients
  explicit polynomial(std::initializer_list<double> coeffs)
    : coefficients(coeffs) {}

  // Degree of this polynomial
  std::size_t degree() const { return coefficients.size() - 1; }

  // Read-only access to the coefficients
  std::vector<double> const& coeffs() const { return coefficients; }

  // Evaluate P_n(x)
  double operator()(double x) const {
    double sum = 0;
    double x_n = 1;
    for(double c : coefficients) {
      sum += c * x_n;
      x_n *= x;
    }
    return sum;
  }

  // Addition of polynomials
  polynomial operator+(polynomial const& p) const {
    auto max_degree = std::max(degree(), p.degree());
    polynomial result(max_degree);
    for(std::size_t i = 0; i <= max_degree; ++i) {
      result.coefficients[i] = (i <= degree() ? coefficients[i] : 0) +
                               (i <= p.degree() ? p.coefficients[i] : 0);
    }
    return result;
  }

  // Subtraction of polynomials
  polynomial operator-(polynomial const& p) const {
    auto max_degree = std::max(degree(), p.degree());
    polynomial result(max_degree);
    for(std::size_t i = 0; i <= max_degree; ++i) {
      result.coefficients[i] = (i <= degree() ? coefficients[i] : 0) -
                               (i <= p.degree() ? p.coefficients[i] : 0);
    }
    return result;
  }

  // Multiplication of polynomials
  polynomial operator*(polynomial const& p) const {
    polynomial result(degree() + p.degree());
    for(std::size_t i = 0; i <= degree(); ++i) {
      for(std::size_t j = 0; j <= p.degree(); ++j) {
        result.coefficients[i + j] += coefficients[i] * p.coefficients[j];
      }
    }
    return result;
  }

  // Multiplication by a constant
  polynomial operator*(double x) const {
    polynomial res(*this);
    std::for_each(res.coefficients.begin(),
                  res.coefficients.end(),
                  [x](double& c) { c *= x; });
    return res;
  }

  // Multiplication by a constant pre-factor
  friend polynomial operator*(double x, polynomial const& p) { return p * x; }
};

// Specialize struct scalar_traits to let libcommute know about our polynomial
// type so that it can be used as a coefficient type of libcommute::expression.
namespace libcommute {

template <> struct scalar_traits<polynomial> {
  // Zero value test
  static bool is_zero(polynomial const& p) {
    return std::all_of(p.coeffs().begin(), p.coeffs().end(), [](double c) {
      return c == 0;
    });
  }
  // Make a constant polynomial from a double value
  static polynomial make_const(double x) { return polynomial{x}; }
};

} // namespace libcommute

int main() {

  //
  // Displaced harmonic oscillator
  //

  // Oscillator frequency \omega_0 = 2.
  polynomial w0{2};
  // Displacement as a linear function of an external parameter,
  // g(\lambda) = 3\lambda.
  polynomial g{0, 3};

  // Hamiltonian of the oscillator
  using namespace static_indices; // For a_dag() and a()
  auto H = 0.5 * w0 + w0 * (a_dag<polynomial>() - g) * (a<polynomial>() - g);

  // H will act in the Hilbert space 'hs'. Since we are working with a boson,
  // we have to truncate the Hilbert space and allow a finite number of
  // excitations (here, N = 0, 1, ..., 2^4 - 1).
  auto hs = make_hilbert_space(H, boson_es_constructor(4));

  // Parametric quantum operator (depends on parameter \lambda)
  auto Hop = make_param_loperator(H, hs);

  //
  // Check that
  // H|0> = \omega_0 (1/2 + g(\lambda)^2) |0> - \omega_0 g(\lambda) |1>
  // for a few values of \lambda.
  //

  // Initial state |0>
  std::vector<double> ket0(hs.dim());
  ket0[0] = 1.0;

  for(double lambda : {.0, 1.0, 2.0, 3.0}) {
    // Act with H(\lambda) upon |0>
    auto ket = Hop(ket0, lambda);

    // Print all non-zero elements of |ket>
    std::cout << "\\lambda = " << lambda << " :|ket> = ";
    for(unsigned int i = 0; i < ket.size(); ++i) {
      if(ket[i] != 0) std::cout << " +(" << ket[i] << ")|" << i << ">";
    }
    std::cout << '\n';
  }

  return 0;
}
