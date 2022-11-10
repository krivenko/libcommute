/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

//
// Implementation of the Virasoro algebra in libcommute
//
// Among other things, function main() checks recurrence relations (5)-(7) from
//
//   "A presentation for the Virasoro and super-Virasoro algebras",
//   D. B. Fairlie, J. Nuyts and C. K. Zachos ,
//   Commun. Math. Phys. 117, pp. 595–614 (1988),
//   https://doi.org/10.1007/BF01218387
//

#include <iostream>
#include <memory>
#include <tuple>

#include <libcommute/libcommute.hpp>

namespace libcommute {

// Define a new generator type by deriving from the abstract base class
// 'libcommute::generator<IndexTypes...>'. Generators of the Virasoro algebra
// L_n carry one integer index 'n'.
// Canonical order of generators in a monomial L_n * L_m * L_k * ... is such
// that n < m < k < ...
class generator_virasoro : public generator<int> {

  using base = generator<int>;
  using linear_function_t = typename base::linear_function_t;

  // For the sake of simplicity, we fix the central charge once and for all.
  static constexpr double central_charge = 2.0;

public:
  // This function must return a unique algebra ID shared by all generators
  // of a particular algebra.
  int algebra_id() const override {
    // Use the lowest algebra ID available to user-defined algebras
    return min_user_defined_algebra_id;
  }

  // Construct generator with a given index 'n'
  explicit generator_virasoro(int n) : base(n) {}
  // Standard constructors, assignments and destructor
  generator_virasoro(generator_virasoro const&) = default;
  generator_virasoro(generator_virasoro&&) noexcept = default;
  generator_virasoro& operator=(generator_virasoro const&) = default;
  generator_virasoro& operator=(generator_virasoro&&) noexcept = default;
  ~generator_virasoro() override = default;

  // Virtual copy-constructor: Make a smart pointer managing
  // a copy of this generator
  std::unique_ptr<base> clone() const override {
    return make_unique<generator_virasoro>(*this);
  }

  // Given a product L_m * L_n with m > n, transform it to the canonically
  // ordered form L_n * L_m plus some linear function of generators.
  // For the Virasoro algebra the transformation is
  //
  // L_m * L_n -> L_n * L_m + (m-n)*L_{m+n} + c(m^3 - m)\delta(m,-n)
  //
  // L_m will be passed to swap_with() as *this, i.e. L_m.swap_with(L_n, f).
  double swap_with(base const& L_n, linear_function_t& f) const override {
    // Ensure that L_m > L_n, or equivalently m > n.
    assert(*this > L_n);

    auto const& L_n_ = dynamic_cast<generator_virasoro const&>(L_n);

    // Extract indices from L_m and L_n
    int m = std::get<0>(base::indices());
    int n = std::get<0>(L_n_.indices());

    // Write linear terms of the transformed expressions into 'f'
    f.set(m == -n ? (central_charge * (m * m * m - m)) : 0, // Constant term
          make_unique<generator_virasoro>(m + n),           // L_{m+n}
          m - n // Coefficient in front of L_{m+n}
    );

    // Return coefficient in front of L_n * L_m in the transformed expression
    return 1;
  }

  // Given a product L_m * L_n with m <= n, optionally transform it to
  // some linear function of generators. For the Virasoro algebra such a
  // transformation exists for L_m = L_n,
  //
  // L_m * L_n -> 0
  //
  // L_m will be passed to simplify_prod() as *this,
  // i.e. L_m.simplify_prod(L_n, f).
  bool simplify_prod(base const& L_n, linear_function_t& f) const override {
    // Ensure that L_m <= L_n, or equivalently m <= n.
    assert(!(*this > L_n));

    if(*this == L_n) {
      // The transformed product is identically zero
      f.set(0);
      return true;
    } else
      return false; // No suitable transformation can be applied
  }

  // Hermitian conjugate: (L_n)^\dagger = L_{-n}
  void conj(linear_function_t& f) const override {
    int conj_n = -std::get<0>(base::indices());
    f.set(0, make_unique<generator_virasoro>(conj_n), 1);
  }

  // Print L_n to stream
  std::ostream& print(std::ostream& os) const override {
    int n = std::get<0>(base::indices());
    return os << "L(" << n << ")";
  }
};

// Convenience factory function to create expressions made of one
// monomial L_n.
expression<double, int> L(int n) {
  using ret_t = expression<double, int>;
  return ret_t(1.0, ret_t::monomial_t(generator_virasoro(n)));
}

} // namespace libcommute

using namespace libcommute;

int main() {

  // Check that L(0) is Hermitian and (L_n)^\dagger = L_{-n}
  std::cout << "conj(L(0)) = " << conj(L(0)) << std::endl;
  std::cout << "conj(L(1)) = " << conj(L(1)) << std::endl;

  // Check that L(n)^2 = 0 for a few n
  std::cout << "L(0) * L(0) = " << L(0) * L(0) << std::endl;
  std::cout << "L(1) * L(1) = " << L(1) * L(1) << std::endl;
  std::cout << "L(-1) * L(-1) = " << L(-1) * L(-1) << std::endl;

  // Check recurrence relations from Eq. (5)
  std::cout << "L_1 - (1/5)[L_3, L_{-2}] = "
            << (L(1) - (1.0 / 5) * (L(3) * L(-2) - L(-2) * L(3))) << std::endl;
  std::cout << "L_{-1} - (1/3)[L_1, L_{-2}] = "
            << (L(-1) - (1.0 / 3) * (L(1) * L(-2) - L(-2) * L(1))) << std::endl;
  std::cout << "L_2 - (1/4)[L_3, L_{-1}] = "
            << (L(2) - (1.0 / 4) * (L(3) * L(-1) - L(-1) * L(3))) << std::endl;
  std::cout << "L_0 - (1/2)[L_1, L_{-1}] = "
            << (L(0) - (1.0 / 2) * (L(1) * L(-1) - L(-1) * L(1))) << std::endl;

  // Check recurrence relation Eq. (6) for some higher positive n
  for(int n = 3; n < 10; ++n) {
    std::cout << "L_" << (n + 1) << " - (1/" << (n - 1) << ")"
              << "[L_" << n << ", L_1] = "
              << (L(n + 1) - (1.0 / (n - 1)) * (L(n) * L(1) - L(1) * L(n)))
              << std::endl;
  }
  // Check recurrence relation Eq. (7) for some higher negative n
  for(int n = 2; n < 10; ++n) {
    std::cout << "L_" << (-n - 1) << " - (1/(" << (1 - n) << "))"
              << "[L_" << -n << ", L_{-1}] = "
              << (L(-n - 1) - (1.0 / (1 - n)) * (L(-n) * L(-1) - L(-1) * L(-n)))
              << std::endl;
  }

  return 0;
}
