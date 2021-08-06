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
// Custom algebra of Dirac gamma matrices: linear operator representation
//

#include "gamma.hpp" // Definition of 'gamma' algebra ID

#include <libcommute/loperator/elementary_space.hpp>
#include <libcommute/loperator/loperator.hpp>
#include <libcommute/loperator/monomial_action.hpp>

#include <algorithm>
#include <array>
#include <vector>

namespace libcommute {

//
// 4-dimensional elementary space gamma matrices act in
//

class elementary_space_gamma : public elementary_space<int> {

  using base = elementary_space<int>;

public:

  // Since all 4 gamma matrices act in the same elementary space,
  // we can initialize the base class with any number
  elementary_space_gamma() : base(0) {}
  elementary_space_gamma(elementary_space_gamma const&) = default;
  elementary_space_gamma(elementary_space_gamma &&) noexcept = default;
  elementary_space_gamma& operator=(elementary_space_gamma const&) = default;
  elementary_space_gamma& operator=(elementary_space_gamma &&) noexcept
    = default;
  ~elementary_space_gamma() final = default;

  // Virtual copy-constructor.
  // Make a smart pointer that manages a copy of this elementary space
  std::unique_ptr<base> clone() const final {
    return make_unique<elementary_space_gamma>(*this);
  }

  // Algebra ID, must be the same with generator_gamma::algebra_id()
  int algebra_id() const final { return libcommute::gamma; }

  // We need 2 bits to enumerate all 4 = 2^2 basis vectors
  int n_bits() const final { return 2; }
};

} // namespace libcommute

//
// Action of a product of gamma matrices (a gamma-matrix monomial)
// on a basis vector.
//

namespace libcommute {

template<> class monomial_action<libcommute::gamma> {

  // Indices of matrices in the product, right to left.
  // This order is chosen because the matrix on the right acts on a state first
  std::vector<int> index_sequence;

public:

  // m_range is a pair of iterators over a list of generator_gamma objects.
  // This range represents the product of gamma matrices we want to act with.
  monomial_action(monomial<int>::range_type const& m_range,
                  hilbert_space<int> const& hs) {
    // Iterate over matrices in the product, left to right
    for(auto it = m_range.first; it != m_range.second; ++it)
      // Collect indices of the matrices
      index_sequence.push_back(std::get<0>(it->indices()));
    // Reverse the order
    std::reverse(index_sequence.begin(), index_sequence.end());
  }

  //
  // Act on a basis state
  //
  // In the present implementation, libcommute requires all algebra generators
  // to be represented by generalized permutation matrices [1], i.e. matrices
  // with exactly one non-zero element in each row and each column.
  // Equivalently, any generator acting on a basis state must give a single
  // basis state multiplied by a constant.
  //
  // [1] https://en.wikipedia.org/wiki/Generalized_permutation_matrix
  //
  // 'index' is the index (a 64-bit unsigned integer) of the basis state we are
  // acting upon. It must also receive the index of the resulting basis state.
  // 'coeff' must be multiplied by the overall constant factor acquired as
  // a result of monomial action.
  //
  template<typename ScalarType>
  inline bool act(sv_index_type & index,
                  ScalarType & coeff) const {

    const std::complex<double> I(0,1);

    // Act with all matrices in the product, right to left
    for(int i : index_sequence) {
      switch(i) {
        case 0:
          // Action of \gamma^0
          // It is diagonal => 'index' does not change
          coeff *= std::array<double, 4>{1,1,-1,-1}[index];
          break;
        case 1:
          // Action of \gamma^1
          coeff *= std::array<double, 4>{-1,-1,1,1}[index];
          index = std::array<sv_index_type, 4>{3,2,1,0}[index];
          break;
        case 2:
          // Action of \gamma^2
          coeff *= std::array<std::complex<double>, 4>{-I,I,I,-I}[index];
          index = std::array<sv_index_type, 4>{3,2,1,0}[index];
          break;
        case 3:
          // Action of \gamma^3
          coeff *= std::array<std::complex<double>, 4>{-1,1,1,-1}[index];
          index = std::array<sv_index_type, 4>{2,3,0,1}[index];
          break;
      }
    }
    // This 'true' signals that the action result is not identically zero.
    // Returning 'false' in cases when it is zero would help improve
    // performance.
    return true;
  }
};

} // namespace libcommute

//
// Check new linear operator functionality
//

int main() {

  using namespace libcommute;

  // Hilbert space made of one elementary space for gamma matrices.
  hilbert_space<int> hs{elementary_space_gamma()};

  const std::complex<double> I(0,1);

  // Expression for \gamma^5
  auto gamma5 = I * make_gamma(0)
                  * make_gamma(1)
                  * make_gamma(2)
                  * make_gamma(3);

  // Linear operator representation of \gamma^5
  auto gamma5op = loperator<std::complex<double>,
                            libcommute::gamma>(gamma5, hs);

  //
  // Build the explicit matrix representation of \gamma^5
  //

  // Preallocate state vectors.
  std::vector<std::complex<double>> psi(4);
  std::vector<std::complex<double>> phi(4);

  // Iterate over basis states
  for(int i = 0; i < 4; ++i) {
    // Reset vectors |\psi> and |\phi> to zero.
    std::fill(psi.begin(), psi.end(), 0);
    std::fill(phi.begin(), phi.end(), 0);
    psi[i] = 1; // 'psi' is i-th basis vector now

    phi = gamma5op(psi);

    // Print the result
    std::cout << "gamma5|" << i << "> = ";
    for(int j = 0; j < 4; ++j) {
      std::cout << "+" << phi[j] << "|" << j << ">";
    }
    std::cout << std::endl;
  }

  return 0;
}
