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
#ifndef LIBCOMMUTE_EXPRESSION_GENERATOR_HPP_
#define LIBCOMMUTE_EXPRESSION_GENERATOR_HPP_

#include "../utility.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

template<typename... IndexTypes> class elementary_space;

//
// Abstract algebra generator
//

template<typename... IndexTypes>
class generator {

public:

  using index_types = std::tuple<IndexTypes...>;

  // Linear combination of generators
  using linear_function_t = linear_function<std::unique_ptr<generator>>;

  template<typename... Args>
  generator(Args&&... indices) : indices_(std::forward<Args>(indices)...) {}
  generator(generator const&) = default;
  generator(generator&&) noexcept = default;
  generator& operator=(generator const&) = default;
  generator& operator=(generator&&) noexcept = default;
  virtual ~generator() = default;

  // Get ID of the algebra this generator belongs to
  virtual int algebra_id() const = 0;

  // Make a smart pointer that manages a copy of this generator
  virtual std::unique_ptr<generator> clone() const = 0;

  // Comparisons
  friend bool operator==(generator const& g1, generator const& g2) {
    return g1.algebra_id() == g2.algebra_id() && g1.equal(g2);
  }
  friend bool operator!=(generator const& g1, generator const& g2) {
    return !operator==(g1, g2);
  }
  friend bool operator<(generator const& g1, generator const& g2) {
    if(g1.algebra_id() != g2.algebra_id())
      return g1.algebra_id() < g2.algebra_id();
    else
      return g1.less(g2);
  }
  friend bool operator>(generator const& g1, generator const& g2) {
    if(g1.algebra_id() != g2.algebra_id())
      return g1.algebra_id() > g2.algebra_id();
    else
      return g1.greater(g2);
  }

  // Accessor
  inline index_types const& indices() const { return indices_; }

  // The following methods are called to simplify products of generators
  //
  // We assume that any pair of generators g1 and g2 satisfy
  //
  // g1 * g2 = c * g2 * g1 + f(g),
  //
  // where c is a real constant, and f(g) is a linear function of generators of
  // the same algebra. In particular, this condition is fulfilled by generators
  // of the Lie and Clifford algebras.

  // Given a pair g1 = *this and g2 such that g1 > g2, swap_with() must signal
  // what transformation g1 * g2 -> c * g2 * g1 + f(g) should be applied
  // to the product g1 * g2 to put it into the canonical order.
  // swap_with() returns the constant 'c' and writes the linear function f(g)
  // into its second argument. 'c' is allowed to be zero.
  virtual double
  swap_with(generator const& g2, linear_function_t & f) const = 0;

  // Given a pair g1 = *this and g2 such that g1 * g2 is in the canonical order
  // (g1 <= g2), optionally apply a simplifying transformation
  // g1 * g2 -> f(g).
  // If a simplification is actually possible, simplified_prod() must return
  // true and write the linear function f(g) into its second argument. Otherwise
  // return false.
  virtual bool
  simplify_prod(generator const& g2, linear_function_t & f) const {
    assert(!(*this > g2));
    return false;
  }

  // Given a generator g1 = *this and a power > 2, optionally apply
  // a simplifying transformation g1^power -> f(g).
  // If a simplification is actually possible, reduce_power() must return
  // true and write the linear function f(g) into its second argument. Otherwise
  // return false.
  //
  // N.B. Simplifications for power = 2 must be carried out by simplify_prod().
  virtual bool reduce_power(int power, linear_function_t & f) const {
    assert(power > 2);
    return false;
  }

  // Return the Hermitian conjugate of this generator via f
  virtual void conj(linear_function_t & f) const {
    f.set(0, clone(), 1.0);
  }

  // Stream output
  friend std::ostream & operator<<(std::ostream & os, generator const& g) {
    return g.print(os);
  }

protected:

  index_types indices_;

  // Check two generators of the same algebra for equality
  virtual bool equal(generator const& g) const {
    return indices_ == g.indices_;
  }
  // Ordering
  virtual bool less(generator const& g) const {
    return indices_ < g.indices_;
  }
  virtual bool greater(generator const& g) const {
    return indices_ > g.indices_;
  }
  // Print to stream
  virtual std::ostream & print(std::ostream & os) const {
    os << "g^" << algebra_id() << "(";
    print_tuple(os, this->indices_);
    return os << ")";
  }
};

// Check if g1 and g2 belong to the same algebra
// and call g1.swap_with(g2, f) accordingly
template<typename... IndexTypes>
inline double
swap_with(generator<IndexTypes...> const& g1,
          generator<IndexTypes...> const& g2,
          linear_function<std::unique_ptr<generator<IndexTypes...>>> & f) {
  if(g1.algebra_id() == g2.algebra_id()) {
    return g1.swap_with(g2, f);
  } else {
    // ** Generators of different algebras always commute **
    f.set(0);
    return 1;
  }
}

// Check if g1 and g2 belong to the same algebra
// and call g1.simplify_prod(g2, f) accordingly
template<typename... IndexTypes>
inline bool
simplify_prod(generator<IndexTypes...> const& g1,
              generator<IndexTypes...> const& g2,
              linear_function<std::unique_ptr<generator<IndexTypes...>>> & f) {
  if(g1.algebra_id() == g2.algebra_id()) {
    return g1.simplify_prod(g2, f);
  } else {
    // ** No simplification possible in a heterogeneous product **
    return false;
  }
}

} // namespace libcommute

#endif
