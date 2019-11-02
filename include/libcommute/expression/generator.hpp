/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_EXPRESSION_GENERATOR_HPP_
#define LIBCOMMUTE_EXPRESSION_GENERATOR_HPP_

#include "../utility.hpp"

#include <stdexcept>
#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

template<typename... IndexTypes> class basis_space;

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
  virtual ~generator() {}

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

  // We assume that any pair of generators g1 and g2 satisfy
  //
  // g1 * g2 = c * g2 * g1 + f(g),
  //
  // where c is a constant, and f(g) is a linear function of generators of
  // the same algebra. In particular, this condition is fulfilled by generators
  // of the Lie and Clifford algebras.
  //
  // Given a pair g1 = *this and g2, commute() must return c and write
  // f(g) into its second argument. For the sake of optimization, it is
  // required that g1 > g2.
  virtual double commute(generator const& g2, linear_function_t & f) const = 0;

  // If the given power of this generator can be represented as a linear
  // function of generators, return true and write the linear function into f.
  // The power must be greater than one.
  virtual bool collapse_power(int power, linear_function_t & f) const {
    assert(power >= 2);
    return false;
  }

  // Does the given power of this generator vanish?
  // The power must be greater than one.
  virtual bool has_vanishing_power(int power) const {
    assert(power >= 2);
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
// and call g1.commute(g2) accordingly
template<typename... IndexTypes>
double commute(generator<IndexTypes...> const& g1,
               generator<IndexTypes...> const& g2,
               linear_function<std::unique_ptr<generator<IndexTypes...>>> & f) {
  // ** Generators of different algebras always commute **
  if(g1.algebra_id() != g2.algebra_id()) {
    f.set(0);
    return 1;
  } else {
    return g1.commute(g2, f);
  }
}

} // namespace libcommute

#endif
