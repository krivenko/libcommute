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
#ifndef LIBCOMMUTE_GENERATOR_BOSON_HPP_
#define LIBCOMMUTE_GENERATOR_BOSON_HPP_

#include "algebra_ids.hpp"
#include "generator.hpp"
#include "../qoperator/basis_space_boson.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"

#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

//
// Generator of the bosonic algebra
//

template<typename... IndexTypes>
class generator_boson : public generator<IndexTypes...> {

  using base = generator<IndexTypes...>;
  using linear_function_t = typename base::linear_function_t;

public:

  // Get ID of the algebra this generator belongs to
  virtual int algebra_id() const override { return BOSON_ALGEBRA_ID; }

  // Value semantics
  template<typename... Args>
  generator_boson(bool dagger, Args&&... indices) :
    base(std::forward<Args>(indices)...), dagger_(dagger) {}
  generator_boson(generator_boson const&) = default;
  generator_boson(generator_boson&&) noexcept = default;
  generator_boson& operator=(generator_boson const&) = default;
  generator_boson& operator=(generator_boson&&) noexcept = default;
  virtual ~generator_boson() {}

  // Make a smart pointer that manages a copy of this generator
  virtual std::unique_ptr<base> clone() const override {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    return make_unique<generator_boson>(*this);
  }

  // c = 1, f(g) = \delta(g1, g2)
  virtual double
  commute(base const& g2, linear_function_t & f) const override {
    assert(*this > g2);
    auto const& g2_ = dynamic_cast<generator_boson const&>(g2);
    double delta = base::equal(g2) && dagger_ != g2_.dagger_;
    f.set(delta);
    return 1;
  }

  virtual bool collapse_power(int power, linear_function_t &) const override {
    assert(power >= 2);
    return false;
  }

  // Any power of this generator can be non-vanising
  virtual bool has_vanishing_power(int power) const override {
    assert(power >= 2);
    return false;
  }

  // Accessor
  inline int dagger() const { return dagger_; }

  // Return the Hermitian conjugate of this generator via f
  virtual void conj(linear_function_t & f) const override {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    f.set(0, make_unique<generator_boson>(!dagger_, base::indices_), 1);
  }

  // Make a basis space for bosons
  virtual std::unique_ptr<basis_space<IndexTypes...>>
  make_basis_space() const override {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    return make_unique<basis_space_boson<IndexTypes...>>(0, base::indices_);
  }

protected:
  // Creation or annihilation operator?
  bool dagger_;

  // Check two generators of the same algebra for equality
  virtual bool equal(base const& g) const override {
    auto const& b_g =  dynamic_cast<generator_boson const&>(g);
    return dagger_ == b_g.dagger_ && base::equal(g);
  }

  // Ordering
  virtual bool less(base const& g) const override {
    auto const& b_g =  dynamic_cast<generator_boson const&>(g);
    // Example: a+_1 < a+_2 < a+_3 < a_3 < a_2 < a_1
    if(this->dagger_ != b_g.dagger_)
      return (this->dagger_ > b_g.dagger_);
    else
      return this->dagger_ ? base::less(g) : base::greater(g);
  }
  virtual bool greater(base const& g) const override {
    auto const& b_g =  dynamic_cast<generator_boson const&>(g);
    // Example: a_1 > a_2 > a_3 > a+_3 > a+_2 > a+_1
    if(this->dagger_ != b_g.dagger_)
      return (this->dagger_ < b_g.dagger_);
    else
      return this->dagger_ ? base::greater(g) : base::less(g);
  }

  // Print to stream
  virtual std::ostream & print(std::ostream & os) const override {
    os << "A" << (this->dagger_ ? "+" : "") << "(";
    print_tuple(os, this->indices_);
    return os << ")";
  }
};

// Check if generator belongs to the bosonic algebra
template<typename... IndexTypes>
inline bool is_boson(generator<IndexTypes...> const& gen) {
  return gen.algebra_id() == BOSON_ALGEBRA_ID;
}

namespace static_indices {

// Convenience factory function
template<typename... IndexTypes>
inline generator_boson<c_str_to_string_t<IndexTypes>...>
make_boson(bool dagger, IndexTypes&&... indices) {
  return {dagger, std::forward<IndexTypes>(indices)...};
}

} // namespace libcommute::static_indices
} // namespace libcommute

#if __cplusplus >= 201703L
#include "dyn_indices.hpp"

namespace libcommute {
namespace dynamic_indices {

// Convenience factory functions for dynamic indices
template<typename... IndexTypes>
inline generator_boson<dyn_indices>
make_boson(bool dagger, IndexTypes&&... indices) {
  return {dagger, dyn_indices(std::forward<IndexTypes>(indices)...)};
}

} // namespace libcommute::dynamic_indices
} // namespace libcommute
#endif

#endif
