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

#include "generator.hpp"
#include "../utility.hpp"

#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

// ID of the bosonic algebra
static constexpr int BOSON_ALGEBRA_ID = -2;

//
// Generator of the bosonic algebra
//

template<typename... IndexTypes>
class generator_boson : public generator<IndexTypes...> {

  using base = generator<IndexTypes...>;

public:

  // Get ID of the algebra this generator belongs to
  virtual int algebra_id() const override { return BOSON_ALGEBRA_ID; }

  // Value semantics
  generator_boson(bool dagger, IndexTypes&&... indices) :
    base(std::forward<IndexTypes>(indices)...), dagger_(dagger) {}
  generator_boson() = delete;
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

protected:
  // Creation or annihilation operator?
  bool dagger_;

  // Check two generators of the same algebra for equality
  virtual bool equal(base const& g) const override {
    auto const& b_g =  dynamic_cast<generator_boson const&>(g);
    return dagger_ == b_g.dagger_ && this->indices_ == b_g.indices_;
  }

  // Ordering
  virtual bool less(base const& g) const override {
    auto const& b_g =  dynamic_cast<generator_boson const&>(g);
    // Example: a+_1 < a+_2 < a+_3 < a_3 < a_2 < a_1
    if(this->dagger_ != b_g.dagger_)
      return (this->dagger_ > b_g.dagger_);
    else
      return this->dagger_ ? this->indices_ < b_g.indices_ :
                             this->indices_ > b_g.indices_;
  }
  virtual bool greater(base const& g) const override {
    auto const& b_g =  dynamic_cast<generator_boson const&>(g);
    // Example: a_1 > a_2 > a_3 > a+_3 > a+_2 > a+_1
    if(this->dagger_ != b_g.dagger_)
      return (this->dagger_ < b_g.dagger_);
    else
      return this->dagger_ ? this->indices_ > b_g.indices_ :
                             this->indices_ < b_g.indices_;
  }

  // Print to stream
  virtual std::ostream & print(std::ostream & os) const override {
    os << "A" << (this->dagger_ ? "+" : "") << "(";
    print_tuple(os, this->indices_);
    return os << ")";
  }
};

// Convenience factory function
template<typename... IndexTypes>
inline generator_boson<typename c_str_to_string_t<IndexTypes>::type...>
make_boson(bool dagger, IndexTypes&&... indices) {
  return {dagger, std::forward<IndexTypes>(indices)...};
}

// Check if generator belongs to the bosonic algebra
template <typename... IndexTypes>
bool is_boson(generator<IndexTypes...> const& gen) {
  return gen.algebra_id() == BOSON_ALGEBRA_ID;
}

} // namespace libcommute

#endif
