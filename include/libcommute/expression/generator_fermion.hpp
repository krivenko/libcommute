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
#ifndef LIBCOMMUTE_GENERATOR_FERMION_HPP_
#define LIBCOMMUTE_GENERATOR_FERMION_HPP_

#include "generator.hpp"
#include "../utility.hpp"

#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

// ID of the fermionic algebra
static constexpr int FERMION_ALGEBRA_ID = -3;

//
// Generator of the fermionic algebra
//

template<typename... IndexTypes>
class generator_fermion : public generator<IndexTypes...> {

  using base = generator<IndexTypes...>;

public:

  // Get ID of the algebra this generator belongs to
  virtual int algebra_id() const override { return FERMION_ALGEBRA_ID; }

  // Value semantics
  template<typename... Args>
  generator_fermion(bool dagger, Args&&... indices) :
    base(std::forward<Args>(indices)...), dagger_(dagger) {}
  generator_fermion() = delete;
  generator_fermion(generator_fermion const&) = default;
  generator_fermion(generator_fermion&&) noexcept = default;
  generator_fermion& operator=(generator_fermion const&) = default;
  generator_fermion& operator=(generator_fermion&&) noexcept = default;
  virtual ~generator_fermion() {}

  // Make a smart pointer that manages a copy of this generator
  virtual std::unique_ptr<base> clone() const override {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    return make_unique<generator_fermion>(*this);
  }

  // A square of C^+/C is vanising
  virtual int nilpotent_power() const override { return 2; }

  // c = -1, f(g) = \delta(g1, g2)
  virtual double
  commute(base const& g2,
          linear_function<std::unique_ptr<base>> & f) const override {
    assert(*this > g2);
    auto const& g2_ = dynamic_cast<generator_fermion const&>(g2);
    f.const_term = (this->indices_ == g2_.indices_ &&
                    dagger_ != g2_.dagger_) ?
                    1 : 0;
    f.terms.clear();
    return -1;
  }

  // Accessor
  inline int dagger() const { return dagger_; }

protected:
  // Creation or annihilation operator?
  bool dagger_;

  // Check two generators of the same algebra for equality
  virtual bool equal(base const& g) const override {
    auto const& f_g =  dynamic_cast<generator_fermion const&>(g);
    return dagger_ == f_g.dagger_ && this->indices_ == f_g.indices_;
  }

  // Ordering
  virtual bool less(base const& g) const override {
    auto const& f_g =  dynamic_cast<generator_fermion const&>(g);
    // Example: c+_1 < c+_2 < c+_3 < c_3 < c_2 < c_1
    if(this->dagger_ != f_g.dagger_)
      return (this->dagger_ > f_g.dagger_);
    else
      return this->dagger_ ? this->indices_ < f_g.indices_ :
                             this->indices_ > f_g.indices_;
  }
  virtual bool greater(base const& g) const override {
    auto const& f_g =  dynamic_cast<generator_fermion const&>(g);
    // Example: c_1 > c_2 > c_3 > c+_3 > c+_2 > c+_1
    if(this->dagger_ != f_g.dagger_)
      return (this->dagger_ < f_g.dagger_);
    else
      return this->dagger_ ? this->indices_ > f_g.indices_ :
                             this->indices_ < f_g.indices_;
  }

  // Print to stream
  virtual std::ostream & print(std::ostream & os) const override {
    os << "C" << (this->dagger_ ? "+" : "") << "(";
    print_tuple(os, this->indices_);
    return os << ")";
  }
};

// Convenience factory function
template<typename... IndexTypes>
inline generator_fermion<c_str_to_string_t<IndexTypes>...>
make_fermion(bool dagger, IndexTypes&&... indices) {
  return {dagger, std::forward<IndexTypes>(indices)...};
}

// Check if generator belongs to the fermionic algebra
template<typename... IndexTypes>
inline bool is_fermion(generator<IndexTypes...> const& gen) {
  return gen.algebra_id() == FERMION_ALGEBRA_ID;
}

} // namespace libcommute

#if __cplusplus >= 201703L
#include "dyn_indices.hpp"
namespace libcommute {
namespace dyn {

// Convenience factory functions for dynamic indices
template<typename... IndexTypes>
inline generator_fermion<dyn_indices>
make_fermion(bool dagger, IndexTypes&&... indices) {
  return {dagger, dyn_indices(std::forward<IndexTypes>(indices)...)};
}

} // namespace dyn
} // namespace libcommute
#endif

#endif
