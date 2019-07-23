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
#ifndef LIBCOMMUTE_GENERATOR_SPIN_HPP_
#define LIBCOMMUTE_GENERATOR_SPIN_HPP_

#include "generator.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"

#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

// ID of the spin/angular momentum algebra
static constexpr int SPIN_ALGEBRA_ID = -1;

//
// Generator of the spin algebra
//

// S_+, S_- or S_z
enum spin_component : int {plus = 0, minus = 1, z = 2};

template<typename... IndexTypes>
class generator_spin : public generator<IndexTypes...> {

  using base = generator<IndexTypes...>;
  using linear_function_t = typename base::linear_function_t;

public:

  // Get ID of the algebra this generator belongs to
  virtual int algebra_id() const override { return SPIN_ALGEBRA_ID; }

  // Value semantics
  template<typename... Args>
  generator_spin(spin_component c, Args&&... indices) :
    base(std::forward<Args>(indices)...), multiplicity_(2), c_(c) {}
  template<typename... Args>
  generator_spin(double spin, spin_component c, Args&&... indices) :
    base(std::forward<Args>(indices)...), multiplicity_(2*spin+1), c_(c) {
      // Multiplicity has to be integer
      assert(2*spin == int(spin*2));
    }

  generator_spin(generator_spin const&) = default;
  generator_spin(generator_spin&&) noexcept = default;
  generator_spin& operator=(generator_spin const&) = default;
  generator_spin& operator=(generator_spin&&) noexcept = default;
  virtual ~generator_spin() {}

  // Make a smart pointer that manages a copy of this generator
  virtual std::unique_ptr<base> clone() const override {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    return make_unique<generator_spin>(*this);
  }

  // Raising and lowering operators are nilpotent
  virtual std::pair<bool, double> has_constant_power(int power) const override {
    if(c_ == spin_component::z) {
      if(multiplicity_ == 2 && power > 0 && (power%2) == 0)
        return std::make_pair(true, double(1) / (std::uint64_t(1) << power));
      else
        return std::make_pair(false, 0);
    } else {
      return std::make_pair(power >= multiplicity_, 0);
    }
  }

  // Generators with different indices or multiplicities commute.
  // Same indices and multiplicity:
  //  S_- * S_+ = S_+ * S_- - 2*S_z
  //  S_z * S_+ = S_+ * S_z + S_+
  //  S_z * S_- = S_- * S_z - S_-
  virtual double commute(base const& g2, linear_function_t & f) const override {
    assert(*this > g2);
    auto const& g2_ = dynamic_cast<generator_spin const&>(g2);
    if(base::equal(g2) && this->multiplicity_ == g2_.multiplicity_) {
      if(c_ == spin_component::z) {
        if(g2_.c_ == spin_component::plus) {
          f.set(0, std::make_pair(g2_.clone(), 1));
        } else { /// g2.c_ == spin_component::minus
          f.set(0, std::make_pair(g2_.clone(), -1));
        }
      } else { // c_ == spin_component::minus && g2.c_ == spin_component::plus
        f.set(0, std::make_pair(g2_.clone(), -2));
        dynamic_cast<generator_spin&>(*f.terms.back().first).c_ =
          spin_component::z;
      }
    } else
      f.set(0);
    return 1;
  }

  // Accessors
  inline int multiplicity() const { return multiplicity_; }
  inline int spin() const { return (multiplicity_-1)/2; }
  inline spin_component component() const { return c_; }

  // Return the Hermitian conjugate of this generator via f
  virtual void conj(linear_function_t & f) const override {
    spin_component new_c = (c_ == spin_component::z ? spin_component::z :
    (c_ == spin_component::plus ? spin_component::minus : spin_component::plus)
    );
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    f.set(0, std::make_pair(make_unique<generator_spin>((multiplicity_-1)/2.0,
                                                        new_c,
                                                        base::indices_),
                            1)
    );
  }

protected:
  // Multiplicity, 2S+1
  int multiplicity_;

  // Creation or annihilation operator?
  spin_component c_;

  // Check two generators of the same algebra for equality
  virtual bool equal(base const& g) const override {
    auto const& b_g =  dynamic_cast<generator_spin const&>(g);
    return multiplicity_ == b_g.multiplicity_ &&
                      c_ == b_g.c_ &&
                      base::equal(g);
  }

  // Ordering
  virtual bool less(base const& g) const override {
    auto const& s_g =  dynamic_cast<generator_spin const&>(g);
    // Example: S1/2+_1 < S1/2-_1 < S1/2z_1 < S1/2+_2 < S1/2-_2 < S1/2z_2 <
    //          S3/2+_1 < S3/2-_1 < S3/2z_1 < S3/2+_2 < S3/2-_2 < S3/2z_2
    if(this->multiplicity_ != s_g.multiplicity_)
      return (this->multiplicity_ < s_g.multiplicity_);
    else if(!base::equal(g))
      return base::less(g);
    else
      return this->c_ < s_g.c_;
  }
  virtual bool greater(base const& g) const override {
    auto const& s_g =  dynamic_cast<generator_spin const&>(g);
    // Example: S3/2z_2 > S3/2-_2 > S3/2+_2 > S3/2z_1 > S3/2-_1 > S3/2+_1 >
    //          S1/2z_2 > S1/2-_2 > S1/2+_2 > S1/2z_1 > S1/2-_1 > S1/2+_1
    if(this->multiplicity_ != s_g.multiplicity_)
      return (this->multiplicity_ > s_g.multiplicity_);
    if(!base::equal(g))
      return base::greater(g);
    else
      return this->c_ > s_g.c_;
  }

  // Print to stream
  virtual std::ostream & print(std::ostream & os) const override {
    os << "S";
    if(multiplicity_ != 2) {
      if(multiplicity_ % 2 == 0)
        os << (multiplicity_ - 1) << "/2";
      else
        os << ((multiplicity_ - 1) / 2);
    }
    switch(this->c_) {
      case spin_component::plus: os << "+"; break;
      case spin_component::minus: os << "-"; break;
      case spin_component::z: os << "z"; break;
    }
    os << "(";
    print_tuple(os, this->indices_);
    return os << ")";
  }
};

// Check if generator belongs to the spin algebra
template<typename... IndexTypes>
inline bool is_spin(generator<IndexTypes...> const& gen) {
  return gen.algebra_id() == SPIN_ALGEBRA_ID;
}

namespace static_indices {

// Convenience factory functions
template<typename... IndexTypes>
inline generator_spin<c_str_to_string_t<IndexTypes>...>
make_spin(spin_component c, IndexTypes&&... indices) {
  return {c, std::forward<IndexTypes>(indices)...};
}

template<typename... IndexTypes>
inline generator_spin<c_str_to_string_t<IndexTypes>...>
make_spin(double spin, spin_component c, IndexTypes&&... indices) {
  return {spin, c, std::forward<IndexTypes>(indices)...};
}

} // namespace libcommute::static_indices
} // namespace libcommute

#if __cplusplus >= 201703L
#include "dyn_indices.hpp"

namespace libcommute {
namespace dynamic_indices {

// Convenience factory functions for dynamic indices
template<typename... IndexTypes>
inline generator_spin<dyn_indices>
make_spin(spin_component c, IndexTypes&&... indices) {
  return {c, dyn_indices(std::forward<IndexTypes>(indices)...)};
}

template<typename... IndexTypes>
inline generator_spin<dyn_indices>
make_spin(double spin, spin_component c, IndexTypes&&... indices) {
  return {spin, c, dyn_indices(std::forward<IndexTypes>(indices)...)};
}

} // namespace libcommute::dynamic_indices
} // namespace libcommute
#endif

#endif
