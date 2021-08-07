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
#ifndef LIBCOMMUTE_EXPRESSION_GENERATOR_SPIN_HPP_
#define LIBCOMMUTE_EXPRESSION_GENERATOR_SPIN_HPP_

#include "generator.hpp"
#include "../algebra_ids.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

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
  int algebra_id() const override { return libcommute::spin; }

  // Value semantics
  template<typename... Args>
  generator_spin(spin_component c, Args&&... indices) :
    base(std::forward<Args>(indices)...), multiplicity_(2), c_(c) {}
  template<typename... Args>
  generator_spin(double spin, spin_component c, Args&&... indices) :
    base(std::forward<Args>(indices)...),
    multiplicity_(static_cast<int>(2*spin+1)),
    c_(c) {
      // Multiplicity has to be integer
      assert(2*spin == int(spin*2));
    }

  generator_spin(generator_spin const&) = default;
  generator_spin(generator_spin&&) noexcept = default;
  generator_spin& operator=(generator_spin const&) = default;
  generator_spin& operator=(generator_spin&&) noexcept = default;
  ~generator_spin() override = default;

  // Make a smart pointer that manages a copy of this generator
  std::unique_ptr<base> clone() const override {
    return make_unique<generator_spin>(*this);
  }

  // Generators with different indices or multiplicities commute.
  double
  swap_with(base const& g2, linear_function_t & f) const override {
    assert(*this > g2);
    auto const& g2_ = dynamic_cast<generator_spin const&>(g2);
    if(base::equal(g2) && this->multiplicity_ == g2_.multiplicity_) {
      if(this->multiplicity_ == 2)
        return swap_with_spin_one_half(g2_, f);
      else
        return swap_with_higher_spin(g2_, f);
    } else {
      f.set(0);
      return 1;
    }
  }

  // Simplifications are possible for S=1/2
  // S_z * S_z = 1/4
  // S_+ * S_+ = 0
  // S_- * S_- = 0
  // S_+ * S_z = -1/2 S_+
  // S_- * S_z = 1/2 S_-
  // S_+ * S_- = 1/2 + S_z
  bool
  simplify_prod(base const& g2, linear_function_t & f) const override {
    assert(!(*this > g2));
    auto const& g2_ = dynamic_cast<generator_spin const&>(g2);
    if(!base::equal(g2) || this->multiplicity_ != 2) return false;

    if(this->c_ == g2_.c_) {
      f.set(this->c_ == spin_component::z ? 0.25 : 0);
    } else {
      if(g2_.c_ == spin_component::z) {
        if(this->c_ == spin_component::plus) {
          f.set(0, clone(), -0.5);
        } else { /// c_ == spin_component::minus
          f.set(0, clone(), 0.5);
        }
      } else { // c_ == spin_component::plus && g2.c_ == spin_component::minus
        f.set(0.5, g2_.clone(), 1);
        dynamic_cast<generator_spin&>(*f.terms.back().first).c_ =
          spin_component::z;
      }
    }

    return true;
  }

  bool reduce_power(int power, linear_function_t & f) const override {
    assert(power > 2);
    // Raising and lowering operators are nilpotent
    return c_ != spin_component::z && power >= multiplicity_;
  }

  // Accessors
  inline int multiplicity() const { return multiplicity_; }
  inline double spin() const { return (multiplicity_-1)/2.0; }
  inline spin_component component() const { return c_; }

  // Return the Hermitian conjugate of this generator via f
  void conj(linear_function_t & f) const override {
    spin_component new_c = (c_ == spin_component::z ? spin_component::z :
    (c_ == spin_component::plus ? spin_component::minus : spin_component::plus)
    );
    double s = (multiplicity_-1)/2.0;
    f.set(0, make_unique<generator_spin>(s, new_c, base::indices()), 1);
  }

private:

  // Multiplicity, 2S+1
  int multiplicity_;

  // Creation or annihilation operator?
  spin_component c_;

protected:
  // Check two generators of the same algebra for equality
  bool equal(base const& g) const override {
    auto const& b_g = dynamic_cast<generator_spin const&>(g);
    return multiplicity_ == b_g.multiplicity_ &&
                      c_ == b_g.c_ &&
                      base::equal(g);
  }

  // Ordering
  bool less(base const& g) const override {
    auto const& s_g = dynamic_cast<generator_spin const&>(g);
    // Example: S1/2+_1 < S1/2-_1 < S1/2z_1 < S1/2+_2 < S1/2-_2 < S1/2z_2 <
    //          S3/2+_1 < S3/2-_1 < S3/2z_1 < S3/2+_2 < S3/2-_2 < S3/2z_2
    if(this->multiplicity_ != s_g.multiplicity_)
      return (this->multiplicity_ < s_g.multiplicity_);
    else if(!base::equal(g))
      return base::less(g);
    else
      return this->c_ < s_g.c_;
  }
  bool greater(base const& g) const override {
    auto const& s_g = dynamic_cast<generator_spin const&>(g);
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
  std::ostream & print(std::ostream & os) const override {
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
    print_tuple(os, this->indices());
    return os << ")";
  }

private:

  //  S_- * S_+ = 1/2 - S_z
  //  S_z * S_+ = 1/2 S_+
  //  S_z * S_- = -1/2 S_-
  double swap_with_spin_one_half(generator_spin const& g2_,
                                 linear_function_t & f) const {
    if(c_ == spin_component::z) {
      if(g2_.c_ == spin_component::plus) {
        f.set(0, g2_.clone(), 0.5);
      } else { /// g2.c_ == spin_component::minus
        f.set(0, g2_.clone(), -0.5);
      }
    } else { // c_ == spin_component::minus && g2.c_ == spin_component::plus
      f.set(0.5, g2_.clone(), -1);
      dynamic_cast<generator_spin&>(*f.terms.back().first).c_ =
        spin_component::z;
    }
    return 0;
  }

  //  S_- * S_+ = S_+ * S_- - 2*S_z
  //  S_z * S_+ = S_+ * S_z + S_+
  //  S_z * S_- = S_- * S_z - S_-
  double swap_with_higher_spin(generator_spin const& g2_,
                               linear_function_t & f) const {

    if(c_ == spin_component::z) {
      if(g2_.c_ == spin_component::plus) {
        f.set(0, g2_.clone(), 1);
      } else { /// g2.c_ == spin_component::minus
        f.set(0, g2_.clone(), -1);
      }
    } else { // c_ == spin_component::minus && g2.c_ == spin_component::plus
      f.set(0, g2_.clone(), -2);
      dynamic_cast<generator_spin&>(*f.terms.back().first).c_ =
        spin_component::z;
    }
    return 1;
  }

};

// Check if generator belongs to the spin algebra
template<typename... IndexTypes>
inline bool is_spin(generator<IndexTypes...> const& gen) {
  return gen.algebra_id() == libcommute::spin;
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
