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
#ifndef LIBCOMMUTE_LOPERATOR_ELEMENTARY_SPACE_SPIN_HPP_
#define LIBCOMMUTE_LOPERATOR_ELEMENTARY_SPACE_SPIN_HPP_

#include "../algebra_ids.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"

#include "elementary_space.hpp"

#include <cassert>
#include <cmath>

namespace libcommute {

//
// Elementary space generated by one spin degree of freedom
//

template<typename... IndexTypes>
class elementary_space_spin : public elementary_space<IndexTypes...> {

  using base = elementary_space<IndexTypes...>;

public:

  // Value semantics
  elementary_space_spin() = delete;
  template<typename... Args>
  elementary_space_spin(double spin, Args&&... indices)
    : base(std::forward<Args>(indices)...), multiplicity_(2*spin+1) {
    // Multiplicity has to be integer
    assert(2*spin == int(spin*2));
    n_bits_ = std::ceil(std::log2(multiplicity_));
  }
  elementary_space_spin(elementary_space_spin const&) = default;
  elementary_space_spin(elementary_space_spin&&) noexcept = default;
  elementary_space_spin& operator=(elementary_space_spin const&) = default;
  elementary_space_spin& operator=(elementary_space_spin&&) noexcept = default;

  // Make a smart pointer that manages a copy of this elementary space
  virtual std::unique_ptr<base> clone() const override {
    return make_unique<elementary_space_spin>(*this);
  }

  // ID of the algebra this elementary space is associated with
  virtual int algebra_id() const override { return spin; }

  // The minimal number of binary digits needed to represent any state
  // in this elementary space
  virtual int n_bits() const override { return n_bits_; }

protected:

  // Multiplicity, 2S+1
  int multiplicity_;

  // n_bits_ is the smallest positive number such that 2^{n_bits_} - 1
  // exceed multiplicity of the respective spin algebra.
  int n_bits_;

  // Equality
  virtual bool equal(base const& es) const override {
    auto const& es_s = dynamic_cast<elementary_space_spin const&>(es);
    return this->multiplicity_ == es_s.multiplicity_ && base::equal(es);
  }

  // Ordering
  virtual bool less(base const& es) const override {
    auto const& es_s =  dynamic_cast<elementary_space_spin const&>(es);
    if(this->multiplicity_ != es_s.multiplicity_)
      return this->multiplicity_ < es_s.multiplicity_;
    else
      return base::less(es);
  }
  virtual bool greater(base const& es) const override {
    auto const& es_s =  dynamic_cast<elementary_space_spin const&>(es);
    if(this->multiplicity_ != es_s.multiplicity_)
      return this->multiplicity_ > es_s.multiplicity_;
    else
      return base::greater(es);
  }
};

namespace static_indices {

// Convenience factory function
template<typename... IndexTypes>
inline elementary_space_spin<c_str_to_string_t<IndexTypes>...>
make_space_spin(double spin, IndexTypes&&... indices) {
  return {spin, std::forward<IndexTypes>(indices)...};
}

} // namespace libcommute::static_indices
} // namespace libcommute

#if __cplusplus >= 201703L
#include "../expression/dyn_indices.hpp"

namespace libcommute {
namespace dynamic_indices {

// Convenience factory function for dynamic indices
template<typename... IndexTypes>
inline elementary_space_spin<dyn_indices>
make_space_spin(double spin, IndexTypes&&... indices) {
  return {spin, dyn_indices(std::forward<IndexTypes>(indices)...)};
}

} // namespace libcommute::dynamic_indices
} // namespace libcommute
#endif

#endif
