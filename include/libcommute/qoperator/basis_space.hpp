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
#ifndef LIBCOMMUTE_BASIS_SPACE_HPP_
#define LIBCOMMUTE_BASIS_SPACE_HPP_

#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

//
// qoperator objects act on states in a Hilbert space. That Hilbert space is
// formed as a direct product of smaller basis spaces. For example, in the case
// of fermions one has as many basis spaces as there are fermionic degrees of
// freedom, and each of those basis spaces is 2-dimensional (|0> and |1>).
//

template<typename... IndexTypes>
class basis_space {

public:

  using index_types = std::tuple<IndexTypes...>;

  // Value symantics
  basis_space() = delete;
  basis_space(IndexTypes const&... indices) : indices_(indices...) {}
  basis_space(IndexTypes&&... indices) : indices_(std::move(indices)...) {}
  basis_space(basis_space const&) = default;
  basis_space(basis_space&&) noexcept = default;
  basis_space& operator=(basis_space const&) = default;
  basis_space& operator=(basis_space&&) noexcept = default;

  // Make a smart pointer that manages a copy of this basis space
  virtual std::unique_ptr<basis_space> clone() const = 0;

  // Comparisons
  friend bool operator==(basis_space const& bs1, basis_space const& bs2) {
    return bs1.algebra_id() == bs2.algebra_id() && bs1.equal(bs2);
  }
  friend bool operator!=(basis_space const& bs1, basis_space const& bs2) {
    return !operator==(bs1, bs2);
  }
  friend bool operator<(basis_space const& bs1, basis_space const& bs2) {
    if(bs1.algebra_id() != bs2.algebra_id())
      return bs1.algebra_id() < bs2.algebra_id();
    else
      return bs1.less(bs2);
  }
  friend bool operator>(basis_space const& bs1, basis_space const& bs2) {
    if(bs1.algebra_id() != bs2.algebra_id())
      return bs1.algebra_id() > bs2.algebra_id();
    else
      return bs1.greater(bs2);
  }

  // ID of the algebra this basis space is associated with
  virtual int algebra_id() const = 0;

  // The minimal number of binary digits needed to represent any state
  // in this basis space
  virtual int n_bits() const = 0;

  // Indices accessor
  index_types const& indices() const { return indices_; }

protected:

  index_types indices_;

  // Check two basis spaces for equality
  virtual bool equal(basis_space const& bs) const {
    return indices_ == bs.indices_;
  }
  // Ordering
  virtual bool less(basis_space const& bs) const {
    return indices_ < bs.indices_;
  }
  virtual bool greater(basis_space const& bs) const {
    return indices_ > bs.indices_;
  }
};

}

#endif
