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
#ifndef LIBCOMMUTE_QOPERATOR_ELEMENTARY_SPACE_HPP_
#define LIBCOMMUTE_QOPERATOR_ELEMENTARY_SPACE_HPP_

#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

//
// qoperator objects act on states in a Hilbert space. That Hilbert space is
// formed as a direct product of smaller elementary spaces. For example,
// in the case of fermions one has as many elementary spaces as there are
// fermionic degrees of freedom, and each of those elementary spaces is
// 2-dimensional (|0> and |1>).
//

template<typename... IndexTypes>
class elementary_space {

public:

  using index_types = std::tuple<IndexTypes...>;

  // Value semantics
  elementary_space(IndexTypes const&... indices) : indices_(indices...) {}
  elementary_space(index_types const& indices) : indices_(indices) {}
  elementary_space(index_types && indices) : indices_(std::move(indices)) {}
  elementary_space(elementary_space const&) = default;
  elementary_space(elementary_space&&) noexcept = default;
  elementary_space& operator=(elementary_space const&) = default;
  elementary_space& operator=(elementary_space&&) noexcept = default;
  virtual ~elementary_space() {}

  // Make a smart pointer that manages a copy of this elementary space
  virtual std::unique_ptr<elementary_space> clone() const = 0;

  // Comparisons
  friend bool operator==(elementary_space const& es1,
                         elementary_space const& es2) {
    return es1.algebra_id() == es2.algebra_id() && es1.equal(es2);
  }
  friend bool operator!=(elementary_space const& es1,
                         elementary_space const& es2) {
    return !operator==(es1, es2);
  }
  friend bool operator<(elementary_space const& es1,
                        elementary_space const& es2) {
    if(es1.algebra_id() != es2.algebra_id())
      return es1.algebra_id() < es2.algebra_id();
    else
      return es1.less(es2);
  }
  friend bool operator>(elementary_space const& es1,
                        elementary_space const& es2) {
    if(es1.algebra_id() != es2.algebra_id())
      return es1.algebra_id() > es2.algebra_id();
    else
      return es1.greater(es2);
  }

  // ID of the algebra this elementary space is associated with
  virtual int algebra_id() const = 0;

  // The minimal number of binary digits needed to represent any state
  // in this elementary space
  virtual int n_bits() const = 0;

  // Indices accessor
  index_types const& indices() const { return indices_; }

protected:

  index_types indices_;

  // Check two elementary spaces for equality
  virtual bool equal(elementary_space const& es) const {
    return indices_ == es.indices_;
  }
  // Ordering
  virtual bool less(elementary_space const& es) const {
    return indices_ < es.indices_;
  }
  virtual bool greater(elementary_space const& es) const {
    return indices_ > es.indices_;
  }
};

}

#endif
