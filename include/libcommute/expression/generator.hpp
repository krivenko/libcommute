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
#ifndef LIBCOMMUTE_GENERATOR_HPP_
#define LIBCOMMUTE_GENERATOR_HPP_

#include "../utility.hpp"

#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

//
// Abstract algebra generator
//

template<typename... IndexTypes>
class generator {

public:

  using index_types = std::tuple<IndexTypes...>;

  template<typename... Args>
  generator(Args&&... indices) : indices_(std::forward<Args>(indices)...) {}
  generator() = delete;
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

  // Stream output
  friend std::ostream & operator<<(std::ostream & os, generator const& g) {
    return g.print(os);
  }

protected:

  index_types indices_;

  // Check two generators of the same algebra for equality
  virtual bool equal(generator const& g) const = 0;
  // Ordering
  virtual bool less(generator const& g) const = 0;
  virtual bool greater(generator const& g) const = 0;
  // Print to stream
  virtual std::ostream & print(std::ostream & os) const = 0;
};

} // namespace libcommute

#endif
