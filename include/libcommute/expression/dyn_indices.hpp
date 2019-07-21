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
#ifndef LIBCOMMUTE_DYN_INDICES_HPP_
#define LIBCOMMUTE_DYN_INDICES_HPP_

#include <string>
#include <utility>
#include <variant>
#include <vector>

#if __cplusplus < 201703L
#error "This header file requires a C++-17 compliant compiler"
#endif

namespace libcommute {
namespace dynamic_indices {

//
// Dynamic sequence of indices
// Both length of the sequence and type of each index is defined at runtime
// with the index type restricted to be one of `IndexTypes`.
//

template<typename... IndexTypes> class dyn_indices_generic {
  std::vector<std::variant<IndexTypes...>> indices_;

public:

  // Value semantics
  dyn_indices_generic() = default;

  template<typename... Args> dyn_indices_generic(Args&&... args) {
    (indices_.emplace_back(std::forward<Args>(args)), ...);
  }

  dyn_indices_generic(dyn_indices_generic const&) = default;
  dyn_indices_generic(dyn_indices_generic&&) noexcept = default;
  dyn_indices_generic& operator=(dyn_indices_generic const&) = default;
  dyn_indices_generic& operator=(dyn_indices_generic&&) noexcept = default;

  // Equality
  inline friend bool operator==(dyn_indices_generic const& ind1,
                                dyn_indices_generic const& ind2) {
    return ind1.indices_ == ind2.indices_;
  }
  inline friend bool operator!=(dyn_indices_generic const& ind1,
                                dyn_indices_generic const& ind2) {
    return !operator==(ind1, ind2);
  }

  // Ordering
  inline friend bool operator<(dyn_indices_generic const& ind1,
                               dyn_indices_generic const& ind2) {
    if(ind1.indices_.size() != ind2.indices_.size())
      return ind1.indices_.size() < ind2.indices_.size();
    else
      return ind1.indices_ < ind2.indices_;
  }
  inline friend bool operator>(dyn_indices_generic const& ind1,
                               dyn_indices_generic const& ind2) {
    if(ind1.indices_.size() != ind2.indices_.size())
      return ind1.indices_.size() > ind2.indices_.size();
    else
      return ind1.indices_ > ind2.indices_;
  }

  // Stream output
  friend std::ostream & operator<<(std::ostream & os,
                                   dyn_indices_generic const& ind) {
    const size_t N = ind.indices_.size();
    for(size_t i = 0; i < N; ++i) {
      std::visit([&os](auto const& x) { os << x; }, ind.indices_[i]);
      if(i < N-1) os << ",";
    }
    return os;
  }

};

// Dynamic sequence of integer/string indices
using dyn_indices = dyn_indices_generic<int, std::string>;

} // namespace libcommute::dynamic_indices
} // namespace libcommute

#endif
