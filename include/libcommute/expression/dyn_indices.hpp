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
#ifndef LIBCOMMUTE_EXPRESSION_DYN_INDICES_HPP_
#define LIBCOMMUTE_EXPRESSION_DYN_INDICES_HPP_

#include "../metafunctions.hpp"

#include <iostream>
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

public:

  using indices_t = std::vector<std::variant<IndexTypes...>>;

private:

  indices_t indices_;

public:

  // Value semantics
  dyn_indices_generic() = default;
  explicit dyn_indices_generic(indices_t && indices) :
    indices_(std::move(indices)) {}
  explicit dyn_indices_generic(indices_t const& indices) : indices_(indices) {}

  template<typename Arg,
           typename = std::enable_if_t<
             (!std::is_same_v<remove_cvref_t<Arg>, indices_t>) &&
             (!std::is_same_v<remove_cvref_t<Arg>, dyn_indices_generic>)
             >
           >
  dyn_indices_generic(Arg&& arg) {
    indices_.emplace_back(std::forward<Arg>(arg));
  }

  template<typename... Args, typename = std::enable_if_t<(sizeof...(Args) > 1)>>
  dyn_indices_generic(Args&&... args) {
    (indices_.emplace_back(std::forward<Args>(args)), ...);
  }

  dyn_indices_generic(dyn_indices_generic const&) = default;
  dyn_indices_generic(dyn_indices_generic&&) noexcept = default;
  dyn_indices_generic& operator=(dyn_indices_generic const&) = default;
  dyn_indices_generic& operator=(dyn_indices_generic&&) noexcept = default;
  ~dyn_indices_generic() = default;

  // Index sequence length
  [[nodiscard]] size_t size() const { return indices_.size(); }

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

  // Reference to underlying sequence
  explicit operator indices_t const& () const {
    return indices_;
  }

  // Stream output
  friend std::ostream & operator<<(std::ostream & os,
                                   dyn_indices_generic const& ind) {
    const size_t N = ind.indices_.size();
    for(size_t i = 0; i < N; ++i) {
      std::visit([&os](auto const& x) { os << x; }, ind.indices_[i]);
      if(i + 1 < N) os << ",";
    }
    return os;
  }

};

// Dynamic sequence of integer/string indices
using dyn_indices = dyn_indices_generic<int, std::string>;

} // namespace libcommute::dynamic_indices
} // namespace libcommute

#endif
