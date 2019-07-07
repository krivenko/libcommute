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
#ifndef LIBCOMMUTE_MONOMIAL_HPP_
#define LIBCOMMUTE_MONOMIAL_HPP_

#include "generator.hpp"
#include "../utility"

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace libcommute {

//
// Monomial: a product of algebra generators
//

template <typename... IndexTypes>
class monomial {

  // Helper method for one of constructors
  template<typename GenType1, typename... GenTypesTail>
  void constructor_impl(GenType1 && generator, GenTypesTail&&... more_gens) {
    using gen1_t = typename std::remove_reference<GenType1>::type;
    generators_.emplace_back(new gen1_t(std::forward<GenType1>(generator)));
    constructor_impl(std::forward<GenTypesTail>(more_gens)...);
  }
  void constructor_impl() {}

public:

  using index_types = std::tuple<IndexTypes...>;
  using generator_type = generator<IndexTypes...>;
  using gen_ptr_type = std::unique_ptr<generator_type>;

  monomial() = default;

  // Construct from a list of generators
  template<typename... GenTypes>
  monomial(GenTypes&&... generators) {
    constructor_impl(std::forward<GenTypes>(generators)...);
  }

  // Construct from a list of pointers to generators
  monomial(std::initializer_list<generator_type*> generators) {
    for(auto p : generators) generators_.emplace_back(p->clone());
  }

  // Value semantics
  monomial(monomial const&) = default;
  monomial(monomial&&) noexcept = default;
  monomial& operator=(monomial const&) = default;
  monomial& operator=(monomial&&) noexcept = default;

  // Number of generators in this monomial
  inline size_t size() const { return generators_.size(); }

  // Iteration interface
  using const_iterator = typename std::vector<gen_ptr_type>::const_iterator;

  inline const_iterator begin() const noexcept { return generators_.begin(); }
  inline const_iterator cbegin() const noexcept { return generators_.cbegin(); }

  inline const_iterator end() const noexcept { return generators_.end(); }
  inline const_iterator cend() const noexcept { return generators_.cend(); }

  // Element access
  inline gen_ptr_type const& operator[](size_t n) const {
    return generators_[n];
  }

  // Equality
  friend bool operator==(monomial const& m1, monomial const& m2) {
    return m1.size() == m2.size() &&
      std::equal(m1.begin(), m1.end(), m2.begin(),
                 [](gen_ptr_type const& p1, gen_ptr_type const& p2) {
                   return *p1 == *p2;
                  }
      );
  }
  friend bool operator!=(monomial const& m1, monomial const& m2) {
    return !operator==(m1, m2);
  }

  // Ordering
  friend bool operator<(monomial const& m1, monomial const& m2) {
    if(m1.size() != m2.size())
      return m1.size() < m2.size();
    else {
      auto cmp = [](gen_ptr_type const& p1, gen_ptr_type const& p2) {
        return *p1 < *p2;
      };
      return std::lexicographical_compare(m1.begin(),
                                          m1.end(),
                                          m2.begin(),
                                          m2.end(), cmp);
    }
  }
  friend bool operator>(monomial const& m1, monomial const& m2) {
    if(m1.size() != m2.size())
      return m1.size() > m2.size();
    else {
      auto cmp = [](gen_ptr_type const& p1, gen_ptr_type const& p2) {
        return *p1 > *p2;
      };
      return std::lexicographical_compare(m1.begin(),
                                          m1.end(),
                                          m2.begin(),
                                          m2.end(), cmp);
    }
  }

  // Print to stream
  friend std::ostream& operator<<(std::ostream& os, monomial const& m) {
    auto it = m.begin(), end_it = m.end();
    auto next_it = it + 1;
    int power = 1;
    for(;it != end_it; ++it, ++next_it) {
      if(next_it == end_it || **next_it != **it) {
      os << (power > 1 ? "[" : "") << **it
         << (power > 1 ? "]^" + std::to_string(power) : "");
      power = 1;
    } else
      ++power;
    }
    return os;
  }

private:

  std::vector<gen_ptr_type> generators_;
};

} // namespace libcommute

#endif
