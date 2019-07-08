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
#include <iterator>
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
    generators_.emplace_back(make_unique<gen1_t>(generator));
    constructor_impl(std::forward<GenTypesTail>(more_gens)...);
  }
  void constructor_impl() {}

public:

  using index_types = std::tuple<IndexTypes...>;
  using generator_type = generator<IndexTypes...>;
  using gen_ptr_type = std::unique_ptr<generator_type>;

  // Construct from a list of generators
  template<typename... GenTypes>
  explicit monomial(GenTypes&&... generators) {
    constructor_impl(std::forward<GenTypes>(generators)...);
  }
  template<> explicit monomial<>() {};

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
  class const_iterator {
    using vector_it = typename std::vector<gen_ptr_type>::const_iterator;
    vector_it v_it_;

  public:

    using iterator_category = std::random_access_iterator_tag;
    using value_type = generator_type const&;
    using difference_type = std::ptrdiff_t;
    using pointer = std::unique_ptr<generator_type> const&;
    using reference = generator_type const&;

    explicit const_iterator(vector_it const& v_it) : v_it_(v_it) {}

    const_iterator() = default;
    const_iterator(const_iterator const&) = default;
    const_iterator(const_iterator&&) noexcept = default;
    const_iterator& operator=(const_iterator const&) = default;
    const_iterator& operator=(const_iterator&&) noexcept = default;

    // Increments
    const_iterator& operator++() { ++v_it_; return *this;}
    const_iterator operator++(int) {
      const_iterator retval = *this;
      ++(*this);
      return retval;
    }

    // Decrements
    const_iterator& operator--() { --v_it_; return *this;}
    const_iterator operator--(int) {
      const_iterator retval = *this;
      --(*this);
      return retval;
    }

    // Arithmetics
    const_iterator& operator+=(size_t n) { v_it_ += n ; return *this; }
    const_iterator& operator-=(size_t n) { v_it_ -= n ; return *this; }
    friend const_iterator operator+(const_iterator const& it, size_t n) {
      return const_iterator(it.v_it_ + n);
    }
    friend const_iterator operator+(size_t n, const_iterator const& it) {
      return const_iterator(it.v_it_ + n);
    }
    friend const_iterator operator-(const_iterator const& it, size_t n) {
      return const_iterator(it.v_it_ - n);
    }
    friend difference_type operator-(const_iterator const& it1,
                                     const_iterator const& it2) {
      return it1.v_it_ - it2.v_it_;
    }

    // Equality
    bool operator==(const_iterator const& it) const {return v_it_ == it.v_it_;}
    bool operator!=(const_iterator const& it) const {return !(*this == it);}

    // Ordering
    bool operator<(const_iterator const& it) const {return v_it_ < it.v_it_;}
    bool operator>(const_iterator const& it) const {return v_it_ > it.v_it_;}
    bool operator<=(const_iterator const& it) const {return v_it_ <= it.v_it_;}
    bool operator>=(const_iterator const& it) const {return v_it_ >= it.v_it_;}

    // Dereference
    reference operator*() const { return **v_it_; }
    pointer operator->() const { return *v_it_; }
    reference operator[](size_t n) const { return *v_it_[n]; }

    // swap()
    friend void swap(const_iterator& lhs, const_iterator& rhs) {
      std::swap(lhs.v_it_, rhs.v_it_);
    }
  };

  inline const_iterator begin() const noexcept {
    return const_iterator(generators_.begin());
  }
  inline const_iterator cbegin() const noexcept {
    return const_iterator(generators_.cbegin());
  }

  inline const_iterator end() const noexcept {
    return const_iterator(generators_.end());
  }
  inline const_iterator cend() const noexcept {
    return const_iterator(generators_.cend());
  }

  // Element access
  inline generator_type const& operator[](size_t n) const {
    return *generators_[n];
  }

  // Equality
  friend bool operator==(monomial const& m1, monomial const& m2) {
    return m1.size() == m2.size() &&
           std::equal(m1.begin(), m1.end(), m2.begin());
  }
  friend bool operator!=(monomial const& m1, monomial const& m2) {
    return !operator==(m1, m2);
  }

  // Ordering
  friend bool operator<(monomial const& m1, monomial const& m2) {
    if(m1.size() != m2.size())
      return m1.size() < m2.size();
    else {
      return std::lexicographical_compare(m1.begin(),
                                          m1.end(),
                                          m2.begin(),
                                          m2.end());
    }
  }
  friend bool operator>(monomial const& m1, monomial const& m2) {
    if(m1.size() != m2.size())
      return m1.size() > m2.size();
    else {
      return std::lexicographical_compare(m1.begin(),
                                          m1.end(),
                                          m2.begin(),
                                          m2.end(),
                                          std::greater<generator_type const&>()
      );
    }
  }

  // Print to stream
  friend std::ostream& operator<<(std::ostream& os, monomial const& m) {
    auto it = m.begin(), end_it = m.end();
    auto next_it = it + 1;
    int power = 1;
    for(;it != end_it; ++it, ++next_it) {
      if(next_it == end_it || *next_it != *it) {
        os << (power > 1 ? "[" : "") << *it
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
