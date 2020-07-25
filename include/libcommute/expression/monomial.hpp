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
#ifndef LIBCOMMUTE_EXPRESSION_MONOMIAL_HPP_
#define LIBCOMMUTE_EXPRESSION_MONOMIAL_HPP_

#include "generator.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"

#include <algorithm>
#include <cassert>
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

template<typename... IndexTypes>
class monomial {

  // Helper method for one of constructors
  template<typename GenType1, typename... GenTypesTail>
  void constructor_impl(GenType1 && generator, GenTypesTail&&... more_gens) {
    using gen1_t = typename std::remove_reference<GenType1>::type;
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    generators_.emplace_back(make_unique<gen1_t>(generator));
    constructor_impl(std::forward<GenTypesTail>(more_gens)...);
  }
  void constructor_impl() {}

  // Check if all provided types are derived from generator<IndexTypes...>
  template<typename... Types>
  using all_are_generators = all_derived_from<generator<IndexTypes...>,
                                              Types...>;

public:

  using index_types = std::tuple<IndexTypes...>;
  using generator_type = generator<IndexTypes...>;
  using gen_ptr_type = std::unique_ptr<generator_type>;

  // Construct empty (constant) monomial
  monomial() = default;

  // Construct from a list of >=1 generators
  template<typename... GenTypes,
           typename = typename std::enable_if<
              all_are_generators<GenTypes...>::value>::type
           >
  explicit monomial(GenTypes&&... generators) {
    constructor_impl(std::forward<GenTypes>(generators)...);
  }

  // Construct from a list of pointers to generators
  monomial(std::initializer_list<generator_type*> generators) {
    for(auto p : generators) generators_.emplace_back(p->clone());
  }

  // Construct from a list of smart pointers to generators
  monomial(std::initializer_list<gen_ptr_type> generators) {
    for(auto const& p : generators) generators_.emplace_back(p->clone());
  }

  // Value semantics
  monomial(monomial const& m) {
    generators_.reserve(m.generators_.size());
    for(gen_ptr_type const& g : m.generators_)
      generators_.emplace_back(g->clone());
  }
  monomial(monomial&&) noexcept = default;
  monomial& operator=(monomial const& m) {
    generators_.clear();
    for(gen_ptr_type const& g : m.generators_)
      generators_.emplace_back(g->clone());
    return *this;
  }
  monomial& operator=(monomial&&) noexcept = default;

  // Number of generators in this monomial
  inline size_t size() const { return generators_.size(); }

  // Is this monomial empty?
  inline bool empty() const { return generators_.empty(); }

  //
  // Iteration interface
  //

  // Constant iterator over generators comprising this monomial
  class const_iterator;

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

  using range_type = std::pair<const_iterator, const_iterator>;

  // Reverse iterator
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  inline const_reverse_iterator rbegin() const noexcept {
    return const_reverse_iterator(generators_.rbegin());
  }
  inline const_reverse_iterator crbegin() const noexcept {
    return const_reverse_iterator(generators_.crbegin());
  }

  inline const_reverse_iterator rend() const noexcept {
    return const_reverse_iterator(generators_.rend());
  }
  inline const_reverse_iterator crend() const noexcept {
    return const_reverse_iterator(generators_.crend());
  }

  // Element access
  inline generator_type const& operator[](size_t n) const {
    return *generators_[n];
  }
  inline generator_type& operator[](size_t n) { return *generators_[n]; }

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

  // Is this monomial sorted?
  inline bool is_sorted() const { return std::is_sorted(begin(), end()); }

  // Concatenate monomials, generators and ranges specified
  // by a pair of monomial iterators
  template<typename... PartTypes>
  friend monomial concatenate(PartTypes&&... parts) {
    monomial res;
    res.generators_.reserve(concat_parts_total_size(parts...));
    res.concat_impl(parts...);
    return res;
  }

  // Swap a pair of generators in this monomial
  void swap_generators(size_t n1, size_t n2) {
    assert(n1 < size());
    assert(n2 < size());
    std::swap(generators_[n1], generators_[n2]);
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

  //
  // Implementation bits of concatenate()
  //

  // Size of a generator is 1
  static constexpr size_t monomial_part_size(generator_type const& m) {
    return 1;
  }
  // Size of a complete monomial
  static size_t monomial_part_size(monomial const& m) { return m.size(); }
  // Size of a range within a monomial
  static size_t
  monomial_part_size(range_type const& range) {
    return std::distance(range.first, range.second);
  }

  // Total number of generators in a mixed list of generators,
  // monomials and monomial ranges
  template<typename P1, typename... PTail>
  static size_t concat_parts_total_size(P1&& p1, PTail&&... p_tail) {
    return monomial_part_size(std::forward<P1>(p1)) +
           concat_parts_total_size(std::forward<PTail>(p_tail)...);
  }
  template<typename P1>
  static size_t concat_parts_total_size(P1&& p1) {
    return monomial_part_size(std::forward<P1>(p1));
  }

  // Append generator
  void append_generators(generator_type const& g) {
    generators_.emplace_back(g.clone());
  }
  // Append generators from a monomial
  void append_generators(monomial const& m) {
    for(auto const& g : m.generators_)
      generators_.emplace_back(g->clone());
  }
  // Append generators from a monomial range
  void append_generators(range_type const& r) {
    for(auto it = r.first; it != r.second; ++it)
      generators_.emplace_back(it->clone());
  }

  // Append generators from a mixed list of monomials and monomial ranges
  template<typename P1, typename... PTail>
  void concat_impl(P1&& p1, PTail&&... p_tail) {
    append_generators(p1);
    concat_impl(std::forward<PTail>(p_tail)...);
  }
  template<typename P1>
  void concat_impl(P1&& p1) { append_generators(p1); }

  std::vector<gen_ptr_type> generators_;
};

// Constant iterator over generators comprising a monomial
template<typename... IndexTypes>
class monomial<IndexTypes...>::const_iterator {

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

} // namespace libcommute

#endif
