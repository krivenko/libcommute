/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2025 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_LOPERATOR_HILBERT_SPACE_HPP_
#define LIBCOMMUTE_LOPERATOR_HILBERT_SPACE_HPP_

#include "../expression/expression.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"
#include "elementary_space.hpp"
#include "es_constructor.hpp"
#include "state_vector.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace libcommute {

// Range of bits in a bit string, [start;end]
using bit_range_t = std::pair<int, int>;

//
// Hilbert space as the ordered product of base spaces
//

template <typename... IndexTypes> class hilbert_space {

  using elementary_space_t = elementary_space<IndexTypes...>;
  using es_ptr_type = std::unique_ptr<elementary_space_t>;

  // Size of a Hilbert space must be representable by an integer with
  // at most this number of bits.
  static constexpr int max_n_bits =
      std::numeric_limits<sv_index_type>::digits - 1;

  // Compare two elementary_space_t objects wrapped in std::unique_ptr<T>
  struct less {
    inline bool operator()(es_ptr_type const& p1, es_ptr_type const& p2) const {
      return *p1 < *p2;
    }
  };

  // Helper method for one of constructors
  template <typename ESType1, typename... ESTypesTail>
  inline void constructor_impl(ESType1&& es, ESTypesTail&&... more_es) {
    add(std::forward<ESType1>(es));
    constructor_impl(std::forward<ESTypesTail>(more_es)...);
  }
  void constructor_impl() {}

  // Check if all provided types are derived from generator<IndexTypes...>
  template <typename... Types>
  using all_are_elementary_spaces =
      all_derived_from<elementary_space_t, Types...>;

public:
  using index_types = std::tuple<IndexTypes...>;

  //
  // Exception types
  //

  struct elementary_space_exists : public std::runtime_error {
    es_ptr_type elementary_space_ptr;
    explicit elementary_space_exists(elementary_space_t const& es)
      : std::runtime_error("Elementary space already exists"),
        elementary_space_ptr(es.clone()) {}
  };

  struct elementary_space_not_found : public std::runtime_error {
    es_ptr_type elementary_space_ptr;
    explicit elementary_space_not_found(elementary_space_t const& es)
      : std::runtime_error("Elementary space not found"),
        elementary_space_ptr(es.clone()) {}
  };

  struct hilbert_space_too_big : public std::runtime_error {
    int n_bits;
    explicit hilbert_space_too_big(int n_bits)
      : std::runtime_error("Hilbert space size is not representable "
                           "by a " +
                           std::to_string(max_n_bits + 1) +
                           "-bit integer (n_bits = " + std::to_string(n_bits) +
                           ")"),
        n_bits(n_bits) {}
  };

  // Construct empty space
  hilbert_space() = default;

  // Construct from individual elementary spaces
  template <typename... Args,
            typename = typename std::enable_if<
                all_are_elementary_spaces<Args...>::value>::type>
  explicit hilbert_space(Args&&... args) {
    constructor_impl(std::forward<Args>(args)...);
  }

  // Construct from a vector of pointers to elementary spaces
  explicit hilbert_space(
      std::vector<elementary_space_t*> const& elementary_spaces) {
    for(auto const* p : elementary_spaces)
      add(*p);
  }

  // Inspect polynomial expression `expr` and collect elementary spaces
  // associated to every algebra generator found in `expr`. Construction
  // of elementary spaces is performed by `es_constr` functor.
  template <typename ScalarType,
            typename ESConstructor = default_es_constructor>
  explicit hilbert_space(expression<ScalarType, IndexTypes...> const& expr,
                         ESConstructor&& es_constr = {}) {
    for(auto const& m : expr) {
      for(auto const& g : m.monomial) {
        elementary_spaces_.emplace(std::forward<ESConstructor>(es_constr)(g),
                                   bit_range_t(0, 0));
      }
    }
    recompute_bit_ranges();

    using value_type = typename decltype(elementary_spaces_)::value_type;
    dim_ *= std::accumulate(elementary_spaces_.begin(),
                            elementary_spaces_.end(),
                            sv_index_type(1),
                            [](sv_index_type d, value_type const& es) {
                              return d * es.first->dim();
                            });
  }

  // Value semantics
  hilbert_space(hilbert_space const& hs)
    : bit_range_end_(hs.bit_range_end_), dim_(hs.dim_) {
    for(auto const& es : hs.elementary_spaces_) {
      elementary_spaces_.emplace_hint(elementary_spaces_.end(),
                                      es.first->clone(),
                                      es.second);
    }
  }
  hilbert_space(hilbert_space&&) noexcept = default;
  hilbert_space& operator=(hilbert_space const& hs) {
    bit_range_end_ = hs.bit_range_end_;
    elementary_spaces_.clear();
    for(auto const& es : hs.elementary_spaces_) {
      elementary_spaces_.emplace_hint(elementary_spaces_.end(),
                                      es.first->clone(),
                                      es.second);
    }
    dim_ = hs.dim_;
    return *this;
  }
  hilbert_space& operator=(hilbert_space&&) noexcept = default;

  ~hilbert_space() = default;

  // Equality
  inline friend bool operator==(hilbert_space const& hs1,
                                hilbert_space const& hs2) {
    using value_type = typename decltype(elementary_spaces_)::value_type;
    return hs1.size() == hs2.size() &&
           std::equal(hs1.elementary_spaces_.begin(),
                      hs1.elementary_spaces_.end(),
                      hs2.elementary_spaces_.begin(),
                      [](value_type const& v1, value_type const& v2) {
                        return *v1.first == *v2.first && v1.second == v2.second;
                      });
  }
  inline friend bool operator!=(hilbert_space const& hs1,
                                hilbert_space const& hs2) {
    return !operator==(hs1, hs2);
  }

  // Append a new elementary space to the ordered product
  void add(elementary_space_t const& es) {
    auto r = elementary_spaces_.emplace(es.clone(), bit_range_t(0, 0));
    if(!r.second) throw elementary_space_exists(es);
    recompute_bit_ranges();
    dim_ *= es.dim();
  }

  // Is a given elementary space found in this Hilbert space?
  bool has(elementary_space_t const& es) const {
    return elementary_spaces_.count(es.clone()) == 1;
  }

  // Linear index of a given elementary space in this Hilbert space
  int index(elementary_space_t const& es) const {
    auto it = elementary_spaces_.find(es.clone());
    if(it == elementary_spaces_.end())
      throw elementary_space_not_found(es);
    else
      return std::distance(elementary_spaces_.begin(), it);
  }

  // Bit range spanned by elementary space
  bit_range_t bit_range(elementary_space_t const& es) const {
    auto it = elementary_spaces_.find(es.clone());
    if(it == elementary_spaces_.end())
      throw elementary_space_not_found(es);
    else
      return it->second;
  }

  // Is an elementary space with a given algebra ID found in this Hilbert space?
  bool has_algebra(int algebra_id) const {
    return algebra_bit_ranges_.count(algebra_id);
  }

  // Bit range spanned by algebra ID
  bit_range_t const& algebra_bit_range(int algebra_id) const {
    auto it = algebra_bit_ranges_.find(algebra_id);
    if(it == algebra_bit_ranges_.end())
      throw std::runtime_error("No elementary spaces with algebra ID " +
                               std::to_string(algebra_id));
    else
      return it->second;
  }

  // Number of elementary spaces
  std::size_t size() const { return elementary_spaces_.size(); }

  // The minimal number of binary digits needed to represent any state
  // in this Hilbert space
  int total_n_bits() const { return bit_range_end_ + 1; }

  // Minimal size of a StateVector object compatible with this Hilbert space
  sv_index_type vec_size() const { return sv_index_type(1) << total_n_bits(); }
  friend sv_index_type get_vec_size(hilbert_space const& hs) {
    return hs.vec_size();
  }

  // Dimension of this Hilbert space
  sv_index_type dim() const { return dim_; }
  friend sv_index_type get_dim(hilbert_space const& hs) { return hs.dim(); }

  // Does this Hilbert space have a non-power-of-two dimension?
  bool is_sparse() const { return vec_size() > dim(); }

  // Apply functor `f` to all basis state indices
  template <typename Functor>
  inline friend void foreach(hilbert_space const& hs, Functor&& f) {
    sv_index_type d = hs.vec_size();
    for(sv_index_type i = 0; i < d; ++i) {
      std::forward<Functor>(f)(i);
    }
  }

  // Return index of the product basis state, which decomposes as
  // |0> |0> ... |0> |n>_{es} |0> ... |0>.
  sv_index_type basis_state_index(elementary_space_t const& es,
                                  sv_index_type n) {
    assert(n < (sv_index_type(1) << es.n_bits()));
    return n << bit_range(es).first;
  }

private:
  // Recompute bit ranges in elementary_spaces_
  void recompute_bit_ranges() {
    algebra_bit_ranges_.clear();
    bit_range_end_ = -1;
    for(auto& es : elementary_spaces_) {
      int n_bits = es.first->n_bits();
      bit_range_t range(bit_range_end_ + 1, bit_range_end_ + n_bits);
      es.second = range;

      int algebra_id = es.first->algebra_id();
      auto algebra_range_it = algebra_bit_ranges_.find(algebra_id);
      if(algebra_range_it == algebra_bit_ranges_.end())
        algebra_bit_ranges_.emplace(algebra_id, range);
      else {
        if(range.first < algebra_range_it->second.first)
          algebra_range_it->second.first = range.first;
        if(range.second > algebra_range_it->second.second)
          algebra_range_it->second.second = range.second;
      }

      bit_range_end_ += n_bits;
    }
    if(bit_range_end_ >= max_n_bits)
      throw hilbert_space_too_big(bit_range_end_ + 1);
  }

  // List of base spaces in the product and their corresponding bit ranges
  std::map<es_ptr_type, bit_range_t, less> elementary_spaces_;

  // Bit ranges spanned by all elementary spaces
  // associated with the same algebra
  std::map<int, bit_range_t> algebra_bit_ranges_;

  // End of the bit range spanned by this Hilbert space
  int bit_range_end_ = -1;

  // Dimension of this Hilbert space
  sv_index_type dim_ = 1;
};

// Convenience factory function
template <typename ScalarType,
          typename... IndexTypes,
          typename ESConstructor = default_es_constructor>
inline hilbert_space<IndexTypes...>
make_hilbert_space(expression<ScalarType, IndexTypes...> const& expr,
                   ESConstructor&& es_constr = {}) {
  return hilbert_space<IndexTypes...>(expr,
                                      std::forward<ESConstructor>(es_constr));
}

} // namespace libcommute

#endif
