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
#ifndef LIBCOMMUTE_HILBERT_SPACE_HPP_
#define LIBCOMMUTE_HILBERT_SPACE_HPP_

#include "basis_space.hpp"
#include "../expression/expression.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"

#include <algorithm>
#include <exception>
#include <map>
#include <utility>
#include <vector>

namespace libcommute {

//
// Hilbert space as the ordered product of base spaces
//

template<typename... IndexTypes>
class hilbert_space {

  using basis_space_t = basis_space<IndexTypes...>;
  using bs_ptr_type = std::unique_ptr<basis_space_t>;

public:

  // Compare two basis_space_t objects wrapped in std::unique_ptr<T>
  struct less {
    inline bool operator()(bs_ptr_type const& p1, bs_ptr_type const& p2) const {
      return *p1 < *p2;
    }
  };

  // Helper method for one of constructors
  template<typename BSType1, typename... BSTypesTail>
  inline void add_impl(BSType1 && bs, BSTypesTail&&... more_bs) {
    using bs1_t = typename std::remove_reference<BSType1>::type;
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif

    int n_bits = bs.n_bits();
    auto r = basis_spaces_.emplace(make_unique<bs1_t>(bs),
                                   std::make_pair(bit_range_end_ + 1,
                                                  bit_range_end_ + n_bits)
                                  );
    if(!r.second) throw basis_space_exists(bs);

    bit_range_end_ += n_bits;
    add_impl(std::forward<BSTypesTail>(more_bs)...);
  }
  void add_impl() {}

  // Check if all provided types are derived from generator<IndexTypes...>
  template<typename... Types>
  using all_are_basis_spaces = all_derived_from<basis_space_t, Types...>;

public:

  //
  // Exception types
  //

  struct basis_space_exists : public std::runtime_error {
    bs_ptr_type basis_space;
    basis_space_exists(basis_space_t const& bs) :
      std::runtime_error("Basis space already exists"),
      basis_space(bs.clone())
    {}
  };

  struct basis_space_not_found : public std::runtime_error {
    bs_ptr_type basis_space;
    basis_space_not_found(basis_space_t const& bs) :
      std::runtime_error("Basis space not found"),
      basis_space(bs.clone())
    {}
  };

  // Construct empty space
  hilbert_space() = default;

  // Construct from individual basis spaces
  template<typename... Args,
           typename = typename std::enable_if<
              all_are_basis_spaces<Args...>::value>::type
           >
  explicit hilbert_space(Args&&... args) {
    add_impl(std::forward<Args>(args)...);
  }

  // Collect all generators found in `expr` and add all respective basis spaces
  template<typename ScalarType>
  hilbert_space(expression<ScalarType, IndexTypes...> const& expr) {
    // TODO
  }

  // Copy all basis spaces from `hs` and add those found in `expr` and missing
  // from `hs`.
  template<typename ScalarType>
  hilbert_space(expression<ScalarType, IndexTypes...> const& expr,
                hilbert_space const& hs) {
    // TODO
  }

  // Value semantics
  hilbert_space(hilbert_space const& hs) : bit_range_end_(hs.bit_range_end_) {
    for(auto const& bs : hs.basis_spaces_) {
      basis_spaces_.emplace_hint(basis_spaces_.end(),
                                 bs.first->clone(),
                                 bs.second);
    }

  }
  hilbert_space(hilbert_space&&) noexcept = default;
  hilbert_space& operator=(hilbert_space const& hs) {
    bit_range_end_ = hs.bit_range_end_;
    basis_spaces_.clear();
    for(auto const& bs : hs.basis_spaces_) {
      auto r = basis_spaces_.emplace_hint(basis_spaces_.end(),
                                          bs.first->clone(),
                                          bs.second);
    }
    return *this;
  }
  hilbert_space& operator=(hilbert_space&&) noexcept = default;

  // Equality
  inline friend bool operator==(hilbert_space const& hs1,
                                hilbert_space const& hs2) {
    using value_type = typename decltype(basis_spaces_)::value_type;
    return hs1.size() == hs2.size() &&
           std::equal(hs1.basis_spaces_.begin(),
                      hs1.basis_spaces_.end(),
                      hs2.basis_spaces_.begin(),
                      [](value_type const& v1, value_type const& v2) {
                        return *v1.first == *v2.first && v1.second == v2.second;
                      }
          );
  }
  inline friend bool operator!=(hilbert_space const& hs1,
                                hilbert_space const& hs2) {
    return !operator==(hs1, hs2);
  }

  // Append a new basis space to the ordered product
  void add(basis_space_t const& bs) { add_impl(bs); }
  void add(basis_space_t && bs) { add_impl(std::move(bs));  }

  // Is a given basis space found in this Hilbert space
  bool has(basis_space_t const& bs) const {
    return basis_spaces_.count(bs.clone()) == 1;
  }

  // Bit range spanned by a basis space
  std::pair<int, int> bit_range(basis_space_t const& bs) const {
    auto it = basis_spaces_.find(bs.clone());
    if(it == basis_spaces_.end())
      throw basis_space_not_found(bs);
    else
      return it->second;
  }

  // Number of basis spaces
  size_t size() const { return basis_spaces_.size(); }

  // The minimal number of binary digits needed to represent any state
  // in this Hilbert space
  int total_n_bits() const { return bit_range_end_ + 1; }

private:

  // List of base spaces in the product and their corresponding bit ranges
  std::map<bs_ptr_type, std::pair<int, int>, less> basis_spaces_;

  // End of the bit range spanned by this Hilbert space
  int bit_range_end_ = -1;
};

}

#endif
