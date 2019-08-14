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
#ifndef LIBCOMMUTE_BASIS_SPACE_BOSON_HPP_
#define LIBCOMMUTE_BASIS_SPACE_BOSON_HPP_

#include "../algebra_tags.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"

#include "basis_space.hpp"

namespace libcommute {

//
// 2^n-dimensional basis space generated by one bosonic degree of freedom
//

template<typename... IndexTypes>
class basis_space_boson : public basis_space<IndexTypes...> {

  using base = basis_space<IndexTypes...>;

public:

  // Value symantics
  basis_space_boson() = delete;
  template<typename... Args>
  basis_space_boson(int n_bits, Args&&... indices)
    : base(std::forward<Args>(indices)...), n_bits_(n_bits) {
  }
  basis_space_boson(basis_space_boson const&) = default;
  basis_space_boson(basis_space_boson&&) noexcept = default;
  basis_space_boson& operator=(basis_space_boson const&) = default;
  basis_space_boson& operator=(basis_space_boson&&) noexcept = default;

  // Make a smart pointer that manages a copy of this basis space
  virtual std::unique_ptr<base> clone() const override {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    return make_unique<basis_space_boson>(*this);
  }

  // ID of the algebra this basis space is associated with
  virtual int algebra_id() const override { return boson::algebra_id(); }

  // The minimal number of binary digits needed to represent any state
  // in this basis space
  virtual int n_bits() const override { return n_bits_; }

protected:

  // This space is spanned by bosonic states |0>, |1>, ..., |2^{n_bits}-1>
  int n_bits_;
};

namespace static_indices {

// Convenience factory function
template<typename... IndexTypes>
inline basis_space_boson<c_str_to_string_t<IndexTypes>...>
make_space_boson(int n_bits, IndexTypes&&... indices) {
  return {n_bits, std::forward<IndexTypes>(indices)...};
}

} // namespace libcommute::static_indices
} // namespace libcommute

#if __cplusplus >= 201703L
#include "../expression/dyn_indices.hpp"

namespace libcommute {
namespace dynamic_indices {

// Convenience factory function for dynamic indices
template<typename... IndexTypes>
inline basis_space_boson<dyn_indices>
make_space_boson(int n_bits, IndexTypes&&... indices) {
  assert(n_bits > 0);
  return {n_bits, dyn_indices(std::forward<IndexTypes>(indices)...)};
}

} // namespace libcommute::dynamic_indices
} // namespace libcommute
#endif

#endif
