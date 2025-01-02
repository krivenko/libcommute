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
#ifndef LIBCOMMUTE_ALGEBRA_IDS_HPP_
#define LIBCOMMUTE_ALGEBRA_IDS_HPP_

#include <type_traits>

namespace libcommute {

//
// Predefined algebra IDs
//

static constexpr int fermion = -3;
static constexpr int boson = -2;
static constexpr int spin = -1;

// The lowest algebra ID available to user-defined algebras
static constexpr int min_user_defined_algebra_id = 0;

//
// Check that a list of algebra IDs is ordered
//

namespace detail {

template <int ID1, int ID2, int... IDsTail> struct algebra_ids_ordered_impl {
  static constexpr bool value =
      (ID1 < ID2) && algebra_ids_ordered_impl<ID2, IDsTail...>::value;
};
template <int ID1, int ID2> struct algebra_ids_ordered_impl<ID1, ID2> {
  static constexpr bool value = ID1 < ID2;
};

} // namespace detail

template <int... AlgebraIDs>
struct algebra_ids_ordered : detail::algebra_ids_ordered_impl<AlgebraIDs...> {};
// List of 1 ID is always ordered
template <int AlgebraID>
struct algebra_ids_ordered<AlgebraID> : std::true_type {};
// List of zero IDs is always ordered
template <> struct algebra_ids_ordered<> : std::true_type {};

} // namespace libcommute

#endif
