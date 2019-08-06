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
#ifndef LIBCOMMUTE_ALGEBRA_TAGS_HPP_
#define LIBCOMMUTE_ALGEBRA_TAGS_HPP_

#include <type_traits>

namespace libcommute {

//
// Predefined algebra tag types and respective algebra IDs
//

struct fermion { static constexpr int algebra_id() { return -3; } };
struct boson { static constexpr int algebra_id() { return -2; } };
struct spin { static constexpr int algebra_id() { return -1; } };

// The lowest algebra ID available to user-defined algebras
#define LIBCOMMUTE_MIN_USER_DEFINED_ALGEBRA_ID 0

//
// Check that a list of algebra tags is ordered according to their IDs
//

namespace detail {

template<typename Tag1, typename Tag2, typename... TagsTail>
struct algebra_tags_ordered_impl {
  static constexpr bool value = (Tag1::algebra_id() < Tag2::algebra_id()) &&
                            algebra_tags_ordered_impl<Tag2, TagsTail...>::value;
};
template<typename Tag1, typename Tag2>
struct algebra_tags_ordered_impl<Tag1, Tag2> {
  static constexpr bool value = Tag1::algebra_id() < Tag2::algebra_id();
};

} // namespace libcommute::detail

template<typename... AlgebraTags>
struct algebra_tags_ordered : detail::algebra_tags_ordered_impl<AlgebraTags...>
{};
// List of 1 tag is always ordered
template<typename AlgebraTag>
struct algebra_tags_ordered<AlgebraTag> : std::true_type {};
// List of zero tags is always ordered
template<> struct algebra_tags_ordered<> : std::true_type {};

} // namespace libcommute

#endif
