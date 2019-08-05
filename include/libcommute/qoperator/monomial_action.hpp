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
#ifndef LIBCOMMUTE_MONOMIAL_ACTION_HPP_
#define LIBCOMMUTE_MONOMIAL_ACTION_HPP_

#include "hilbert_space.hpp"
#include "state_vector.hpp"
#include "../expression/monomial.hpp"
#include "../utility.hpp"

//
// Action of a monomial on a basis state
//

// In the present implementation we assume that all algebra generators
// are represented by generalized permutation matrices [1], i.e. matrices
// with exactly one non-zero element in each row and each column.
//
// [1] https://en.wikipedia.org/wiki/Generalized_permutation_matrix

namespace libcommute {

template<typename... AlgebraTags> class monomial_action;

namespace detail {

template<typename... IndexTypes>
using monomial_range_t = typename monomial<IndexTypes...>::range_type;

template<typename AlgebraTag1, typename... AlgebraTagsTail>
class monomial_action_impl : public monomial_action<AlgebraTag1>,
                             public monomial_action_impl<AlgebraTagsTail...> {

  using base_head = monomial_action<AlgebraTag1>;
  using base_tail = monomial_action_impl<AlgebraTagsTail...>;

  template<typename... IndexTypes>
  static monomial_range_t<IndexTypes...>
  find_algebra_monomial_subrange(monomial_range_t<IndexTypes...> & m_range) {
    // TODO
  }

public:

  template<typename... IndexTypes>
  monomial_action_impl(monomial_range_t<IndexTypes...> m_range,
                       hilbert_space<IndexTypes...> const& hs)
    : base_head(find_algebra_monomial_subrange<IndexTypes...>(m_range), hs),
      base_tail(m_range, hs) {
  }

  inline bool act(sv_index_type in_index,
                  sv_index_type & out_index,
                  double & coeff) {
    if(base_head::act(in_index, out_index, coeff))
      return base_tail::act(in_index, out_index, coeff);
  }
};

// Specialization of monomial_action_impl: end of algebra tag chain
template<typename AlgebraTag>
class monomial_action_impl<AlgebraTag> : public monomial_action<AlgebraTag> {

  using base = monomial_action<AlgebraTag>;

public:

  template<typename... IndexTypes>
  monomial_action_impl(monomial_range_t<IndexTypes...> const& m_range,
                       hilbert_space<IndexTypes...> const& hs)
    : base(m_range, hs) {
  }

  inline bool act(sv_index_type in_index,
                  sv_index_type & out_index,
                  double & coeff) {
    return base::act(in_index, out_index, coeff);
  }
};

} // namespace libcommute::detail

template<typename... AlgebraTags>
class monomial_action : public detail::monomial_action_impl<AlgebraTags...> {

  static_assert(all_types_different<AlgebraTags...>::value,
                "All algebra tags must be different");

  using base = detail::monomial_action_impl<AlgebraTags...>;

public:

  template<typename... IndexTypes>
  monomial_action(detail::monomial_range_t<IndexTypes...> const& m_range,
                  hilbert_space<IndexTypes...> const& hs)
    : base(m_range, hs) {
  }
};

} // namespace libcommute

#endif
