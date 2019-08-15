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
#include "../expression/generator.hpp"
#include "../expression/monomial.hpp"
#include "../utility.hpp"

#include <memory>
#include <sstream>
#include <stdexcept>

//
// Action of a monomial on a basis state
//

// In the present implementation we assume that all algebra generators
// are represented by generalized permutation matrices [1], i.e. matrices
// with exactly one non-zero element in each row and each column.
//
// [1] https://en.wikipedia.org/wiki/Generalized_permutation_matrix

namespace libcommute {

//
// Exception: Cannot construct a basis space
//
template<typename... IndexTypes>
struct unknown_generator : public std::runtime_error {
  std::unique_ptr<generator<IndexTypes...>> generator_ptr;
  inline static std::string make_what(generator<IndexTypes...> const& g) {
    std::stringstream ss;
    ss << "Action of generator " << g << " on a state vector is undefined";
    return ss.str();
  }
  unknown_generator(generator<IndexTypes...> const& g) :
    std::runtime_error(make_what(g)),
    generator_ptr(g.clone())
  {}
};


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
    if(m_range.first == m_range.second)
      return m_range;

    auto res_range = std::make_pair(m_range.first, m_range.first);

    for(auto it = m_range.first; it != m_range.second; ++it) {
      if(it->algebra_id() < AlgebraTag1::algebra_id())
        throw unknown_generator<IndexTypes...>(*it);
      else if(it->algebra_id() == AlgebraTag1::algebra_id()) {
        ++res_range.second;
        ++m_range.first;
      } else // it->algebra_id() > AlgebraTag1::algebra_id()
        break;
    }

    return res_range;
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
                  double & coeff) const {
    return base_head::act(in_index, out_index, coeff) &&
           base_tail::act(sv_index_type(out_index), out_index, coeff);
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
                  double & coeff) const {
    return base::act(in_index, out_index, coeff);
  }
};

} // namespace libcommute::detail

template<typename... AlgebraTags>
class monomial_action : public detail::monomial_action_impl<AlgebraTags...> {

  static_assert(algebra_tags_ordered<AlgebraTags...>::value,
                "Algebra tags must be ordered according to their IDs");

  using base = detail::monomial_action_impl<AlgebraTags...>;

public:

  template<typename... IndexTypes>
  monomial_action(detail::monomial_range_t<IndexTypes...> const& m_range,
                  hilbert_space<IndexTypes...> const& hs)
    : base(m_range, hs) {
  }
};

// Specialization for the case when no algebra tags have been provided
template<> class monomial_action<> {

public:

  template<typename... IndexTypes>
  monomial_action(detail::monomial_range_t<IndexTypes...> const& m_range,
                  hilbert_space<IndexTypes...> const& hs) {
    // Without algebra tags, we support only the constant (empty) monomial
    if(m_range.first != m_range.second)
      throw unknown_generator<IndexTypes...>(*m_range.first);
  }

  inline bool act(sv_index_type in_index,
                  sv_index_type & out_index,
                  double & coeff) const {
    out_index = in_index;
    return true;
  }
};

} // namespace libcommute

#endif
