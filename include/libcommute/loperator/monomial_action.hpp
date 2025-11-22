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
#ifndef LIBCOMMUTE_LOPERATOR_MONOMIAL_ACTION_HPP_
#define LIBCOMMUTE_LOPERATOR_MONOMIAL_ACTION_HPP_

#include "../expression/generator.hpp"
#include "../expression/monomial.hpp"
#include "../utility.hpp"
#include "hilbert_space.hpp"
#include "state_vector.hpp"

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
// Exception: Cannot construct an elementary space
//
template <typename... IndexTypes>
struct unknown_generator : public std::runtime_error {
  std::shared_ptr<const generator<IndexTypes...>> generator_ptr;
  inline static std::string make_what(generator<IndexTypes...> const& g) {
    std::stringstream ss;
    ss << "Action of generator " << g << " on a state vector is undefined";
    return ss.str();
  }
  explicit unknown_generator(generator<IndexTypes...> const& g)
    : std::runtime_error(make_what(g)), generator_ptr(g.shared_from_this()) {}
};

template <int... AlgebraIDs> class monomial_action;

namespace detail {

template <typename... IndexTypes>
using monomial_range_t = typename monomial<IndexTypes...>::range_type;

template <int AlgebraID1, int... AlgebraIDsTail>
class monomial_action_impl : public monomial_action<AlgebraID1>,
                             public monomial_action_impl<AlgebraIDsTail...> {

  using base_head = monomial_action<AlgebraID1>;
  using base_tail = monomial_action_impl<AlgebraIDsTail...>;

  template <typename... IndexTypes>
  static monomial_range_t<IndexTypes...>
  find_algebra_monomial_subrange(monomial_range_t<IndexTypes...>& m_range) {
    if(m_range.first == m_range.second) return m_range;

    auto res_range = std::make_pair(m_range.first, m_range.first);

    for(auto it = m_range.first; it != m_range.second; ++it) {
      if(it->algebra_id() < AlgebraID1)
        throw unknown_generator<IndexTypes...>(*it);
      else if(it->algebra_id() == AlgebraID1) {
        ++res_range.second;
        ++m_range.first;
      } else // it->algebra_id() > AlgebraID1
        break;
    }

    return res_range;
  }

public:
  template <typename... IndexTypes>
  monomial_action_impl(monomial_range_t<IndexTypes...> m_range,
                       hilbert_space<IndexTypes...> const& hs)
    : base_head(find_algebra_monomial_subrange<IndexTypes...>(m_range), hs),
      base_tail(m_range, hs) {}

  template <typename ScalarType>
  // cppcheck-suppress duplInheritedMember
  inline bool act(sv_index_type& index, ScalarType& coeff) const {
    return base_head::act(index, coeff) && base_tail::act(index, coeff);
  }
};

// Specialization of monomial_action_impl: end of algebra ID chain
template <int AlgebraID>
class monomial_action_impl<AlgebraID> : public monomial_action<AlgebraID> {

  using base = monomial_action<AlgebraID>;

public:
  template <typename... IndexTypes>
  monomial_action_impl(monomial_range_t<IndexTypes...> const& m_range,
                       hilbert_space<IndexTypes...> const& hs)
    : base(m_range, hs) {}

  template <typename ScalarType>
  // cppcheck-suppress duplInheritedMember
  inline bool act(sv_index_type& index, ScalarType& coeff) const {
    return base::act(index, coeff);
  }
};

} // namespace detail

template <int... AlgebraIDs>
class monomial_action : public detail::monomial_action_impl<AlgebraIDs...> {

  static_assert(algebra_ids_ordered<AlgebraIDs...>::value,
                "Algebra IDs must be ordered according to their IDs");

  using base = detail::monomial_action_impl<AlgebraIDs...>;

public:
  template <typename... IndexTypes>
  monomial_action(detail::monomial_range_t<IndexTypes...> const& m_range,
                  hilbert_space<IndexTypes...> const& hs)
    : base(m_range, hs) {}
};

// Specialization for the case when no algebra IDs have been provided
template <> class monomial_action<> {

public:
  template <typename... IndexTypes>
  monomial_action(detail::monomial_range_t<IndexTypes...> const& m_range,
                  hilbert_space<IndexTypes...> const& hs) {
    // Without algebra IDs, we support only the constant (empty) monomial
    if(m_range.first != m_range.second)
      throw unknown_generator<IndexTypes...>(*m_range.first);
  }

  template <typename ScalarType>
  // cppcheck-suppress duplInheritedMember
  inline bool act(sv_index_type& index, ScalarType& coeff) const {
    return true;
  }
};

} // namespace libcommute

#endif
