/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_LOPERATOR_ES_CONSTRUCTOR_HPP_
#define LIBCOMMUTE_LOPERATOR_ES_CONSTRUCTOR_HPP_

#include "../algebra_ids.hpp"
#include "../expression/generator_spin.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"
#include "elementary_space.hpp"
#include "elementary_space_boson.hpp"
#include "elementary_space_fermion.hpp"
#include "elementary_space_spin.hpp"

#include <memory>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace libcommute {

//
// Exception: Cannot construct an elementary space
//
template <typename... IndexTypes>
struct es_construction_failure : public std::runtime_error {
  std::unique_ptr<generator<IndexTypes...>> generator_ptr;
  inline static std::string make_what(generator<IndexTypes...> const& g) {
    std::stringstream ss;
    ss << "Cannot construct elementary space associated with algebra generator "
       << g;
    return ss.str();
  }
  explicit es_construction_failure(generator<IndexTypes...> const& g)
    : std::runtime_error(make_what(g)), generator_ptr(g.clone()) {}
};

//
// Functors to construct elementary spaces associated with algebra generators.
//

template <int... AlgebraIDs> class es_constructor;

namespace detail {

template <int AlgebraID1, int... AlgebraIDsTail>
class es_constructor_impl : public es_constructor<AlgebraID1>,
                            public es_constructor_impl<AlgebraIDsTail...> {

  using base_head = es_constructor<AlgebraID1>;
  using base_tail = es_constructor_impl<AlgebraIDsTail...>;

public:
  es_constructor_impl() = default;

  template <
      typename Arg1,
      typename... ArgsTail,
      typename std::enable_if<std::is_constructible<base_head, Arg1>::value,
                              void*>::type = nullptr>
  inline es_constructor_impl(Arg1&& arg1, ArgsTail&&... args_tail)
    : base_head(std::forward<Arg1>(arg1)),
      base_tail(std::forward<ArgsTail>(args_tail)...) {}

  template <
      typename Arg1,
      typename... ArgsTail,
      typename std::enable_if<!std::is_constructible<base_head, Arg1>::value,
                              void*>::type = nullptr>
  inline es_constructor_impl(Arg1&& arg1, ArgsTail&&... args_tail)
    : base_tail(std::forward<Arg1>(arg1),
                std::forward<ArgsTail>(args_tail)...) {}

  template <typename... IndexTypes>
  inline std::unique_ptr<elementary_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
    if(g.algebra_id() == AlgebraID1)
      return base_head::operator()(g);
    else
      return base_tail::operator()(g);
  }
};

// Specialization of es_constructor_impl: end of algebra ID chain
template <int AlgebraID>
class es_constructor_impl<AlgebraID> : public es_constructor<AlgebraID> {

  using base = es_constructor<AlgebraID>;

public:
  es_constructor_impl() = default;
  template <typename Arg>
  inline explicit es_constructor_impl(Arg&& arg)
    : base(std::forward<Arg>(arg)) {}

  template <typename... IndexTypes>
  inline std::unique_ptr<elementary_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
    if(g.algebra_id() == AlgebraID)
      return base::operator()(g);
    else
      throw es_construction_failure<IndexTypes...>(g);
  }
};

} // namespace detail

template <int... AlgebraIDs>
class es_constructor : public detail::es_constructor_impl<AlgebraIDs...> {

  static_assert(sizeof...(AlgebraIDs) > 0,
                "There must be at least one algebra ID");
  static_assert(algebra_ids_ordered<AlgebraIDs...>::value,
                "Algebra IDs must be ordered");

  using base = detail::es_constructor_impl<AlgebraIDs...>;

public:
  template <typename... Args>
  // cppcheck-suppress noExplicitConstructor
  inline es_constructor(Args&&... args) : base(std::forward<Args>(args)...) {}
};

//
// Elementary space constructors for specific algebras
//

template <> class es_constructor<fermion> {
public:
  es_constructor() = default;

  template <typename... IndexTypes>
  inline std::unique_ptr<elementary_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
    return make_unique<elementary_space_fermion<IndexTypes...>>(g.indices());
  }
};

template <> class es_constructor<boson> {

  int bits_per_boson_;

public:
  es_constructor() = delete;
  explicit es_constructor(int bits_per_boson)
    : bits_per_boson_(bits_per_boson) {}

  template <typename... IndexTypes>
  inline std::unique_ptr<elementary_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
    return make_unique<elementary_space_boson<IndexTypes...>>(bits_per_boson_,
                                                              g.indices());
  }
};

template <> class es_constructor<spin> {
public:
  es_constructor() = default;

  template <typename... IndexTypes>
  inline std::unique_ptr<elementary_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
    double s = dynamic_cast<generator_spin<IndexTypes...> const&>(g).spin();
    return make_unique<elementary_space_spin<IndexTypes...>>(s, g.indices());
  }
};

//
// Useful type aliases
//

using default_es_constructor = es_constructor<fermion, spin>;
using boson_es_constructor = es_constructor<fermion, boson, spin>;

} // namespace libcommute

#endif
