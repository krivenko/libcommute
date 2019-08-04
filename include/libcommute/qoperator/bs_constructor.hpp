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
#ifndef LIBCOMMUTE_BS_CONSTRUCTOR_HPP_
#define LIBCOMMUTE_BS_CONSTRUCTOR_HPP_

#include "basis_space.hpp"
#include "basis_space_fermion.hpp"
#include "basis_space_boson.hpp"
#include "basis_space_spin.hpp"
#include "../algebra_tags.hpp"
#include "../metafunctions.hpp"
#include "../expression/generator_spin.hpp"
#include "../utility.hpp"

#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace libcommute {

//
// Exception: Cannot construct a basis space
//
template<typename... IndexTypes>
struct bs_construction_failure : public std::runtime_error {
  std::unique_ptr<generator<IndexTypes...>> generator_ptr;
  inline static std::string make_what(generator<IndexTypes...> const& g) {
    std::stringstream ss;
    ss << "Cannot construct basis space associated with algebra generator "
       << g;
    return ss.str();
  }
  bs_construction_failure(generator<IndexTypes...> const& g) :
    std::runtime_error(make_what(g)),
    generator_ptr(g.clone())
  {}
};

//
// Functors to construct basis spaces associated with algebra generators.
//

template<typename... AlgebraTags> class bs_constructor;

namespace detail {

template<typename AlgebraTag1, typename... AlgebraTagsTail>
class bs_constructor_impl : public bs_constructor<AlgebraTag1>,
                            public bs_constructor_impl<AlgebraTagsTail...> {

  using base_head = bs_constructor<AlgebraTag1>;
  using base_tail = bs_constructor_impl<AlgebraTagsTail...>;

public:

  bs_constructor_impl() = default;

  template<typename Arg1,
           typename... ArgsTail,
           typename std::enable_if<
             std::is_constructible<base_head, Arg1>::value,
             void*
           >::type = nullptr>
  inline bs_constructor_impl(Arg1&& arg1, ArgsTail&&... args_tail)
    : base_head(std::forward<Arg1>(arg1)),
      base_tail(std::forward<ArgsTail>(args_tail)...) {}

  template<typename Arg1,
           typename... ArgsTail,
           typename std::enable_if<
             !std::is_constructible<base_head, Arg1>::value,
             void*
           >::type = nullptr>
  inline bs_constructor_impl(Arg1&& arg1, ArgsTail&&... args_tail)
    : base_tail(std::forward<Arg1>(arg1),
                std::forward<ArgsTail>(args_tail)...) {}

  template<typename... IndexTypes>
  inline  std::unique_ptr<basis_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
    if(g.algebra_id() == AlgebraTag1::algebra_id())
      return base_head::operator()(g);
    else
      return base_tail::operator()(g);
  }
};

// Specialization of bs_constructor_impl: end of algebra tag chain
template<typename AlgebraTag>
class bs_constructor_impl<AlgebraTag> : public bs_constructor<AlgebraTag> {

  using base = bs_constructor<AlgebraTag>;

public:

  bs_constructor_impl() = default;
  template<typename Arg>
  inline bs_constructor_impl(Arg&& arg) : base(std::forward<Arg>(arg)) {}

  template<typename... IndexTypes>
  inline  std::unique_ptr<basis_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
    if(g.algebra_id() == AlgebraTag::algebra_id())
      return base::operator()(g);
    else
      throw bs_construction_failure<IndexTypes...>(g);
  }
};

} // namespace libcommute::detail

template<typename... AlgebraTags>
class bs_constructor : public detail::bs_constructor_impl<AlgebraTags...> {

  static_assert(all_types_different<AlgebraTags...>::value,
                "All algebra tags have to be different");

  using base = detail::bs_constructor_impl<AlgebraTags...>;

public:

  template<typename... Args>
  inline bs_constructor(Args&&... args) : base(std::forward<Args>(args)...) {}
};

//
// Basis space constructors for specific algebras
//

template<> class bs_constructor<fermion> {
public:
  bs_constructor() = default;

  template<typename... IndexTypes>
  inline std::unique_ptr<basis_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    return make_unique<basis_space_fermion<IndexTypes...>>(g.indices());
  }
};

template<> class bs_constructor<boson> {

  int bits_per_boson_;

public:
  bs_constructor() = delete;
  bs_constructor(int bits_per_boson) : bits_per_boson_(bits_per_boson) {}

  template<typename... IndexTypes>
  inline std::unique_ptr<basis_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    return make_unique<basis_space_boson<IndexTypes...>>(bits_per_boson_,
                                                         g.indices());
  }
};

template<> class bs_constructor<spin> {
public:
  bs_constructor() = default;

  template<typename... IndexTypes>
  inline std::unique_ptr<basis_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    double spin = dynamic_cast<generator_spin<IndexTypes...> const&>(g).spin();
    return make_unique<basis_space_spin<IndexTypes...>>(spin, g.indices());
  }
};

//
// Useful type aliases
//

using default_bs_constructor = bs_constructor<fermion, spin>;
using boson_bs_constructor = bs_constructor<fermion, boson, spin>;

} // namespace libcommute

#endif
