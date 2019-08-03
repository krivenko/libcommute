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
#include "../metafunctions.hpp"
#include "../algebra_ids.hpp"
#include "../expression/generator_spin.hpp"

#include <memory>
#include <sstream>
#include <stdexcept>

namespace libcommute {

//
// Functor to construct basis spaces associated with algebra generators.
// This CRTP-base class only knows how to construct spaces for fermions and
// spins. It calls the `construct()` method of the derived class when another
// algebra ID is met.
//

template<typename Derived> class bs_constructor {

public:

  // Cannot construct a basis space
  template<typename... IndexTypes>
  struct construction_failure : public std::runtime_error {
    std::unique_ptr<generator<IndexTypes...>> generator_ptr;
    inline static std::string make_what(generator<IndexTypes...> const& g) {
      std::stringstream ss;
      ss << "Cannot construct basis space associated with generator " << g;
      return ss.str();
    }
    construction_failure(generator<IndexTypes...> const& g) :
      std::runtime_error(make_what(g)),
      generator_ptr(g.clone())
    {}
  };

  template<typename... IndexTypes>
  std::unique_ptr<basis_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    switch(g.algebra_id()) {
      case FERMION_ALGEBRA_ID:
        return make_unique<basis_space_fermion<IndexTypes...>>(g.indices());
      case SPIN_ALGEBRA_ID: {
        double spin =
          dynamic_cast<generator_spin<IndexTypes...> const&>(g).spin();
        return make_unique<basis_space_spin<IndexTypes...>>(spin, g.indices());
      };
      default:
        return static_cast<Derived const&>(*this).construct(g);
    }
  }
};

//
// Always fails for unknown algebra IDs (including BOSON_ALGEBRA_ID)
//
class default_bs_constructor : public bs_constructor<default_bs_constructor> {
public:
  template<typename... IndexTypes>
  std::unique_ptr<basis_space<IndexTypes...>>
  construct(generator<IndexTypes...> const& g) const {
    throw construction_failure<IndexTypes...>(g);
  }
};

//
// Constructs all bosonic spaces with the same dimension.
//
class boson_bs_constructor : public bs_constructor<boson_bs_constructor> {

  int bits_per_boson_;

public:

  boson_bs_constructor(int bits_per_boson) : bits_per_boson_(bits_per_boson) {}

  template<typename... IndexTypes>
  std::unique_ptr<basis_space<IndexTypes...>>
  construct(generator<IndexTypes...> const& g) const {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    if(g.algebra_id() == BOSON_ALGEBRA_ID) {
      return make_unique<basis_space_boson<IndexTypes...>>(bits_per_boson_,
                                                           g.indices());
    } else {
      throw construction_failure<IndexTypes...>(g);
    }
  }
};

} // namespace libcommute

#endif
