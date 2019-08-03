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
// Functors to construct basis spaces associated with algebra generators.
//

// Exception: Cannot construct a basis space
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
// This basis space constructor can handle only fermionic and spin generators
//
struct default_bs_constructor {
  template<typename... IndexTypes> std::unique_ptr<basis_space<IndexTypes...>>
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
        throw bs_construction_failure<IndexTypes...>(g);
    }
  }
};

//
// Constructs all bosonic spaces with equal dimensions
//
class bs_constructor_boson : public default_bs_constructor {

  int bits_per_boson_;

public:

  bs_constructor_boson(int bits_per_boson) : bits_per_boson_(bits_per_boson) {}

  template<typename... IndexTypes>
  std::unique_ptr<basis_space<IndexTypes...>>
  operator()(generator<IndexTypes...> const& g) const {
    if(g.algebra_id() == BOSON_ALGEBRA_ID) {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
      using std::make_unique;
#endif
      return make_unique<basis_space_boson<IndexTypes...>>(bits_per_boson_,
                                                           g.indices());
    } else
      return default_bs_constructor::operator()(g);
  }
};

} // namespace libcommute

#endif
