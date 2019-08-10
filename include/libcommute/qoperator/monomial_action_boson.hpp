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
#ifndef LIBCOMMUTE_MONOMIAL_ACTION_BOSON_HPP_
#define LIBCOMMUTE_MONOMIAL_ACTION_BOSON_HPP_

#include "hilbert_space.hpp"
#include "monomial_action.hpp"
#include "../expression/generator_boson.hpp"

//
// Action of a monomial comprised of bosonic algebra generators
//

namespace libcommute {

template<> class monomial_action<boson> {

  // TODO

public:

  template<typename... IndexTypes>
  monomial_action(detail::monomial_range_t<IndexTypes...> const& m_range,
                  hilbert_space<IndexTypes...> const& hs) {
    for(auto it = m_range.first; it != m_range.second; ++it) {
      if(it->algebra_id() != boson::algebra_id())
        throw unknown_generator<IndexTypes...>(*it);

      // TODO
    }

    // TODO
  }

  inline bool act(sv_index_type in_index,
                  sv_index_type & out_index,
                  double & coeff) {
    // TODO
    return true;
  }
};

} // namespace libcommute

#endif
