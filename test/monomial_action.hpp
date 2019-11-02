/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_TEST_MONOMIAL_ACTION_HPP_
#define LIBCOMMUTE_TEST_MONOMIAL_ACTION_HPP_

#include <libcommute/expression/monomial.hpp>
#include <libcommute/qoperator/basis_space.hpp>
#include <libcommute/qoperator/hilbert_space.hpp>
#include <libcommute/qoperator/state_vector.hpp>

#include <memory>
#include <utility>

namespace libcommute {

//
// Padding basis space: used to 'pad' a Hilbert space with unused bits
//   algebra_id() = fermion::algebra_id()-2
//   n_bits() = 2
//
template<typename... IndexTypes>
class basis_space_padding : public basis_space<IndexTypes...> {
public:
  basis_space_padding() = delete;
  template<typename... Args>
  basis_space_padding(Args&&... indices) :
    basis_space<IndexTypes...>(std::forward<Args>(indices)...) {}

  virtual std::unique_ptr<basis_space<IndexTypes...>> clone() const override {
    return make_unique<basis_space_padding>(*this);
  }
  virtual int algebra_id() const override { return fermion::algebra_id()-2; }
  virtual int n_bits() const override { return 2; }
};

//
// Reference action of a monomial on a basis state with a given index
//
template<typename RefGenAction, typename... IndexTypes>
bool ref_monomial_action(monomial<IndexTypes...> const& mon,
                         RefGenAction const& ga,
                         sv_index_type & index,
                         double & coeff
                         ) {
  for(auto it = mon.rbegin(); it != mon.rend(); ++it) {
    bool nonzero = ga(*it, index, coeff);
    if(!nonzero) return false;
  }
  return true;
}

//
// Check action of a monomial
//
template<typename MAType, typename RefGenAction, typename... IndexTypes>
void check_monomial_action(monomial<IndexTypes...> const& mon,
                           hilbert_space<IndexTypes...> const& hs,
                           RefGenAction const& ga,
                           std::vector<sv_index_type> const& in_index_list
                          ) {
  MAType ma(std::make_pair(mon.begin(), mon.end()), hs);

  for(auto in_index : in_index_list) {
    double coeff = 2;
    sv_index_type out_index = in_index;
    bool nonzero = ma.act(out_index, coeff);

    // Reference
    double coeff_ref = 2;
    sv_index_type out_index_ref = in_index;
    bool nonzero_ref = ref_monomial_action(mon,
                                           ga,
                                           out_index_ref,
                                           coeff_ref);
    // Checks
    CHECK(nonzero == nonzero_ref);
    if(nonzero_ref != false) {
      CHECK(out_index == out_index_ref);
      CHECK_THAT(coeff, Catch::WithinAbs(coeff_ref, 1e-10));
    }
  }
}

} // namespace libcommute

#endif
