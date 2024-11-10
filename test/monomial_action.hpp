/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2024 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_TEST_MONOMIAL_ACTION_HPP_
#define LIBCOMMUTE_TEST_MONOMIAL_ACTION_HPP_

#include <libcommute/expression/monomial.hpp>
#include <libcommute/loperator/elementary_space.hpp>
#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/state_vector.hpp>

#include <memory>
#include <utility>

namespace libcommute {

//
// Padding elementary space: used to 'pad' a Hilbert space with unused bits
//   algebra_id() = fermion::algebra_id()-2
//   n_bits() = 2
//
template <typename... IndexTypes>
class elementary_space_padding final : public elementary_space<IndexTypes...> {
public:
  elementary_space_padding() = delete;
  template <typename... Args>
  explicit elementary_space_padding(Args&&... indices)
    : elementary_space<IndexTypes...>(std::forward<Args>(indices)...) {}

  // cppcheck-suppress duplInheritedMember
  std::unique_ptr<elementary_space<IndexTypes...>> clone() const final {
    return make_unique<elementary_space_padding>(*this);
  }
  int algebra_id() const final { return fermion - 2; }
  int n_bits() const final { return 2; }
};

//
// Reference action of a monomial on a basis state with a given index
//
template <typename RefGenAction, typename... IndexTypes>
bool ref_monomial_action(monomial<IndexTypes...> const& mon,
                         RefGenAction const& ga,
                         sv_index_type& index,
                         double& coeff) {
  for(auto it = mon.rbegin(); it != mon.rend(); ++it) {
    bool nonzero = ga(*it, index, coeff);
    if(!nonzero) return false;
  }
  return true;
}

//
// Check action of a monomial
//
template <typename MAType, typename RefGenAction, typename... IndexTypes>
void check_monomial_action(monomial<IndexTypes...> const& mon,
                           hilbert_space<IndexTypes...> const& hs,
                           RefGenAction const& ga,
                           std::vector<sv_index_type> const& in_index_list) {
  MAType ma(std::make_pair(mon.begin(), mon.end()), hs);

  for(auto in_index : in_index_list) {
    double coeff = 2;
    sv_index_type out_index = in_index;
    bool nonzero = ma.act(out_index, coeff);

    // Reference
    double coeff_ref = 2;
    sv_index_type out_index_ref = in_index;
    bool nonzero_ref = ref_monomial_action(mon, ga, out_index_ref, coeff_ref);
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
