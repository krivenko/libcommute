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
#ifndef LIBCOMMUTE_LOPERATOR_LOPERATOR_HPP_
#define LIBCOMMUTE_LOPERATOR_LOPERATOR_HPP_

#include "../expression/expression.hpp"
#include "../metafunctions.hpp"
#include "../scalar_traits.hpp"
#include "../utility.hpp"
#include "hilbert_space.hpp"
#include "monomial_action_boson.hpp"
#include "monomial_action_fermion.hpp"
#include "monomial_action_spin.hpp"
#include "state_vector.hpp"

#include <utility>
#include <vector>

namespace libcommute {

//
// Linear operator acting on a state vector in a Hilbert space
//

template <typename ScalarType, int... AlgebraIDs> class loperator_base {

  using monomial_action_t = monomial_action<AlgebraIDs...>;

public:
  using scalar_type = ScalarType;

  loperator_base() = default;

  template <typename... IndexTypes>
  loperator_base(expression<scalar_type, IndexTypes...> const& expr,
                 hilbert_space<IndexTypes...> const& hs) {
    for(auto const& m : expr) {
      // cppcheck-suppress useStlAlgorithm
      m_actions_.emplace_back(
          monomial_action_t(
              std::make_pair(m.monomial.begin(), m.monomial.end()),
              hs),
          m.coeff);
    }
  }

  // Add monomial to the list of actions
  void add_monomial_action(monomial_action_t const& ma,
                           scalar_type const& coeff) {
    m_actions_.emplace_back(ma, coeff);
  }

  // Value semantics
  loperator_base(loperator_base const&) = default;
  loperator_base(loperator_base&&) noexcept = default;
  loperator_base& operator=(loperator_base const&) = default;
  loperator_base& operator=(loperator_base&&) noexcept = default;
  ~loperator_base() = default;

private:
  std::vector<std::pair<monomial_action_t, scalar_type>> m_actions_;

protected:
  inline std::vector<std::pair<monomial_action_t, scalar_type>> const&
  m_actions() const {
    return m_actions_;
  }
};

//
// Linear operator with constant monomial coefficients
//

template <typename ScalarType, int... AlgebraIDs>
class loperator : public loperator_base<ScalarType, AlgebraIDs...> {

  using base = loperator_base<ScalarType, AlgebraIDs...>;

public:
  loperator() = default;

  using base::base;

  // Value semantics
  loperator(loperator const&) = default;
  loperator(loperator&&) noexcept = default;
  loperator& operator=(loperator const&) = default;
  loperator& operator=(loperator&&) noexcept = default;
  ~loperator() = default;

  // Act on state and return the resulting state.
  template <typename StateVector>
  inline StateVector operator()(StateVector const& sv) const {
    StateVector res = zeros_like(sv);
    act_impl(sv, res);
    return res;
  }

  // Act on state `src` and return the resulting state via `dst`.
  template <typename SrcStateVector, typename DstStateVector>
  inline void operator()(SrcStateVector&& src, DstStateVector&& dst) const {
    set_zeros(dst);
    act_impl(std::forward<SrcStateVector>(src),
             std::forward<DstStateVector>(dst));
  }

  // Act on state and return the resulting state.
  template <typename StateVector>
  inline StateVector operator*(StateVector const& sv) const {
    return operator()(sv);
  }

private:
  // Implementation details of operator()
  template <typename SrcStateVector, typename DstStateVector>
  // NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
  inline void act_impl(SrcStateVector&& src, DstStateVector&& dst) const {

    // A workaround for GCC Bug 58972
    // (base::m_actions() would be inaccessible from within the lambda below)
    auto const& m_act = base::m_actions();

    foreach(std::forward<SrcStateVector>(src),
            [&](sv_index_type in_index,
                element_type_t<remove_cvref_t<SrcStateVector>> const& a) {
              for(auto const& ma : m_act) {
                sv_index_type index = in_index;
                auto coeff = scalar_traits<ScalarType>::make_const(1);
                bool nonzero = ma.first.act(index, coeff);
                if(nonzero)
                  update_add_element(dst, index, ma.second * coeff * a);
              }
            });
  }
};

//
// Linear operator with parameter-dependent (callable) coefficients
//

template <typename ScalarType, int... AlgebraIDs>
class parametric_loperator : public loperator_base<ScalarType, AlgebraIDs...> {

  using base = loperator_base<ScalarType, AlgebraIDs...>;

public:
  parametric_loperator() = delete;

  using base::base;

  // Value semantics
  parametric_loperator(parametric_loperator const&) = default;
  parametric_loperator(parametric_loperator&&) noexcept = default;
  parametric_loperator& operator=(parametric_loperator const&) = default;
  parametric_loperator& operator=(parametric_loperator&&) noexcept = default;
  ~parametric_loperator() = default;

  // Act on state and return the resulting state.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `args` as arguments to produce
  // the actual coefficient values.
  template <typename StateVector, typename... CoeffArgs>
  inline StateVector operator()(StateVector const& sv,
                                CoeffArgs&&... args) const {
    StateVector res = zeros_like(sv);
    std::vector<evaluated_coeff_t<CoeffArgs...>> evaluated_coeffs;
    evaluated_coeffs.reserve(base::m_actions().size());
    act_impl(sv, res, evaluated_coeffs, std::forward<CoeffArgs>(args)...);
    return res;
  }

  template <typename... CoeffArgs>
  using evaluated_coeff_t = invoke_result_t<ScalarType, CoeffArgs...>;

  // Act on state `src` and return the resulting state via `dst`.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `coeff_args` as arguments to produce
  // the actual coefficient values.
  template <typename StateVector, typename... CoeffArgs>
  inline void operator()(StateVector const& src,
                         StateVector& dst,
                         CoeffArgs&&... args) const {
    set_zeros(dst);
    std::vector<evaluated_coeff_t<CoeffArgs...>> evaluated_coeffs;
    evaluated_coeffs.reserve(base::m_actions().size());
    act_impl(src, dst, evaluated_coeffs, std::forward<CoeffArgs>(args)...);
  }

  // Act on state `src` and return the resulting state via `dst`.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `coeff_args` as arguments to produce
  // the actual coefficient values. This method does not allocate memory
  // to store those values and stores them in `evaluated_coeffs` instead.
  template <typename StateVector, typename... CoeffArgs>
  inline void act_and_store_coeffs(
      StateVector const& src,
      StateVector& dst,
      std::vector<evaluated_coeff_t<CoeffArgs...>>& evaluated_coeffs,
      CoeffArgs&&... args) const {

    set_zeros(dst);
    evaluated_coeffs.clear();
    act_impl(src, dst, evaluated_coeffs, std::forward<CoeffArgs>(args)...);
  }

  // Transform this parametric operator into the non-parametric form
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `args` as arguments to produce
  // the output coefficient values.
  template <typename... CoeffArgs>
  inline loperator<evaluated_coeff_t<CoeffArgs...>, AlgebraIDs...>
  at(CoeffArgs&&... args) const {
    loperator<evaluated_coeff_t<CoeffArgs...>, AlgebraIDs...> lop;
    for(auto const& m : base::m_actions()) {
      lop.add_monomial_action(m.first,
                              m.second(std::forward<CoeffArgs>(args)...));
    }
    return lop;
  }

private:
  // Implementation details of operator()
  template <typename StateVector, typename... CoeffArgs>
  inline void
  act_impl(StateVector const& src,
           StateVector& dst,
           std::vector<evaluated_coeff_t<CoeffArgs...>>& evaluated_coeffs,
           CoeffArgs&&... args) const {

    auto const& m_act = base::m_actions();

    // Evaluate coefficients
    evaluated_coeffs.reserve(m_act.size());
    for(auto const& a : m_act)
      // cppcheck-suppress useStlAlgorithm
      evaluated_coeffs.emplace_back(a.second(std::forward<CoeffArgs>(args)...));

    // Apply monomials
    foreach(
        src,
        [&](sv_index_type in_index, element_type_t<StateVector> const& a) {
          for(std::size_t n = 0; n < m_act.size(); ++n) {
            sv_index_type index = in_index;
            auto coeff =
                scalar_traits<evaluated_coeff_t<CoeffArgs...>>::make_const(1);
            bool nz = m_act[n].first.act(index, coeff);
            if(nz) {
              update_add_element(dst, index, evaluated_coeffs[n] * coeff * a);
            }
          }
        });
  }
};

// Factory function for loperator
template <typename ScalarType, typename... IndexTypes>
inline loperator<ScalarType, fermion, boson, spin>
make_loperator(expression<ScalarType, IndexTypes...> const& expr,
               hilbert_space<IndexTypes...> const& hs) {
  return loperator<ScalarType, fermion, boson, spin>(expr, hs);
}

// Factory function for parametric_loperator
template <typename ScalarType, typename... IndexTypes>
inline parametric_loperator<ScalarType, fermion, boson, spin>
make_param_loperator(expression<ScalarType, IndexTypes...> const& expr,
                     hilbert_space<IndexTypes...> const& hs) {
  return parametric_loperator<ScalarType, fermion, boson, spin>(expr, hs);
}

} // namespace libcommute

#endif
