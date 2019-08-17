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
#ifndef LIBCOMMUTE_QOPERATOR_HPP_
#define LIBCOMMUTE_QOPERATOR_HPP_

#include "../expression/expression.hpp"
#include "hilbert_space.hpp"
#include "monomial_action_fermion.hpp"
#include "monomial_action_boson.hpp"
#include "monomial_action_spin.hpp"
#include "state_vector.hpp"
#include "../scalar_traits.hpp"
#include "../metafunctions.hpp"
#include "../utility.hpp"

#include <memory>
#include <utility>
#include <vector>

namespace libcommute {

//
// Quantum-mechanical operator acting on a state vector in a Hilbert space
//

template<typename ScalarType, typename... AlgebraTags>
class qoperator_base {

  using monomial_action_t = monomial_action<AlgebraTags...>;

public:

  using scalar_type = ScalarType;

  qoperator_base() = delete;

  template<typename... IndexTypes>
  qoperator_base(expression<scalar_type, IndexTypes...> const& expr,
                 hilbert_space<IndexTypes...> const& hs) {
    for(auto const& m : expr) {
      m_actions_.emplace_back(
        monomial_action_t(std::make_pair(m.monomial.begin(), m.monomial.end()),
                          hs),
        m.coeff
      );
    }
  }

  // Value semantics
  qoperator_base(qoperator_base const&) = default;
  qoperator_base(qoperator_base&&) noexcept = default;
  qoperator_base& operator=(qoperator_base const&) = default;
  qoperator_base& operator=(qoperator_base&&) noexcept = default;

protected:

  std::vector<std::pair<monomial_action_t, scalar_type>> m_actions_;
};

//
// Quantum-mechanical operator with constant monomial coefficients
//

template<typename ScalarType, typename... AlgebraTags>
class qoperator : public qoperator_base<ScalarType, AlgebraTags...> {

  using base = qoperator_base<ScalarType, AlgebraTags...>;

public:

  using base::base;

  // Value semantics
  qoperator(qoperator const&) = default;
  qoperator(qoperator&&) noexcept = default;
  qoperator& operator=(qoperator const&) = default;
  qoperator& operator=(qoperator&&) noexcept = default;

  // Act on state and return the resulting state.
  template<typename StateVector>
  inline StateVector operator()(StateVector const& sv) const {
    StateVector res = zeros_like(sv);
    act_impl(sv, res);
    return res;
  }

  // Act on state `src` and return the resulting state via `dst`.
  template<typename StateVector>
  inline void operator()(StateVector const& src, StateVector & dst) const {
    assert(get_size(src) == get_size(dst));
    set_zeros(dst);
    act_impl(src, dst);
  }

  // Act on state and return the resulting state.
  template<typename StateVector>
  inline StateVector operator*(StateVector const& sv) const {
    return operator()(sv);
  }

private:

  // Implementation details of operator()
  template<typename StateVector>
  inline void act_impl(StateVector const& src, StateVector & dst) const {
    foreach(src, [&,this](sv_index_type in_index,
                          element_type<StateVector> const& a) {
      for(auto const& ma : this->m_actions_) {
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
// Quantum-mechanical operator with parameter-dependent (callable) coefficients
//

template<typename ScalarType, typename... AlgebraTags>
class parametric_qoperator : public qoperator_base<ScalarType, AlgebraTags...> {

  using base = qoperator_base<ScalarType, AlgebraTags...>;

public:

  parametric_qoperator() = delete;

  template<typename... IndexTypes>
  parametric_qoperator(expression<ScalarType, IndexTypes...> const& expr,
                       hilbert_space<IndexTypes...> const& hs)
    : base(expr, hs), evaluated_coeffs_(base::m_actions_.size()) {
  }

  // Value semantics
  parametric_qoperator(parametric_qoperator const&) = default;
  parametric_qoperator(parametric_qoperator&&) noexcept = default;
  parametric_qoperator& operator=(parametric_qoperator const&) = default;
  parametric_qoperator& operator=(parametric_qoperator&&) noexcept = default;

  // Act on state and return the resulting state.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `args` as arguments to produce
  // the actual coefficient values.
  template<typename StateVector, typename... CoeffArgs>
  inline StateVector
  operator()(StateVector const& sv, CoeffArgs&&... args) const {
    StateVector res = zeros_like(sv);
    act_impl(sv, res, std::forward<CoeffArgs>(args)...);
    return res;
  }

  // Act on state `src` and return the resulting state via `dst`.
  //
  // Coefficients in front of monomials of the corresponding polynomial
  // expression are invoked with the `coeff_args` as arguments to produce
  // the actual coefficient values.
  template<typename StateVector, typename... CoeffArgs>
  inline void operator()(StateVector const& src,
                         StateVector & dst,
                         CoeffArgs&&... args) const {
    set_zeros(dst);
    act_impl(src, dst, std::forward<CoeffArgs>(args)...);
  }

private:

  // Type-erased preallocated container for evaluated values of monomial
  // coefficients
  mutable std::vector<std::shared_ptr<void>> evaluated_coeffs_;

  // Implementation details of operator()
  template<typename StateVector, typename... CoeffArgs>
  inline void act_impl(StateVector const& src,
                       StateVector & dst,
                       CoeffArgs&&... args) const {

#ifndef LIBCOMMUTE_NO_STD_INVOKE_RESULT
    using evaluated_coeff_t = std::invoke_result_t<ScalarType, CoeffArgs...>;
#else
    using evaluated_coeff_t = invoke_result_t<ScalarType, CoeffArgs...>;
#endif

    // Evaluate coefficients
    for(size_t n = 0; n < base::m_actions_.size(); ++n) {
      evaluated_coeffs_[n].reset(new evaluated_coeff_t(
        base::m_actions_[n].second(std::forward<CoeffArgs>(args)...))
      );
    }

    // Apply monomials
    foreach(src, [&,this](sv_index_type in_index,
                          element_type<StateVector> const& a) {
      for(size_t n = 0; n < base::m_actions_.size(); ++n) {
        sv_index_type index = in_index;
        auto coeff = scalar_traits<ScalarType>::make_const(1);
        bool nz = base::m_actions_[n].first.act(index, coeff);
        if(nz) {
          auto eval_coeff_ptr =
            static_cast<evaluated_coeff_t*>(evaluated_coeffs_[n].get());
          update_add_element(dst, index, (*eval_coeff_ptr) * coeff * a);
        }
      }
    });
  }
};

// Factory function for qoperator
template<typename ScalarType, typename... IndexTypes>
qoperator<ScalarType, fermion, boson, spin>
make_qoperator(expression<ScalarType, IndexTypes...> const& expr,
               hilbert_space<IndexTypes...> const& hs) {
  return qoperator<ScalarType, fermion, boson, spin>(expr, hs);
}

// Factory function for parametric_qoperator
template<typename ScalarType, typename... IndexTypes>
parametric_qoperator<ScalarType, fermion, boson, spin>
make_param_qoperator(expression<ScalarType, IndexTypes...> const& expr,
                     hilbert_space<IndexTypes...> const& hs) {
  return parametric_qoperator<ScalarType, fermion, boson, spin>(expr, hs);
}

} // namespace libcommute

#endif
