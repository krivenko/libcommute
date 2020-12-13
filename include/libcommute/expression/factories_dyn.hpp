/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#ifndef LIBCOMMUTE_EXPRESSION_DYN_FACTORIES_HPP_
#define LIBCOMMUTE_EXPRESSION_DYN_FACTORIES_HPP_

#if __cplusplus < 201703L
#error "This header file requires a C++-17 compliant compiler"
#endif

#include "dyn_indices.hpp"
#include "expression.hpp"
#include "generator_fermion.hpp"
#include "generator_boson.hpp"
#include "generator_spin.hpp"
#include "../scalar_traits.hpp"
#include "../utility.hpp"

#include <utility>

//
// Factory functions for `expression<ScalarType, dyn_indices>`
//

namespace libcommute {
namespace dynamic_indices {

#define INDICES std::forward<IndexTypes>(indices)...
#define DEFINE_FACTORY(NAME, ...)                                         \
template<typename ScalarType = double, typename... IndexTypes>            \
inline                                                                    \
expression<ScalarType, dyn_indices>                                       \
NAME(IndexTypes&&... indices) {                                           \
  using ret_t = expression<ScalarType, dyn_indices>;                      \
  return ret_t(scalar_traits<ScalarType>::make_const(1),                  \
               typename ret_t::monomial_t(__VA_ARGS__)                    \
  );                                                                      \
}

#define DEFINE_FACTORY_SPIN(NAME, ...)                                    \
template<int Multiplicity,                                                \
         typename ScalarType = double,                                    \
         typename... IndexTypes>                                          \
inline                                                                    \
expression<ScalarType, dyn_indices>                                       \
NAME(IndexTypes&&... indices) {                                           \
  static_assert(Multiplicity >= 2, "Invalid multiplicity");               \
  using ret_t = expression<ScalarType, dyn_indices>;                      \
  return ret_t(scalar_traits<ScalarType>::make_const(1),                  \
               typename ret_t::monomial_t(__VA_ARGS__)                    \
  );                                                                      \
}

//
// Free functions to make fermionic operators
//

// Creation operator
DEFINE_FACTORY(c_dag, (dynamic_indices::make_fermion(true, INDICES)))
// Annihilation operator
DEFINE_FACTORY(c, (dynamic_indices::make_fermion(false, INDICES)))
// Number of fermions
DEFINE_FACTORY(n, (dynamic_indices::make_fermion(true, indices...)),
                  (dynamic_indices::make_fermion(false, indices...)))

//
// Free functions to make bosonic operators
//

// Creation operator
DEFINE_FACTORY(a_dag, (dynamic_indices::make_boson(true, INDICES)))
// Annihilation operator
DEFINE_FACTORY(a, (dynamic_indices::make_boson(false, INDICES)))

//
// Free functions to make spin operators (spin 1/2)
//

// Raising operator
DEFINE_FACTORY(S_p, (dynamic_indices::make_spin(spin_component::plus, INDICES)))
// Lowering operator
DEFINE_FACTORY(S_m, (dynamic_indices::make_spin(
                                                spin_component::minus,
                                                INDICES
                                               )
                    )
              )
// S_z
DEFINE_FACTORY(S_z, (dynamic_indices::make_spin(spin_component::z, INDICES)))

//
// Free functions to make spin operators (arbitrary spin)
//
// Raising operator
DEFINE_FACTORY_SPIN(S_p, (dynamic_indices::make_spin(
                                                     (Multiplicity-1)/2.0,
                                                     spin_component::plus,
                                                     INDICES
                                                    )
                         )
                   )
// Lowering operator
DEFINE_FACTORY_SPIN(S_m, (dynamic_indices::make_spin(
                                                     (Multiplicity-1)/2.0,
                                                     spin_component::minus,
                                                     INDICES
                                                    )
                         )
                   )
// S_z
DEFINE_FACTORY_SPIN(S_z, (dynamic_indices::make_spin(
                                                     (Multiplicity-1)/2.0,
                                                     spin_component::z,
                                                     INDICES
                                                    )
                         )
                   )

//
// In the complex case, we can additionally define S_x and S_y
//

template<typename... IndexTypes>
inline expression<std::complex<double>, dyn_indices>
S_x(IndexTypes&&... indices) {
  return std::complex<double>(0.5) * (S_p(indices...) + S_m(indices...));
}
template<typename... IndexTypes>
inline expression<std::complex<double>, dyn_indices>
S_y(IndexTypes&&... indices) {
  return std::complex<double>(0, -0.5) * (S_p(indices...) - S_m(indices...));
}

template<int Multiplicity, typename... IndexTypes>
inline expression<std::complex<double>, dyn_indices>
S_x(IndexTypes&&... indices) {
  return std::complex<double>(0.5) * (S_p<Multiplicity>(indices...) +
                                      S_m<Multiplicity>(indices...));
}
template<int Multiplicity, typename... IndexTypes>
inline expression<std::complex<double>, dyn_indices>
S_y(IndexTypes&&... indices) {
  return std::complex<double>(0, -0.5) * (S_p<Multiplicity>(indices...) -
                                          S_m<Multiplicity>(indices...));
}

// Make a complex expression out of a real one
expr_complex make_complex(expr_real const& expr) {
  return expr_complex(expr);
}

#undef DEFINE_FACTORY
#undef DEFINE_FACTORY_SPIN
#undef INDICES

} // namespace libcommute::dynamic_indices
} // namespace libcommute

#endif
