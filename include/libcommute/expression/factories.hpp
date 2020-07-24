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
#ifndef LIBCOMMUTE_EXPRESSION_FACTORIES_HPP_
#define LIBCOMMUTE_EXPRESSION_FACTORIES_HPP_

#include "expression.hpp"
#include "generator_fermion.hpp"
#include "generator_boson.hpp"
#include "generator_spin.hpp"
#include "../scalar_traits.hpp"
#include "../utility.hpp"

#include <utility>

//
// Factory functions for `expression<ScalarType, IndexTypes...>`
//

namespace libcommute {
namespace static_indices {

#define INDICES std::forward<IndexTypes>(indices)...
#define DEFINE_FACTORY(NAME, ...)                                         \
template<typename ScalarType, typename... IndexTypes>                     \
inline                                                                    \
expression<ScalarType, c_str_to_string_t<IndexTypes>...>                  \
NAME(IndexTypes&&... indices) {                                           \
  using ret_t = expression<ScalarType, c_str_to_string_t<IndexTypes>...>; \
  return ret_t(scalar_traits<ScalarType>::make_const(1),                  \
               typename ret_t::monomial_t(__VA_ARGS__)                    \
  );                                                                      \
}

#define DEFINE_FACTORY_SPIN(NAME, ...)                                    \
template<int Multiplicity, typename ScalarType, typename... IndexTypes>   \
inline                                                                    \
expression<ScalarType, c_str_to_string_t<IndexTypes>...>                  \
NAME(IndexTypes&&... indices) {                                           \
  static_assert(Multiplicity >= 2, "Invalid multiplicity");               \
  using ret_t = expression<ScalarType, c_str_to_string_t<IndexTypes>...>; \
  return ret_t(scalar_traits<ScalarType>::make_const(1),                  \
               typename ret_t::monomial_t(__VA_ARGS__)                    \
  );                                                                      \
}

#define DEFINE_FACTORY_SCALAR_TYPE(NAME, S)                             \
template<typename... IndexTypes>                                        \
inline expression<S, c_str_to_string_t<IndexTypes>...>                  \
NAME(IndexTypes&&... indices) {                                         \
  using libcommute::static_indices::NAME;                               \
  return NAME<S>(INDICES);                                              \
}

#define DEFINE_FACTORY_SPIN_SCALAR_TYPE(NAME, S)                        \
template<int Multiplicity, typename... IndexTypes>                      \
inline expression<S, c_str_to_string_t<IndexTypes>...>                  \
NAME(IndexTypes&&... indices) {                                         \
  using libcommute::static_indices::NAME;                               \
  return NAME<Multiplicity, S>(INDICES);                                \
}

//
// Free functions to make fermionic operators
//

// Creation operator
DEFINE_FACTORY(c_dag, (make_fermion(true, INDICES)))
// Annihilation operator
DEFINE_FACTORY(c, (make_fermion(false, INDICES)))
// Number of fermions
DEFINE_FACTORY(n, (make_fermion(true, INDICES)), (make_fermion(false, INDICES)))

//
// Free functions to make bosonic operators
//

// Creation operator
DEFINE_FACTORY(a_dag, (make_boson(true, INDICES)))
// Annihilation operator
DEFINE_FACTORY(a, (make_boson(false, INDICES)))

//
// Free functions to make spin operators (spin 1/2)
//

// Raising operator
DEFINE_FACTORY(S_p, (make_spin(spin_component::plus, INDICES)))
// Lowering operator
DEFINE_FACTORY(S_m, (make_spin(spin_component::minus, INDICES)))
// S_z
DEFINE_FACTORY(S_z, (make_spin(spin_component::z, INDICES)))

//
// Free functions to make spin operators (arbitrary spin)
//
// Raising operator
DEFINE_FACTORY_SPIN(S_p, (make_spin((Multiplicity-1)/2.0,
                                    spin_component::plus, INDICES)))
// Lowering operator
DEFINE_FACTORY_SPIN(S_m, (make_spin((Multiplicity-1)/2.0,
                                    spin_component::minus,INDICES)))
// S_z
DEFINE_FACTORY_SPIN(S_z, (make_spin((Multiplicity-1)/2.0,
                                    spin_component::z,INDICES)))

//
// Specializations for ScalarType = double
//

namespace real {
DEFINE_FACTORY_SCALAR_TYPE(c_dag, double)
DEFINE_FACTORY_SCALAR_TYPE(c, double)
DEFINE_FACTORY_SCALAR_TYPE(n, double)

DEFINE_FACTORY_SCALAR_TYPE(a_dag, double)
DEFINE_FACTORY_SCALAR_TYPE(a, double)

DEFINE_FACTORY_SCALAR_TYPE(S_p, double)
DEFINE_FACTORY_SCALAR_TYPE(S_m, double)
DEFINE_FACTORY_SCALAR_TYPE(S_z, double)

DEFINE_FACTORY_SPIN_SCALAR_TYPE(S_p, double)
DEFINE_FACTORY_SPIN_SCALAR_TYPE(S_m, double)
DEFINE_FACTORY_SPIN_SCALAR_TYPE(S_z, double)
} // namespace libcommute::static_indices::real

//
// Specializations for ScalarType = std::complex<double>
//

namespace complex {
DEFINE_FACTORY_SCALAR_TYPE(c_dag, std::complex<double>)
DEFINE_FACTORY_SCALAR_TYPE(c, std::complex<double>)
DEFINE_FACTORY_SCALAR_TYPE(n, std::complex<double>)

DEFINE_FACTORY_SCALAR_TYPE(a_dag, std::complex<double>)
DEFINE_FACTORY_SCALAR_TYPE(a, std::complex<double>)

DEFINE_FACTORY_SCALAR_TYPE(S_p, std::complex<double>)
DEFINE_FACTORY_SCALAR_TYPE(S_m, std::complex<double>)
DEFINE_FACTORY_SCALAR_TYPE(S_z, std::complex<double>)

DEFINE_FACTORY_SPIN_SCALAR_TYPE(S_p, std::complex<double>)
DEFINE_FACTORY_SPIN_SCALAR_TYPE(S_m, std::complex<double>)
DEFINE_FACTORY_SPIN_SCALAR_TYPE(S_z, std::complex<double>)

//
// In the complex case, we can additionally define S_x and S_y
//

template<typename... IndexTypes>
inline expression<std::complex<double>, c_str_to_string_t<IndexTypes>...>
S_x(IndexTypes&&... indices) {
  return std::complex<double>(0.5) * (S_p(INDICES) + S_m(INDICES));
}
template<typename... IndexTypes>
inline expression<std::complex<double>, c_str_to_string_t<IndexTypes>...>
S_y(IndexTypes&&... indices) {
  return std::complex<double>(0, -0.5) * (S_p(INDICES) - S_m(INDICES));
}

template<int Multiplicity, typename... IndexTypes>
inline expression<std::complex<double>, c_str_to_string_t<IndexTypes>...>
S_x(IndexTypes&&... indices) {
  return std::complex<double>(0.5) * (S_p<Multiplicity>(INDICES) +
                                      S_m<Multiplicity>(INDICES));
}
template<int Multiplicity, typename... IndexTypes>
inline expression<std::complex<double>, c_str_to_string_t<IndexTypes>...>
S_y(IndexTypes&&... indices) {
  return std::complex<double>(0, -0.5) * (S_p<Multiplicity>(INDICES) -
                                          S_m<Multiplicity>(INDICES));
}

} // namespace libcommute::static_indices::complex

} // namespace libcommute::static_indices
} // namespace libcommute

#undef DEFINE_FACTORY
#undef DEFINE_FACTORY_SPIN
#undef DEFINE_FACTORY_SCALAR_TYPE
#undef DEFINE_FACTORY_SPIN_SCALAR_TYPE
#undef INDICES

#endif
