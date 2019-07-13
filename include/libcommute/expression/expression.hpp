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
#ifndef LIBCOMMUTE_EXPRESSION_HPP_
#define LIBCOMMUTE_EXPRESSION_HPP_

#include "generator.hpp"
#include "monomial.hpp"
#include "scalar_traits.hpp"

#include <complex>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

namespace libcommute {

//
// Polynomial expression involving quantum-mechanical operators
//

template <typename ScalarType, typename... IndexTypes>
class expression {
public:

  using scalar_type = ScalarType;
  using index_types = std::tuple<IndexTypes...>;
  using monomial_t = monomial<IndexTypes...>;

private:
  using monomials_map_t = std::map<monomial_t, ScalarType>;

  // List of all monomials in this polynomial expression
  monomials_map_t monomials;

  // Expression with only the IndexTypes fixed
  template<typename S> using expression_t = expression<S, IndexTypes...>;

public:

  // Value semantics
  expression() = default;
  expression(expression const&) = default;
  expression(expression&&) noexcept = default;
  expression& operator=(expression const&) = default;
  expression& operator=(expression&&) noexcept = default;

  // Construct from an expression of a different scalar type
  template<typename S>
  expression(expression<S, IndexTypes...> const& x) {
   static_assert(std::is_constructible<scalar_type, S>::value,
                 "Incompatible scalar type in construction");
   *this = x;
  }

  // Construct expression with one constant term
  template<typename S>
  explicit expression(S const& x) {
    static_assert(std::is_constructible<scalar_type, S>::value,
                  "Incompatible scalar type in construction");
    if(!scalar_traits<S>::is_zero(x)) monomials.emplace(monomial_t{}, x);
  }

  // Construct from a monomial
  template<typename S>
  explicit expression(S const& x, monomial_t const& monomial) {
    static_assert(std::is_constructible<scalar_type, S>::value,
                  "Incompatible scalar type in construction");
    if(!scalar_traits<S>::is_zero(x)) monomials.emplace(monomial, x);
  }

  template<typename S>
  expression& operator=(expression<S, IndexTypes...> const& x) {
    static_assert(std::is_constructible<scalar_type, S>::value,
                  "Incompatible scalar type in assignment");
   monomials.clear();
   for (auto const& y : x.get_monomials())
     monomials.emplace(monomial_t(y.first), scalar_type(y.second));
   return *this;
  }

  //
  // Accessors
  //

  // Monomials container
  inline monomials_map_t const& get_monomials() const { return monomials; }
  inline monomials_map_t & get_monomials() { return monomials; }

  // Number of monomials
  inline size_t size() const { return monomials.size(); }

  //
  // Arithmetics
  //

  // Unary minus
  auto operator-() const -> expression_t<minus_type<ScalarType>> {
    return unary_minus_impl(std::is_same<minus_type<ScalarType>, ScalarType>());
  }

  // Multiplication by scalar (suffix form)
  template<typename S>
  auto operator*(S && alpha) const
    -> expression_t<mul_type<ScalarType, S>> {
    if(scalar_traits<typename std::remove_reference<S>::type>::is_zero(alpha))
      return {};
    else
      return mul_const_suffix_impl(
        std::forward<S>(alpha),
        std::is_same<mul_type<ScalarType, S>, ScalarType>()
      );
  }

  // Multiplication by scalar (prefix form)
  template<typename S>
  friend auto operator*(S && alpha, expression const& expr)
    -> expression_t<mul_type<S, ScalarType>> {
    if(scalar_traits<typename std::remove_reference<S>::type>::is_zero(alpha))
      return {};
    else
      return expr.mul_const_prefix_impl(
        std::forward<S>(alpha),
        std::is_same<mul_type<S, ScalarType>, ScalarType>()
      );
  }

  //
  // In-place arithmetics
  //

  // In-place multiplication by scalar
  template<typename S>
  expression & operator*=(S && alpha) {
    if(scalar_traits<typename std::remove_reference<S>::type>::is_zero(alpha))
      monomials.clear();
    else
      for(auto & p : monomials) p.second *= alpha;
    return *this;
  }

private:

  //
  // Implementations of arithmetic operations
  //

  // Unary minus: minus_type<ScalarType> == ScalarType
  inline expression_t<ScalarType> unary_minus_impl(std::true_type) const {
    expression res(*this);
    for(auto & p : res.monomials) p.second = -p.second;
    return res;
  }

  // Unary minus: minus_type<ScalarType> != ScalarType
  inline expression_t<minus_type<ScalarType>>
  unary_minus_impl(std::false_type) const {
    expression_t<minus_type<ScalarType>> res;
    auto & res_mons = res.get_monomials();
    for(auto const& p : get_monomials())
      res_mons.emplace_hint(res_mons.end(), p.first, -p.second);
    return res;
  }

  // Multiplication by scalar (suffix form)
  // mul_type<ScalarType, S> == ScalarType
  template<typename S>
  inline expression_t<mul_type<ScalarType, S>>
  mul_const_suffix_impl(S&& alpha, std::true_type) const {
    expression res(*this);
    for(auto & p : res.monomials) p.second = p.second * alpha;
    return res;
  }

  // Multiplication by scalar (suffix form)
  // mul_type<ScalarType, S> != ScalarType
  template<typename S>
  inline expression_t<mul_type<ScalarType, S>>
  mul_const_suffix_impl(S&& alpha, std::false_type) const {
    expression_t<mul_type<ScalarType, S>> res;
    auto & res_mons = res.get_monomials();
    for(auto const& p : get_monomials())
      res_mons.emplace_hint(res_mons.end(), p.first, p.second * alpha);
    return res;
  }

  // Multiplication by scalar (prefix form)
  // mul_type<ScalarType, S> == ScalarType
  template<typename S>
  inline expression_t<mul_type<S, ScalarType>>
  mul_const_prefix_impl(S&& alpha, std::true_type) const {
    expression res(*this);
    for(auto & p : res.monomials) p.second = alpha * p.second;
    return res;
  }

  // Multiplication by scalar (prefix form)
  // mul_type<ScalarType, S> != ScalarType
  template<typename S>
  inline expression_t<mul_type<S, ScalarType>>
  mul_const_prefix_impl(S&& alpha, std::false_type) const {
    expression_t<mul_type<S, ScalarType>> res;
    auto & res_mons = res.get_monomials();
    for(auto const& p : get_monomials())
      res_mons.emplace_hint(res_mons.end(), p.first, alpha * p.second);
    return res;
  }

  // TODO

};

// Aliases for specific scalar types
template<typename... IndexTypes>
using expression_real = expression<double, IndexTypes...>;
template<typename... IndexTypes>
using expression_complex = expression<std::complex<double>, IndexTypes...>;

} // namespace libcommute

#endif
