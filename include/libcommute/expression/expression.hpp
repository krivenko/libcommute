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
// Instances of ScalarType are assumed to form a field under addition and
// multiplication
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
  monomials_map_t monomials_;

  // Expression with only the IndexTypes fixed
  template<typename S> using expression_t = expression<S, IndexTypes...>;

  // Is T an instance of expression_t?
  template<typename T> struct is_expression : std::false_type {};
  template<typename S> struct is_expression<expression_t<S>>
    : std::true_type {};

  // Disable overload for expressions
  template<typename T> using disable_for_expression =
    typename std::enable_if<!is_expression<T>::value>::type;

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
    if(!scalar_traits<S>::is_zero(x)) monomials_.emplace(monomial_t{}, x);
  }

  // Construct from a monomial
  template<typename S>
  explicit expression(S const& x, monomial_t const& monomial) {
    static_assert(std::is_constructible<scalar_type, S>::value,
                  "Incompatible scalar type in construction");
    if(!scalar_traits<S>::is_zero(x)) monomials_.emplace(monomial, x);
  }

  template<typename S>
  expression& operator=(expression<S, IndexTypes...> const& x) {
    static_assert(std::is_constructible<scalar_type, S>::value,
                  "Incompatible scalar type in assignment");
   monomials_.clear();
   for (auto const& y : x.get_monomials())
     monomials_.emplace(monomial_t(y.first), scalar_type(y.second));
   return *this;
  }

  //
  // Accessors
  //

  // Monomials container
  inline monomials_map_t const& get_monomials() const { return monomials_; }
  inline monomials_map_t & get_monomials() { return monomials_; }

  // Number of monomials
  inline size_t size() const { return monomials_.size(); }

  // Equality
  inline friend bool operator==(expression const& e1, expression const& e2) {
    return e1.monomials_ == e2.monomials_;
  }
  inline friend bool operator!=(expression const& e1, expression const& e2) {
    return !operator==(e1, e2);
  }

  //
  // Arithmetics
  //

  template<typename S>
  expression_t<sum_type<ScalarType, S>>
  operator+(expression_t<S> const& expr) {
    return add_impl(expr, std::is_same<sum_type<ScalarType, S>, ScalarType>());
  }

  template<typename S>
  expression_t<diff_type<ScalarType, S>>
  operator-(expression_t<S> const& expr) {
    return sub_impl(expr, std::is_same<diff_type<ScalarType, S>, ScalarType>());
  }

  //
  // Compound assignments
  //

  template<typename S>
  expression & operator+=(expression_t<S> const& expr) {
    for(auto const& p : expr.get_monomials()) {
      auto it = monomials_.find(p.first);
      if(it == monomials_.end())
        monomials_.emplace(p);
      else {
        it->second = it->second + p.second;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          monomials_.erase(it);
      }
    }
    return *this;
  }

  template<typename S>
  expression & operator-=(expression_t<S> const& expr) {
    for(auto const& p : expr.get_monomials()) {
      auto it = monomials_.find(p.first);
      if(it == monomials_.end())
        monomials_.emplace(p.first, -p.second);
      else {
        it->second = it->second - p.second;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          monomials_.erase(it);
      }
    }
    return *this;
  }

  // Unary minus
  auto operator-() const -> expression_t<minus_type<ScalarType>> {
    return unary_minus_impl(std::is_same<minus_type<ScalarType>, ScalarType>());
  }

  //
  // Arithmetics involving constants
  //

  // Multiplication by scalar (postfix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  auto operator*(S && alpha) const
    -> expression_t<mul_type<ScalarType, S>> {
    if(scalar_traits<remove_cvref_t<S>>::is_zero(alpha))
      return {};
    else
      return mul_const_postfix_impl(
        std::forward<S>(alpha),
        std::is_same<mul_type<ScalarType, S>, ScalarType>()
      );
  }

  // Multiplication by scalar (prefix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  friend auto operator*(S && alpha, expression const& expr)
    -> expression_t<mul_type<S, ScalarType>> {
    if(scalar_traits<remove_cvref_t<S>>::is_zero(alpha))
      return {};
    else
      return expr.mul_const_prefix_impl(
        std::forward<S>(alpha),
        std::is_same<mul_type<S, ScalarType>, ScalarType>()
      );
  }

  // Addition of scalar (postfix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  auto operator+(S && alpha) const -> expression_t<sum_type<ScalarType, S>> {
    return add_const_postfix_impl(
      std::forward<S>(alpha),
      std::is_same<sum_type<ScalarType, S>, ScalarType>()
    );
  }

  // Addition of scalar (prefix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  friend auto operator+(S && alpha, expression const& expr)
    -> expression_t<sum_type<S, ScalarType>> {
    return expr.add_const_prefix_impl(
      std::forward<S>(alpha),
      std::is_same<sum_type<S, ScalarType>, ScalarType>()
    );
  }

  // Subtraction of scalar (postfix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  auto operator-(S && alpha) const -> expression_t<diff_type<ScalarType, S>> {
    return sub_const_postfix_impl(
      std::forward<S>(alpha),
      std::is_same<diff_type<ScalarType, S>, ScalarType>()
    );
  }

  // Subtraction of scalar (prefix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  friend auto operator-(S && alpha, expression const& expr)
    -> expression_t<diff_type<S, ScalarType>> {
    using res_s_t = diff_type<S, ScalarType>;
    expression_t<res_s_t> res;
    auto z = scalar_traits<S>::make_const(0);
    auto & res_mons = res.get_monomials();
    for(auto const& p : expr.monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, z - p.second);
    if(!scalar_traits<S>::is_zero(alpha)) {
      auto it = res_mons.find(monomial_t{});
      if(it == res_mons.end()) {
        auto z = scalar_traits<ScalarType>::make_const(0);
        res_mons.emplace_hint(res_mons.begin(), monomial_t{}, alpha - z);
      } else {
        it->second = alpha + it->second;
        if(scalar_traits<res_s_t>::is_zero(it->second))
          res_mons.erase(it);
      }
    }
    return res;
  }

  //
  // Compound assignments from constants
  //

  // Compound assignment/multiplication by scalar
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  expression & operator*=(S && alpha) {
    if(scalar_traits<remove_cvref_t<S>>::is_zero(alpha))
      monomials_.clear();
    else
      for(auto & p : monomials_) p.second *= alpha;
    return *this;
  }

  // Compound assignments/addition of scalar
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  expression & operator+=(S && alpha) {
    using s_t = remove_cvref_t<S>;
    if(!scalar_traits<s_t>::is_zero(alpha)) {
      auto it = monomials_.find(monomial_t{});
      if(it == monomials_.end()) {
        auto val = scalar_traits<ScalarType>::make_const(0);
        val += alpha;
        monomials_.emplace_hint(monomials_.begin(), monomial_t{}, val);
      } else {
        it->second += alpha;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          monomials_.erase(it);
      }
    }
    return *this;
  }

  // Compound assignments/subtraction of scalar
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  expression & operator-=(S && alpha) {
    using s_t = remove_cvref_t<S>;
    if(!scalar_traits<s_t>::is_zero(alpha)) {
      auto it = monomials_.find(monomial_t{});
      if(it == monomials_.end()) {
        monomials_.emplace_hint(monomials_.begin(), monomial_t{}, -alpha);
      } else {
        it->second -= alpha;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          monomials_.erase(it);
      }
    }
    return *this;
  }

  // Stream output
  friend std::ostream& operator<<(std::ostream& os, expression const& expr) {
    if(!expr.monomials_.empty()) {
      bool print_plus = false;
      for(auto const& m : expr.monomials_) {
        os << (print_plus ? " + " : "") << m.second;
        if(!m.first.empty()) os << "*";
        os << m.first;
        print_plus = true;
      }
    } else
      os << scalar_traits<ScalarType>::make_const(0);
    return os;
  }

private:

  //
  // Implementations of arithmetic operations
  //

  //
  // Addition
  //

  // Addition
  // sum_type<ScalarType, S> == ScalarType
  template<typename S>
  inline expression add_impl(expression_t<S> const& expr,
                             std::true_type) const {
    expression res(*this);
    auto & res_mons = res.get_monomials();
    for(auto const& p : expr.get_monomials()) {
      auto it = res_mons.find(p.first);
      if(it == res_mons.end())
        res_mons.emplace(p);
      else {
        it->second = it->second + p.second;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          res_mons.erase(it);
      }
    }
    return res;
  }

  // Addition
  // sum_type<ScalarType, S> != ScalarType
  template<typename S>
  inline expression_t<sum_type<ScalarType, S>>
  add_impl(expression_t<S> const& expr, std::false_type) const {
    expression_t<sum_type<ScalarType, S>> res;
    auto & res_mons = res.get_monomials();
    auto const& m1 = monomials_;
    auto const& m2 = expr.get_monomials();
    auto z1 = scalar_traits<ScalarType>::make_const(0);
    auto z2 = scalar_traits<S>::make_const(0);

    auto it1 = m1.begin();
    auto it2 = m2.begin();

    while(it1 != m1.end() && it2 != m2.end()) {
      if(it1->first < it2->first) {
        res_mons.emplace_hint(res_mons.end(), it1->first, it1->second + z2);
        ++it1;
      } else if(it2->first < it1->first) {
        res_mons.emplace_hint(res_mons.end(), it2->first, z1 + it2->second);
        ++it2;
      } else {
        auto val = it1->second + it2->second;
        if(!scalar_traits<sum_type<ScalarType, S>>::is_zero(val)) {
          res_mons.emplace_hint(res_mons.end(), it1->first, val);
        }
        ++it1;
        ++it2;
      }
    }

    if(it1 == m1.end() && it2 != m2.end()) {
      for(; it2 != m2.end(); ++it2)
        res_mons.emplace_hint(res_mons.end(), it2->first, z1 + it2->second);
    }
    if(it1 != m1.end() && it2 == m2.end()) {
      for(; it1 != m1.end(); ++it1)
        res_mons.emplace_hint(res_mons.end(), it1->first, it1->second + z2);
    }

    return res;
  }

  //
  // Subtraction
  //

  // Subtraction
  // diff_type<ScalarType, S> == ScalarType
  template<typename S>
  inline expression sub_impl(expression_t<S> const& expr,
                             std::true_type) const {
    expression res(*this);
    auto & res_mons = res.get_monomials();
    for(auto const& p : expr.get_monomials()) {
      auto it = res_mons.find(p.first);
      if(it == res_mons.end())
        res_mons.emplace(p.first, -p.second);
      else {
        it->second = it->second - p.second;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          res_mons.erase(it);
      }
    }
    return res;
  }

  // Subtraction
  // diff_type<ScalarType, S> != ScalarType
  template<typename S>
  inline expression_t<diff_type<ScalarType, S>>
  sub_impl(expression_t<S> const& expr, std::false_type) const {
    expression_t<diff_type<ScalarType, S>> res;
    auto & res_mons = res.get_monomials();
    auto const& m1 = monomials_;
    auto const& m2 = expr.get_monomials();
    auto z1 = scalar_traits<ScalarType>::make_const(0);
    auto z2 = scalar_traits<S>::make_const(0);

    auto it1 = m1.begin();
    auto it2 = m2.begin();

    while(it1 != m1.end() && it2 != m2.end()) {
      if(it1->first < it2->first) {
        res_mons.emplace_hint(res_mons.end(), it1->first, it1->second - z2);
        ++it1;
      } else if(it2->first < it1->first) {
        res_mons.emplace_hint(res_mons.end(), it2->first, z1 - it2->second);
        ++it2;
      } else {
        auto val = it1->second - it2->second;
        if(!scalar_traits<diff_type<ScalarType, S>>::is_zero(val)) {
          res_mons.emplace_hint(res_mons.end(), it1->first, val);
        }
        ++it1;
        ++it2;
      }
    }

    if(it1 == m1.end() && it2 != m2.end()) {
      for(; it2 != m2.end(); ++it2)
        res_mons.emplace_hint(res_mons.end(), it2->first, z1 - it2->second);
    }
    if(it1 != m1.end() && it2 == m2.end()) {
      for(; it1 != m1.end(); ++it1)
        res_mons.emplace_hint(res_mons.end(), it1->first, it1->second - z2);
    }

    return res;
  }

  //
  // Unary minus
  //

  // Unary minus: minus_type<ScalarType> == ScalarType
  inline expression unary_minus_impl(std::true_type) const {
    expression res(*this);
    for(auto & p : res.monomials_) p.second = -p.second;
    return res;
  }

  // Unary minus: minus_type<ScalarType> != ScalarType
  inline expression_t<minus_type<ScalarType>>
  unary_minus_impl(std::false_type) const {
    expression_t<minus_type<ScalarType>> res;
    auto & res_mons = res.get_monomials();
    for(auto const& p : monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, -p.second);
    return res;
  }

  //
  // Multiplication by scalar
  //

  // Multiplication by scalar (postfix form)
  // mul_type<ScalarType, S> == ScalarType
  template<typename S>
  inline expression mul_const_postfix_impl(S&& alpha, std::true_type) const {
    expression res(*this);
    for(auto & p : res.monomials_) p.second = p.second * alpha;
    return res;
  }

  // Multiplication by scalar (postfix form)
  // mul_type<ScalarType, S> != ScalarType
  template<typename S>
  inline expression_t<mul_type<ScalarType, S>>
  mul_const_postfix_impl(S&& alpha, std::false_type) const {
    expression_t<mul_type<ScalarType, S>> res;
    auto & res_mons = res.get_monomials();
    for(auto const& p : monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, p.second * alpha);
    return res;
  }

  // Multiplication by scalar (prefix form)
  // mul_type<S, ScalarType> == ScalarType
  template<typename S>
  inline expression mul_const_prefix_impl(S&& alpha, std::true_type) const {
    expression res(*this);
    for(auto & p : res.monomials_) p.second = alpha * p.second;
    return res;
  }

  // Multiplication by scalar (prefix form)
  // mul_type<S, ScalarType> != ScalarType
  template<typename S>
  inline expression_t<mul_type<S, ScalarType>>
  mul_const_prefix_impl(S&& alpha, std::false_type) const {
    expression_t<mul_type<S, ScalarType>> res;
    auto & res_mons = res.get_monomials();
    for(auto const& p : monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, alpha * p.second);
    return res;
  }

  //
  // Addition of scalar
  //

  // Addition of scalar (postfix form)
  // sum_type<ScalarType, S> == ScalarType
  template<typename S>
  inline expression add_const_postfix_impl(S&& alpha, std::true_type) const {
    expression res(*this);
    auto & res_mons = res.get_monomials();
    if(!scalar_traits<S>::is_zero(alpha)) {
      auto it = res_mons.find(monomial_t{});
      if(it == res_mons.end()) {
        res_mons.emplace_hint(res_mons.begin(), monomial_t{}, alpha);
      } else {
        it->second = it->second + alpha;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          res_mons.erase(it);
      }
    }
    return res;
  }

  // Addition of scalar (postfix form)
  // sum_type<ScalarType, S> != ScalarType
  template<typename S>
  inline expression_t<sum_type<ScalarType, S>>
  add_const_postfix_impl(S&& alpha, std::false_type) const {
    using res_s_t = sum_type<ScalarType, S>;
    expression_t<res_s_t> res;
    auto z = scalar_traits<S>::make_const(0);
    auto & res_mons = res.get_monomials();
    for(auto const& p : monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, p.second + z);
    if(!scalar_traits<S>::is_zero(alpha)) {
      auto it = res_mons.find(monomial_t{});
      if(it == res_mons.end()) {
        res_mons.emplace_hint(res_mons.begin(), monomial_t{}, alpha);
      } else {
        it->second = it->second + alpha;
        if(scalar_traits<res_s_t>::is_zero(it->second))
          res_mons.erase(it);
      }
    }
    return res;
  }

  // Addition of scalar (prefix form)
  // sum_type<S, ScalarType> == ScalarType
  template<typename S>
  inline expression add_const_prefix_impl(S&& alpha, std::true_type) const {
    expression res(*this);
    auto & res_mons = res.get_monomials();
    if(!scalar_traits<S>::is_zero(alpha)) {
      auto it = res_mons.find(monomial_t{});
      if(it == res_mons.end()) {
        res_mons.emplace_hint(res_mons.begin(), monomial_t{}, alpha);
      } else {
        it->second = alpha + it->second;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          res_mons.erase(it);
      }
    }
    return res;
  }

  // Addition of scalar (prefix form)
  // sum_type<S, ScalarType> != ScalarType
  template<typename S>
  inline expression_t<sum_type<S, ScalarType>>
  add_const_prefix_impl(S&& alpha, std::false_type) const {
    using res_s_t = sum_type<S, ScalarType>;
    expression_t<res_s_t> res;
    auto z = scalar_traits<S>::make_const(0);
    auto & res_mons = res.get_monomials();
    for(auto const& p : monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, z + p.second);
    if(!scalar_traits<S>::is_zero(alpha)) {
      auto it = res_mons.find(monomial_t{});
      if(it == res_mons.end()) {
        res_mons.emplace_hint(res_mons.begin(), monomial_t{}, alpha);
      } else {
        it->second = alpha + it->second;
        if(scalar_traits<res_s_t>::is_zero(it->second))
          res_mons.erase(it);
      }
    }
    return res;
  }

  //
  // Subtraction of scalar
  //

  // Subtraction of scalar (postfix form)
  // diff_type<ScalarType, S> == ScalarType
  template<typename S>
  inline expression sub_const_postfix_impl(S&& alpha, std::true_type) const {
    expression res(*this);
    auto & res_mons = res.get_monomials();
    if(!scalar_traits<S>::is_zero(alpha)) {
      auto it = res_mons.find(monomial_t{});
      if(it == res_mons.end()) {
        res_mons.emplace_hint(res_mons.begin(), monomial_t{}, -alpha);
      } else {
        it->second = it->second - alpha;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          res_mons.erase(it);
      }
    }
    return res;
  }

  // Subtraction of scalar (postfix form)
  // diff_type<ScalarType, S> != ScalarType
  template<typename S>
  inline expression_t<diff_type<ScalarType, S>>
  sub_const_postfix_impl(S&& alpha, std::false_type) const {
    using res_s_t = diff_type<ScalarType, S>;
    expression_t<res_s_t> res;
    auto z = scalar_traits<S>::make_const(0);
    auto & res_mons = res.get_monomials();
    for(auto const& p : monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, p.second - z);
    if(!scalar_traits<S>::is_zero(alpha)) {
      auto it = res_mons.find(monomial_t{});
      if(it == res_mons.end()) {
        auto z = scalar_traits<ScalarType>::make_const(0);
        res_mons.emplace_hint(res_mons.begin(), monomial_t{}, z - alpha);
      } else {
        it->second = it->second - alpha;
        if(scalar_traits<res_s_t>::is_zero(it->second))
          res_mons.erase(it);
      }
    }
    return res;
  }
};

// Aliases for specific scalar types
template<typename... IndexTypes>
using expr_real = expression<double, IndexTypes...>;
template<typename... IndexTypes>
using expr_complex = expression<std::complex<double>, IndexTypes...>;

} // namespace libcommute

#endif
