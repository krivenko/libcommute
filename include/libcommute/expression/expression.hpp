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

  // Unary minus
  auto operator-() const -> expression_t<minus_type<ScalarType>> {
    return unary_minus_impl(std::is_same<minus_type<ScalarType>, ScalarType>());
  }

  // Multiplication by scalar (postfix form)
  template<typename S>
  auto operator*(S && alpha) const
    -> expression_t<mul_type<ScalarType, S>> {
    if(scalar_traits<typename remove_cvref<S>::type>::is_zero(alpha))
      return {};
    else
      return mul_const_postfix_impl(
        std::forward<S>(alpha),
        std::is_same<mul_type<ScalarType, S>, ScalarType>()
      );
  }

  // Multiplication by scalar (prefix form)
  template<typename S>
  friend auto operator*(S && alpha, expression const& expr)
    -> expression_t<mul_type<S, ScalarType>> {
    if(scalar_traits<typename remove_cvref<S>::type>::is_zero(alpha))
      return {};
    else
      return expr.mul_const_prefix_impl(
        std::forward<S>(alpha),
        std::is_same<mul_type<S, ScalarType>, ScalarType>()
      );
  }

  // Addition of scalar (postfix form)
  template<typename S>
  auto operator+(S && alpha) const -> expression_t<sum_type<ScalarType, S>> {
    return add_const_postfix_impl(
      std::forward<S>(alpha),
      std::is_same<sum_type<ScalarType, S>, ScalarType>()
    );
  }

  // Addition of scalar (prefix form)
  template<typename S>
  friend auto operator+(S && alpha, expression const& expr)
    -> expression_t<sum_type<S, ScalarType>> {
    return expr.add_const_prefix_impl(
      std::forward<S>(alpha),
      std::is_same<sum_type<S, ScalarType>, ScalarType>()
    );
  }

  // Subtraction of scalar (postfix form)
  template<typename S>
  auto operator-(S && alpha) const -> expression_t<diff_type<ScalarType, S>> {
    return sub_const_postfix_impl(
      std::forward<S>(alpha),
      std::is_same<diff_type<ScalarType, S>, ScalarType>()
    );
  }

  // Subtraction of scalar (prefix form)
  template<typename S>
  friend auto operator-(S && alpha, expression const& expr)
    -> expression_t<diff_type<S, ScalarType>> {
    using res_s_t = diff_type<S, ScalarType>;
    expression_t<res_s_t> res;
    auto z = scalar_traits<S>::zero();
    auto & res_mons = res.get_monomials();
    for(auto const& p : expr.monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, z - p.second);
    if(!scalar_traits<S>::is_zero(alpha)) {
      auto it = res_mons.find(monomial_t{});
      if(it == res_mons.end()) {
        auto z = scalar_traits<ScalarType>::zero();
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
  template<typename S>
  expression & operator*=(S && alpha) {
    if(scalar_traits<typename remove_cvref<S>::type>::is_zero(alpha))
      monomials_.clear();
    else
      for(auto & p : monomials_) p.second *= alpha;
    return *this;
  }

  // Compound assignments/addition of scalar
  template<typename S>
  expression & operator+=(S && alpha) {
    using s_t = typename remove_cvref<S>::type;
    if(!scalar_traits<s_t>::is_zero(alpha)) {
      auto it = monomials_.find(monomial_t{});
      if(it == monomials_.end()) {
        auto val = scalar_traits<ScalarType>::zero();
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
  template<typename S>
  expression & operator-=(S && alpha) {
    using s_t = typename remove_cvref<S>::type;
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
      os << scalar_traits<ScalarType>::zero();
    return os;
  }

  // TODO

private:

  //
  // Implementations of arithmetic operations
  //

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
  // Multiplication
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
  // Addition
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
    auto z = scalar_traits<S>::zero();
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
    auto z = scalar_traits<S>::zero();
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
  // Subtraction
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
    auto z = scalar_traits<S>::zero();
    auto & res_mons = res.get_monomials();
    for(auto const& p : monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, p.second - z);
    if(!scalar_traits<S>::is_zero(alpha)) {
      auto it = res_mons.find(monomial_t{});
      if(it == res_mons.end()) {
        auto z = scalar_traits<ScalarType>::zero();
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
using expression_real = expression<double, IndexTypes...>;
template<typename... IndexTypes>
using expression_complex = expression<std::complex<double>, IndexTypes...>;

} // namespace libcommute

#endif
