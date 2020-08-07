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
#ifndef LIBCOMMUTE_EXPRESSION_EXPRESSION_HPP_
#define LIBCOMMUTE_EXPRESSION_EXPRESSION_HPP_

#include "generator.hpp"
#include "monomial.hpp"
#include "../scalar_traits.hpp"
#include "../utility.hpp"
#include "../metafunctions.hpp"

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

template<typename ScalarType, typename... IndexTypes>
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
    if(!scalar_traits<S>::is_zero(x))
      normalize_and_store(monomial_t(monomial), x, monomials_);
  }

private:

  // Internal: Construct from a linear function of generators
  expression(typename monomial_t::generator_type::linear_function_t const& f) {
    if(!scalar_traits<ScalarType>::is_zero(f.const_term))
      monomials_.emplace(monomial_t{}, f.const_term);
    for(auto const& g : f.terms) {
      auto v = scalar_traits<ScalarType>::make_const(g.second);
      if(!scalar_traits<ScalarType>::is_zero(v))
        *this += expression(v, monomial_t({g.first->clone()}));
    }
  }

public:

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

  // Set expression to zero
  inline void clear() { monomials_.clear(); }

  // Equality
  inline friend bool operator==(expression const& e1, expression const& e2) {
    return e1.monomials_ == e2.monomials_;
  }
  inline friend bool operator!=(expression const& e1, expression const& e2) {
    return !operator==(e1, e2);
  }

  //
  // Iteration interface
  //

  // Constant iterator over monomials
  class const_iterator;

  inline const_iterator begin() const noexcept {
    return const_iterator(monomials_.begin());
  }
  inline const_iterator cbegin() const noexcept {
    return const_iterator(monomials_.cbegin());
  }

  inline const_iterator end() const noexcept {
    return const_iterator(monomials_.end());
  }
  inline const_iterator cend() const noexcept {
    return const_iterator(monomials_.cend());
  }

  // Apply functor 'f' to all monomial/coefficient pairs
  // and replace the coefficients with the functor return values.
  template<typename F,
#ifndef LIBCOMMUTE_NO_STD_INVOKE_RESULT
           typename NewScalarType = std::invoke_result_t<F,
             expression::monomial_t const&,
             expression::scalar_type const&>>
#else
           typename NewScalarType = invoke_result_t<F,
             expression::monomial_t const&,
             expression::scalar_type const&>>
#endif
  friend expression_t<NewScalarType> transform(expression const& expr, F&& f) {
    expression_t<NewScalarType> res;
    auto & res_mons = res.get_monomials();
    for(auto const& m : expr.monomials_) {
      auto val = f(m.first, m.second);
      if(!scalar_traits<NewScalarType>::is_zero(val))
        res_mons.emplace_hint(res_mons.end(), m.first, val);
    }
    return res;
  }

  // Hermitian conjugate
  friend expression conj(expression const& expr) {
    expression res;
    expression m_contrib;
    typename monomial_t::generator_type::linear_function_t f;
    for(auto const& m : expr.monomials_) {
      auto coeff = scalar_traits<ScalarType>::conj(m.second);
      if(!scalar_traits<ScalarType>::is_zero(coeff)) {
        m_contrib = expression(coeff);
        for(auto m_it = m.first.rbegin(); m_it != m.first.rend(); ++m_it) {
          (*m_it).conj(f);
          m_contrib *= expression(f);
        }
      }
      res += m_contrib;
    }
    return res;
  }

  //
  // Arithmetics
  //

  // Addition
  template<typename S>
  expression_t<sum_type<ScalarType, S>>
  operator+(expression_t<S> const& expr) const {
    return add_impl(expr, std::is_same<sum_type<ScalarType, S>, ScalarType>());
  }

  // Subtraction
  template<typename S>
  expression_t<diff_type<ScalarType, S>>
  operator-(expression_t<S> const& expr) const {
    return sub_impl(expr, std::is_same<diff_type<ScalarType, S>, ScalarType>());
  }

  // Multiplication
  template<typename S>
  expression_t<mul_type<ScalarType, S>>
  operator*(expression_t<S> const& expr) const {
    expression_t<mul_type<ScalarType, S>> res(*this);
    res *= expr;
    return res;
  }

  //
  // Compound assignments
  //

  // Compound assignment/addition
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

  // Compound assignment/subtraction
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

  // Compound assignment/multiplication
  template<typename S>
  expression & operator*=(expression_t<S> const& expr) {
    monomials_map_t res_map;
    for(auto const& m1 : monomials_) {
      for(auto const& m2 : expr.get_monomials()) {
        normalize_and_store(concatenate(m1.first, m2.first),
                            m1.second * m2.second,
                            res_map);
      }
    }
    std::swap(monomials_, res_map);
    return *this;
  }

  // Unary minus
  template<typename S = ScalarType>
  auto operator-() const -> expression_t<minus_type<S>> {
    return unary_minus_impl<S>(std::is_same<minus_type<S>, ScalarType>());
  }

  //
  // Arithmetics involving constants
  //

  // Multiplication by scalar (postfix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  auto operator*(S const& alpha) const
    -> expression_t<mul_type<ScalarType, S>> {
    if(scalar_traits<remove_cvref_t<S>>::is_zero(alpha))
      return {};
    else
      return mul_const_postfix_impl(
        alpha,
        std::is_same<mul_type<ScalarType, S>, ScalarType>()
      );
  }

  // Multiplication by scalar (prefix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  friend auto operator*(S const& alpha, expression const& expr)
    -> expression_t<mul_type<S, ScalarType>> {
    if(scalar_traits<remove_cvref_t<S>>::is_zero(alpha))
      return {};
    else
      return expr.mul_const_prefix_impl(
        alpha,
        std::is_same<mul_type<S, ScalarType>, ScalarType>()
      );
  }

  // Addition of scalar (postfix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  auto operator+(S const& alpha) const ->
    expression_t<sum_type<ScalarType, S>> {
    return add_const_postfix_impl(
      alpha,
      std::is_same<sum_type<ScalarType, S>, ScalarType>()
    );
  }

  // Addition of scalar (prefix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  friend auto operator+(S const& alpha, expression const& expr)
    -> expression_t<sum_type<S, ScalarType>> {
    return expr.add_const_prefix_impl(
      alpha,
      std::is_same<sum_type<S, ScalarType>, ScalarType>()
    );
  }

  // Subtraction of scalar (postfix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  auto operator-(S const& alpha) const ->
    expression_t<diff_type<ScalarType, S>> {
    return sub_const_postfix_impl(
      alpha,
      std::is_same<diff_type<ScalarType, S>, ScalarType>()
    );
  }

  // Subtraction of scalar (prefix form)
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  friend auto operator-(S const& alpha, expression const& expr)
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
  expression & operator*=(S const& alpha) {
    if(scalar_traits<remove_cvref_t<S>>::is_zero(alpha))
      monomials_.clear();
    else
      for(auto & p : monomials_) p.second = p.second * alpha;
    return *this;
  }

  // Compound assignments/addition of scalar
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  expression & operator+=(S const& alpha) {
    using s_t = remove_cvref_t<S>;
    if(!scalar_traits<s_t>::is_zero(alpha)) {
      auto it = monomials_.find(monomial_t{});
      if(it == monomials_.end()) {
        auto val = scalar_traits<ScalarType>::make_const(0);
        val = val + alpha;
        monomials_.emplace_hint(monomials_.begin(), monomial_t{}, val);
      } else {
        it->second = it->second + alpha;
        if(scalar_traits<ScalarType>::is_zero(it->second))
          monomials_.erase(it);
      }
    }
    return *this;
  }

  // Compound assignments/subtraction of scalar
  template<typename S, typename = disable_for_expression<remove_cvref_t<S>>>
  expression & operator-=(S const& alpha) {
    using s_t = remove_cvref_t<S>;
    if(!scalar_traits<s_t>::is_zero(alpha)) {
      auto it = monomials_.find(monomial_t{});
      if(it == monomials_.end()) {
        monomials_.emplace_hint(monomials_.begin(), monomial_t{}, -alpha);
      } else {
        it->second = it->second - alpha;
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
  template<typename S>
  inline expression unary_minus_impl(std::true_type) const {
    expression res(*this);
    for(auto & p : res.monomials_) p.second = -p.second;
    return res;
  }

  // Unary minus: minus_type<ScalarType> != ScalarType
  template<typename S>
  inline expression_t<minus_type<S>> unary_minus_impl(std::false_type) const {
    expression_t<minus_type<S>> res;
    auto & res_mons = res.get_monomials();
    for(auto const& p : monomials_)
      res_mons.emplace_hint(res_mons.end(), p.first, -p.second);
    return res;
  }

   //
  // Multiplication
  //

  // Store monomial in a map while taking care of possible collisions
  static void store_monomial(monomial_t&& m,
                             scalar_type coeff,
                             monomials_map_t& target) {
    bool is_new_monomial;
    typename monomials_map_t::iterator it;
    std::tie(it, is_new_monomial) = target.emplace(m, coeff);
    if(!is_new_monomial) {
      it->second = it->second + coeff;
      if(scalar_traits<ScalarType>::is_zero(it->second))
        target.erase(it);
    }
  }

  // Normalize a monomial and store in a map
  static void normalize_and_store(monomial_t&& m,
                                  scalar_type coeff,
                                  monomials_map_t& target) {
    // Normalization is done by means of a simple bubble sort algorithm.
    // Apart from sorting elements, this function keeps track of the coefficient
    // and recursively calls itself if a permutation of two operators produces a
    // new monomial.

    using generator_t = generator<IndexTypes...>;

    if(m.size() >= 2) {
      bool is_swapped;
      typename generator_t::linear_function_t f;

      // Process linear terms
      auto process_f = [&](int n) {
        // Process monomial generated by the constant term in 'f'
        if(f.const_term != 0) {
          normalize_and_store(concatenate(
                                std::make_pair(m.begin(), m.begin()+n-1),
                                std::make_pair(m.begin()+n+1, m.end())
                              ),
                              coeff * f.const_term,
                              target);
        }

        // Process monomials generated by the rest of the terms in 'f'
        for(auto const& t : f.terms) {
          normalize_and_store(concatenate(
                  std::make_pair(m.begin(), m.begin()+n-1),
                  *t.first,
                  std::make_pair(m.begin()+n+1, m.end())
                ),
                coeff * t.second,
                target);
        }
      };

      do {
        is_swapped = false;
        for(size_t n = 1; n < m.size(); ++n) {
          // Pick a pair of generators in m, m[n-1] and m[n]
          generator_t& prev_gen = m[n - 1];
          generator_t& cur_gen = m[n];

          if(prev_gen > cur_gen) { // Reordering is needed
            double c = swap_with(prev_gen, cur_gen, f);

            process_f(n);

            if(c == 0) { // We have to stop sorting here as all contributions
                         // to the currently processed monomial are already
                         // taken care of by process_f(n).
              return;
            } else { // Swap generators
              coeff = coeff * scalar_traits<ScalarType>::make_const(c);
              m.swap_generators(n - 1, n);
              is_swapped = true;
            }
          } else { // Generators in m[n-1]*m[n] are already in order.
                   // Is a simplification of the product possible?
            if(simplify_prod(prev_gen, cur_gen, f)) {
              process_f(n);
              return;
            }
          }
        }
      } while(is_swapped);

      // Check that coefficient in front of this monomial is not zero
      if(scalar_traits<ScalarType>::is_zero(coeff)) return;
    }

    // Store the result
    reduce_powers_and_store(std::move(m), coeff, target);
  }

  static void reduce_powers_and_store(monomial_t&& m,
                                      scalar_type coeff,
                                      monomials_map_t& target) {
    if(m.empty()) {
      store_monomial(std::move(m), coeff, target);
      return;
    }

    linear_function<std::unique_ptr<generator<IndexTypes...>>> f;
    auto it = m.begin(), end_it = m.end();
    auto next_it = it + 1;
    int power = 1;
    for(;it != end_it; ++it, ++next_it) {
      if(next_it == end_it) {
        store_monomial(std::move(m), coeff, target);
        return;
      } else {
        if(*next_it == *it) {
          ++power;
          if(power > 2 && it->reduce_power(power, f)) {
            if(f.vanishing()) return;

            auto v = scalar_traits<ScalarType>::make_const(f.const_term);
            if(!scalar_traits<ScalarType>::is_zero(v))
              reduce_powers_and_store(concatenate(
                                      std::make_pair(m.begin(), it-power+2),
                                      std::make_pair(next_it+1, end_it)
                                      ),
                                      coeff * v,
                                      target);
            for(auto const& g : f.terms) {
              v = scalar_traits<ScalarType>::make_const(g.second);
              if(!scalar_traits<ScalarType>::is_zero(v))
                normalize_and_store(concatenate(
                                      std::make_pair(m.begin(), it-power+2),
                                      *g.first,
                                      std::make_pair(next_it+1, end_it)
                                    ),
                                    coeff * v,
                                    target);
            }
            return;
          }
        } else {
          power = 1;
        }
      }
    }
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
  inline
  expression add_const_postfix_impl(S const& alpha, std::true_type) const {
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
  add_const_postfix_impl(S const& alpha, std::false_type) const {
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
  inline
  expression add_const_prefix_impl(S const& alpha, std::true_type) const {
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
  add_const_prefix_impl(S const& alpha, std::false_type) const {
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
  inline
  expression sub_const_postfix_impl(S const& alpha, std::true_type) const {
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
  sub_const_postfix_impl(S const& alpha, std::false_type) const {
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

// Constant iterator over monomials
template<typename ScalarType, typename... IndexTypes>
class expression<ScalarType, IndexTypes...>::const_iterator {

  using map_it = typename monomials_map_t::const_iterator;
  map_it m_it_;

public:

  struct value_type {
    monomial_t const& monomial;
    ScalarType const& coeff;
    value_type(monomial_t const& m, ScalarType const& c)
      : monomial(m), coeff(c) {}
  };

  using iterator_category = std::bidirectional_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using pointer = std::unique_ptr<value_type>;
  using reference = value_type;

  explicit const_iterator(map_it const& m_it) : m_it_(m_it) {}

  const_iterator() = default;
  const_iterator(const_iterator const&) = default;
  const_iterator(const_iterator&&) noexcept = default;
  const_iterator& operator=(const_iterator const&) = default;
  const_iterator& operator=(const_iterator&&) noexcept = default;

  // Increments
  const_iterator& operator++() { ++m_it_; return *this;}
  const_iterator operator++(int) {
    const_iterator retval = *this;
    ++(*this);
    return retval;
  }

  // Decrements
  const_iterator& operator--() { --m_it_; return *this;}
  const_iterator operator--(int) {
    const_iterator retval = *this;
    --(*this);
    return retval;
  }

  // Equality
  bool operator==(const_iterator const& it) const {return m_it_ == it.m_it_;}
  bool operator!=(const_iterator const& it) const {return !(*this == it);}

  // Dereference
  reference operator*() const { return {m_it_->first, m_it_->second}; }
  pointer operator->() const {
#ifndef LIBCOMMUTE_NO_STD_MAKE_UNIQUE
    using std::make_unique;
#endif
    return make_unique<value_type>(m_it_->first, m_it_->second);
  }

  // swap()
  friend void swap(const_iterator& lhs, const_iterator& rhs) {
    std::swap(lhs.m_it_, rhs.m_it_);
  }
};

namespace static_indices {

// Aliases for specific scalar types
template<typename... IndexTypes>
using expr_real = expression<double, IndexTypes...>;
template<typename... IndexTypes>
using expr_complex = expression<std::complex<double>, IndexTypes...>;

} // namespace libcommute::static_indices
} // namespace libcommute

#if __cplusplus >= 201703L
#include "dyn_indices.hpp"

namespace libcommute {
namespace dynamic_indices {

// Aliases for specific scalar types
using expr_real = expression<double, dyn_indices>;
using expr_complex = expression<std::complex<double>, dyn_indices>;

} // namespace libcommute::dynamic_indices
} // namespace libcommute

#endif

#endif
