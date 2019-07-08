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
#ifndef LIBCOMMUTE_GENERATOR_HPP_
#define LIBCOMMUTE_GENERATOR_HPP_

#include "../utility.hpp"

#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

namespace libcommute {

//
// List of predefined algebra IDs
//

static constexpr int FERMION_ALGEBRA_ID = -3;
static constexpr int BOSON_ALGEBRA_ID = -2;
static constexpr int SPIN_ALGEBRA_ID = -1;

//
// Abstract algebra generator
//

template<typename... IndexTypes>
class generator {

public:

  using index_types = std::tuple<IndexTypes...>;

  generator(IndexTypes... indices) : indices_(indices...) {}
  generator() = delete;
  generator(generator const&) = default;
  generator(generator&&) noexcept = default;
  generator& operator=(generator const&) = default;
  generator& operator=(generator&&) noexcept = default;
  virtual ~generator() {}

  // Get ID of the algebra this generator belongs to
  virtual int algebra_id() const = 0;

  // Make a smart pointer that manages a copy of this generator
  virtual std::unique_ptr<generator> clone() const = 0;

  // Comparisons
  friend bool operator==(generator const& g1, generator const& g2) {
    return g1.algebra_id() == g2.algebra_id() && g1.equal(g2);
  }
  friend bool operator!=(generator const& g1, generator const& g2) {
    return !operator==(g1, g2);
  }
  friend bool operator<(generator const& g1, generator const& g2) {
    if(g1.algebra_id() != g2.algebra_id())
      return g1.algebra_id() < g2.algebra_id();
    else
      return g1.less(g2);
  }
  friend bool operator>(generator const& g1, generator const& g2) {
    if(g1.algebra_id() != g2.algebra_id())
      return g1.algebra_id() > g2.algebra_id();
    else
      return g1.greater(g2);
  }

  // Stream output
  friend std::ostream & operator<<(std::ostream & os, generator const& g) {
    return g.print(os);
  }

protected:

  index_types indices_;

  // Check two generators of the same algebra for equality
  virtual bool equal(generator const& g) const = 0;
  // Ordering
  virtual bool less(generator const& g) const = 0;
  virtual bool greater(generator const& g) const = 0;
  // Print to stream
  virtual std::ostream & print(std::ostream & os) const = 0;
};

//
// Generator of the fermionic algebra
//

template<typename... IndexTypes>
class fermion_generator : public generator<IndexTypes...> {

  using base = generator<IndexTypes...>;

public:

  // Get ID of the algebra this generator belongs to
  virtual int algebra_id() const override { return FERMION_ALGEBRA_ID; }

  // Value semantics
  fermion_generator(bool dagger, IndexTypes&&... indices) :
    base(std::forward<IndexTypes>(indices)...), dagger_(dagger) {}
  fermion_generator() = delete;
  fermion_generator(fermion_generator const&) = default;
  fermion_generator(fermion_generator&&) noexcept = default;
  fermion_generator& operator=(fermion_generator const&) = default;
  fermion_generator& operator=(fermion_generator&&) noexcept = default;
  virtual ~fermion_generator() {}

  // Make a smart pointer that manages a copy of this generator
  virtual std::unique_ptr<base> clone() const override {
    return make_unique<fermion_generator>(*this);
  }

protected:
  // Creation or annihilation operator?
  bool dagger_;

  // Check two generators of the same algebra for equality
  virtual bool equal(base const& g) const override {
    auto const& f_g =  dynamic_cast<fermion_generator const&>(g);
    return dagger_ == f_g.dagger_ && this->indices_ == f_g.indices_;
  }

  // Ordering
  virtual bool less(base const& g) const override {
    auto const& f_g =  dynamic_cast<fermion_generator const&>(g);
    // Example: c+_1 < c+_2 < c+_3 < c_3 < c_2 < c_1
    if(this->dagger_ != f_g.dagger_)
      return (this->dagger_ > f_g.dagger_);
    else
      return this->dagger_ ? this->indices_ < f_g.indices_ :
                             this->indices_ > f_g.indices_;
  }
  virtual bool greater(base const& g) const override {
    auto const& f_g =  dynamic_cast<fermion_generator const&>(g);
    // Example: c_1 > c_2 > c_3 > c+_3 > c+_2 > c+_1
    if(this->dagger_ != f_g.dagger_)
      return (this->dagger_ < f_g.dagger_);
    else
      return this->dagger_ ? this->indices_ > f_g.indices_ :
                             this->indices_ < f_g.indices_;
  }

  // Print to stream
  virtual std::ostream & print(std::ostream & os) const override {
    os << "C" << (this->dagger_ ? "+" : "") << "(";
    print_tuple(os, this->indices_);
    return os << ")";
  }
};

// Convenience factory function
template<typename... IndexTypes>
inline fermion_generator<typename c_str_to_string_t<IndexTypes>::type...>
make_fermion(bool dagger, IndexTypes&&... indices) {
  return {dagger, std::forward<IndexTypes>(indices)...};
}

// Check if generator belongs to the fermionic algebra
template <typename... IndexTypes>
bool is_fermion(generator<IndexTypes...> const& gen) {
  return gen.algebra_id() == FERMION_ALGEBRA_ID;
}

//
// Generator of the bosonic algebra
//

template<typename... IndexTypes>
class boson_generator : public generator<IndexTypes...> {

  using base = generator<IndexTypes...>;

public:

  // Get ID of the algebra this generator belongs to
  virtual int algebra_id() const override { return BOSON_ALGEBRA_ID; }

  // Value semantics
  boson_generator(bool dagger, IndexTypes&&... indices) :
    base(std::forward<IndexTypes>(indices)...), dagger_(dagger) {}
  boson_generator() = delete;
  boson_generator(boson_generator const&) = default;
  boson_generator(boson_generator&&) noexcept = default;
  boson_generator& operator=(boson_generator const&) = default;
  boson_generator& operator=(boson_generator&&) noexcept = default;
  virtual ~boson_generator() {}

  // Make a smart pointer that manages a copy of this generator
  virtual std::unique_ptr<base> clone() const override {
    return make_unique<boson_generator>(*this);
  }

protected:
  // Creation or annihilation operator?
  bool dagger_;

  // Check two generators of the same algebra for equality
  virtual bool equal(base const& g) const override {
    auto const& b_g =  dynamic_cast<boson_generator const&>(g);
    return dagger_ == b_g.dagger_ && this->indices_ == b_g.indices_;
  }

  // Ordering
  virtual bool less(base const& g) const override {
    auto const& b_g =  dynamic_cast<boson_generator const&>(g);
    // Example: a+_1 < a+_2 < a+_3 < a_3 < a_2 < a_1
    if(this->dagger_ != b_g.dagger_)
      return (this->dagger_ > b_g.dagger_);
    else
      return this->dagger_ ? this->indices_ < b_g.indices_ :
                             this->indices_ > b_g.indices_;
  }
  virtual bool greater(base const& g) const override {
    auto const& b_g =  dynamic_cast<boson_generator const&>(g);
    // Example: a_1 > a_2 > a_3 > a+_3 > a+_2 > a+_1
    if(this->dagger_ != b_g.dagger_)
      return (this->dagger_ < b_g.dagger_);
    else
      return this->dagger_ ? this->indices_ > b_g.indices_ :
                             this->indices_ < b_g.indices_;
  }

  // Print to stream
  virtual std::ostream & print(std::ostream & os) const override {
    os << "A" << (this->dagger_ ? "+" : "") << "(";
    print_tuple(os, this->indices_);
    return os << ")";
  }
};

// Convenience factory function
template<typename... IndexTypes>
inline boson_generator<typename c_str_to_string_t<IndexTypes>::type...>
make_boson(bool dagger, IndexTypes&&... indices) {
  return {dagger, std::forward<IndexTypes>(indices)...};
}

// Check if generator belongs to the bosonic algebra
template <typename... IndexTypes>
bool is_boson(generator<IndexTypes...> const& gen) {
  return gen.algebra_id() == BOSON_ALGEBRA_ID;
}

//
// Generator of the spin algebra
//

// S_+, S_- or S_z
enum spin_component : int {plus = 0, minus = 1, z = 2};

template<typename... IndexTypes>
class spin_generator : public generator<IndexTypes...> {

  using base = generator<IndexTypes...>;

public:

  // Get ID of the algebra this generator belongs to
  virtual int algebra_id() const override { return SPIN_ALGEBRA_ID; }

  // Value semantics
  spin_generator(spin_component c, IndexTypes&&... indices) :
    base(std::forward<IndexTypes>(indices)...), multiplicity_(2), c_(c) {}
  spin_generator(double spin, spin_component c, IndexTypes&&... indices) :
    base(std::forward<IndexTypes>(indices)...), multiplicity_(2*spin+1), c_(c) {
      // Multiplicity has to be integer
      assert(2*spin == int(spin*2));
    }

  spin_generator() = delete;
  spin_generator(spin_generator const&) = default;
  spin_generator(spin_generator&&) noexcept = default;
  spin_generator& operator=(spin_generator const&) = default;
  spin_generator& operator=(spin_generator&&) noexcept = default;
  virtual ~spin_generator() {}

  // Make a smart pointer that manages a copy of this generator
  virtual std::unique_ptr<base> clone() const override {
    return make_unique<spin_generator>(*this);
  }

protected:
  // Multiplicity, 2S+1
  int multiplicity_;

  // Creation or annihilation operator?
  spin_component c_;

  // Check two generators of the same algebra for equality
  virtual bool equal(base const& g) const override {
    auto const& b_g =  dynamic_cast<spin_generator const&>(g);
    return multiplicity_ == b_g.multiplicity_ &&
                      c_ == b_g.c_ && this->indices_ == b_g.indices_;
  }

  // Ordering
  virtual bool less(base const& g) const override {
    auto const& s_g =  dynamic_cast<spin_generator const&>(g);
    // Example: S1/2+_1 < S1/2-_1 < S1/2z_1 < S1/2+_2 < S1/2-_2 < S1/2z_2 <
    //          S3/2+_1 < S3/2-_1 < S3/2z_1 < S3/2+_2 < S3/2-_2 < S3/2z_2
    if(this->multiplicity_ != s_g.multiplicity_)
      return (this->multiplicity_ < s_g.multiplicity_);
    else if(this->indices_ != s_g.indices_)
      return (this->indices_ < s_g.indices_);
    else
      return this->c_ < s_g.c_;
  }
  virtual bool greater(base const& g) const override {
    auto const& s_g =  dynamic_cast<spin_generator const&>(g);
    // Example: S3/2z_2 > S3/2-_2 > S3/2+_2 > S3/2z_1 > S3/2-_1 > S3/2+_1 >
    //          S1/2z_2 > S1/2-_2 > S1/2+_2 > S1/2z_1 > S1/2-_1 > S1/2+_1
    if(this->multiplicity_ != s_g.multiplicity_)
      return (this->multiplicity_ > s_g.multiplicity_);
    if(this->indices_ != s_g.indices_)
      return (this->indices_ > s_g.indices_);
    else
      return this->c_ > s_g.c_;
  }

  // Print to stream
  virtual std::ostream & print(std::ostream & os) const override {
    os << "S";
    if(multiplicity_ != 2) {
      if(multiplicity_ % 2 == 0)
        os << (multiplicity_ - 1) << "/2";
      else
        os << ((multiplicity_ - 1) / 2);
    }
    switch(this->c_) {
      case spin_component::plus: os << "+"; break;
      case spin_component::minus: os << "-"; break;
      case spin_component::z: os << "z"; break;
    }
    os << "(";
    print_tuple(os, this->indices_);
    return os << ")";
  }
};

// Convenience factory function
template<typename... IndexTypes>
inline spin_generator<typename c_str_to_string_t<IndexTypes>::type...>
make_spin(spin_component c, IndexTypes&&... indices) {
  return {c, std::forward<IndexTypes>(indices)...};
}

template<typename... IndexTypes>
inline spin_generator<typename c_str_to_string_t<IndexTypes>::type...>
make_spin(double spin, spin_component c, IndexTypes&&... indices) {
  return {spin, c, std::forward<IndexTypes>(indices)...};
}

// Check if generator belongs to the spin algebra
template <typename... IndexTypes>
bool is_spin(generator<IndexTypes...> const& gen) {
  return gen.algebra_id() == SPIN_ALGEBRA_ID;
}

} // namespace libcommute

#endif
