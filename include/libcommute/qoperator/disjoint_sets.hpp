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
#ifndef LIBCOMMUTE_QOPERATOR_DISJOINT_SETS_HPP_
#define LIBCOMMUTE_QOPERATOR_DISJOINT_SETS_HPP_

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

namespace libcommute {

//
// Disjoint sets data structure
// Adapted from https://github.com/JuliaCollections/DataStructures.jl
//

class disjoint_sets {
  mutable std::vector<std::size_t> parents_;
  std::vector<std::size_t> ranks_;
  std::size_t n_sets_;

public:

  disjoint_sets(std::size_t n_sets)
    : parents_(n_sets), ranks_(n_sets, 0), n_sets_(n_sets) {
    std::iota(parents_.begin(), parents_.end(), 0);
  }

  // Add a new singleton set
  inline std::size_t make_set() {
    parents_.push_back(size());
    ranks_.push_back(0);
    n_sets_ += 1;
    return parents_.back();
  }

  // Add n singleton sets
  inline void make_sets(std::size_t n) {
    parents_.reserve(size() + n);
    ranks_.reserve(size() + n);
    for(std::size_t i = 0; i < n; ++i) {
      parents_.push_back(size());
      ranks_.push_back(0);
    }
    n_sets_ += n;
  }

  // Size of 'parents' array
  inline size_t size() const { return parents_.size(); }

  // Number of disjoint sets
  inline size_t n_sets() const { return n_sets_; }

  // Find representative element of the set containing x (with path compression)
  std::size_t find_root(std::size_t x) const {
    assert(x < size());
    std::size_t p = parents_[x];
    if(parents_[p] != p) {
      parents_[x] = p = find_root(p);
    }
    return p;
  }

  // Do elements x and y belong to the same set?
  inline bool in_same_set(std::size_t x, std::size_t y) const {
    assert(x < size());
    assert(y < size());
    return find_root(x) == find_root(y);
  }

  // Merge sets represented by x and y
  std::size_t root_union(std::size_t x, std::size_t y) {
    assert(x < size());
    assert(y < size());
    assert(parents_[x] == x);
    assert(parents_[y] == y);

    std::size_t x_rank = ranks_[x];
    std::size_t y_rank = ranks_[y];
    if(x_rank < y_rank)
      std::swap(x, y);
    else if(x_rank == y_rank)
      ranks_[x] += 1;

    parents_[y] = x;
    n_sets_ -= 1;
    return x;
  }

  // Merge sets containing x and y
  std::size_t set_union(std::size_t x, std::size_t y) {
    std::size_t x_root = find_root(x);
    std::size_t y_root = find_root(y);
    return x_root != y_root ? root_union(x_root, y_root) : x_root;
  }

  // Compress all sets (make every element represented by its parent)
  void compress_sets() {
    for(std::size_t x : parents_) find_root(x);
  }

  // Normalize all sets (make every representative be the smallest in its set)
  void normalize_sets() {
    for(std::size_t x = 0; x < size(); ++x) {
      std::size_t p = parents_[x];
      if(x > p || parents_[p] != p)
        parents_[x] = parents_[p];
      else {
        parents_[p] = x;
        parents_[x] = x;
      }
    }
  }
};

} // namespace libcommute

#endif
