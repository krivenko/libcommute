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

#include <catch.hpp>

#include <libcommute/expression/factories.hpp>
#include <libcommute/loperator/compressed_state_view.hpp>
#include <libcommute/loperator/elementary_space_boson.hpp>
#include <libcommute/loperator/elementary_space_fermion.hpp>
#include <libcommute/loperator/hilbert_space.hpp>
#include <libcommute/loperator/loperator.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <vector>

using namespace libcommute;

TEST_CASE("View of a compressed state vector", "[compressed_state_view]") {

  using namespace static_indices;

  using hs_type = hilbert_space<int>;
  using state_vector = std::vector<double>;
  using view_type = compressed_state_view<state_vector>;

  SECTION("get_element(), update_add_element() and set_zeros()") {
    hs_type hs{make_space_fermion(0),
               make_space_fermion(1),
               make_space_boson(3, 2),
               make_space_boson(4, 3),
               make_space_boson(4, 4),
               make_space_boson(3, 5),
               make_space_boson(3, 6)};

    state_vector st(hs.dim());
    std::iota(st.begin(), st.end(), 0);

    // 1734 == (01) (10 11 00) (01 1 0)
    // strides = {1, 12, 576}
    sv_index_type const index = 1734;
    sv_index_type const serial_index = 6 + 44 * 12 + 1 * 576;

    SECTION("const") {
      auto view = compressed_state_view<state_vector const>(st, hs);

      CHECK(std::is_same<element_type_t<decltype(view)>, double>::value);
      CHECK(get_element(view, index) == serial_index);
    }

    SECTION("non-const") {
      auto view = compressed_state_view<state_vector>(st, hs);

      CHECK(std::is_same<element_type_t<decltype(view)>, double>::value);
      CHECK(get_element(view, index) == serial_index);

      state_vector st_ref(st.size());
      std::iota(st_ref.begin(), st_ref.end(), 0);
      st_ref[serial_index] = serial_index - 3;
      update_add_element(view, index, -3);
      CHECK(st == st_ref);

      set_zeros(view);
      CHECK(st == state_vector(st.size(), 0));
    }
  }

  SECTION("map_index()") {

    SECTION("Dense Hilbert space") {
      hs_type hs{make_space_fermion(0), make_space_fermion(1)};

      state_vector st(4);
      auto view = view_type(st, hs);
      for(sv_index_type index = 0; index < hs.vec_size(); ++index) {
        CHECK(view.map_index(index) == index);
      }
    }

    SECTION("Sparse Hilbert space") {
      hs_type hs{make_space_fermion(0),
                 make_space_fermion(1),
                 make_space_boson(3, 2),
                 make_space_boson(4, 3),
                 make_space_boson(4, 4),
                 make_space_boson(3, 5),
                 make_space_boson(3, 6)};

      state_vector st(hs.dim());
      auto view = view_type(st, hs);
      sv_index_type serial_index = 0;
      foreach(hs, [&view, &serial_index](sv_index_type index) {
        CHECK(view.map_index(index) == serial_index);
        ++serial_index;
      });
    }
  }

  SECTION("foreach()") {
    state_vector st;

    auto check_foreach = [](view_type const& view, state_vector& st) {
      // Start from 1 so that none of st's elements is vanishing
      std::iota(st.begin(), st.end(), 1);
      sv_index_type count = 0;
      foreach(view, [&](sv_index_type index, double a) {
        CHECK(sv_index_type(a) == view.map_index(index) + 1);
        ++count;
      });
      // Check that all elements of 'st' have been processed
      CHECK(count == st.size());
    };

    // Empty Hilbert space
    hs_type hs;
    foreach(view_type(st, hs), [](sv_index_type, double) { CHECK(false); });

    // Dense Hilbert space
    hs.add(make_space_fermion(0));
    hs.add(make_space_fermion(1));

    st.resize(4);
    REQUIRE(st.size() == hs.dim());
    check_foreach(view_type(st, hs), st);

    // Sparse Hilbert space
    hs.add(make_space_boson(3, 2));
    hs.add(make_space_boson(4, 3));
    hs.add(make_space_boson(4, 4));
    hs.add(make_space_boson(3, 5));
    hs.add(make_space_boson(3, 6));

    st.resize(1728);
    REQUIRE(st.size() == hs.dim());
    check_foreach(view_type(st, hs), st);
  }

  SECTION("make_comp_state_view() and make_const_comp_state_view()") {
    hs_type hs{make_space_fermion(0), make_space_fermion(1)};

    state_vector st{};

    state_vector& st_ref = st;
    state_vector const& st_cref = st;
    state_vector&& st_rref1 = state_vector{1, 2, 3, 4};
    state_vector&& st_rref2 = state_vector{1, 2, 3, 4};

    auto view_st = make_comp_state_view(st, hs);
    auto view_st_ref = make_comp_state_view(st_ref, hs);
    auto view_tmp = make_comp_state_view(state_vector{1, 2, 3, 4}, hs);
    auto view_rref = make_comp_state_view(std::move(st_rref1), hs);

    CHECK(std::is_same<decltype(view_st),
                       compressed_state_view<state_vector, true>>::value);
    CHECK(std::is_same<decltype(view_st_ref),
                       compressed_state_view<state_vector, true>>::value);
    CHECK(std::is_same<decltype(view_tmp),
                       compressed_state_view<state_vector, false>>::value);
    CHECK(std::is_same<decltype(view_rref),
                       compressed_state_view<state_vector, false>>::value);

    auto cview_st = make_const_comp_state_view(st, hs);
    auto cview_st_ref = make_const_comp_state_view(st_ref, hs);
    auto cview_st_cref = make_const_comp_state_view(st_cref, hs);
    auto cview_tmp =
        make_const_comp_state_view(state_vector{1, 2, 3, 4, 5, 6}, hs);
    auto cview_rref = make_const_comp_state_view(std::move(st_rref2), hs);

    CHECK(std::is_same<decltype(cview_st),
                       compressed_state_view<state_vector const, true>>::value);
    CHECK(std::is_same<decltype(cview_st_ref),
                       compressed_state_view<state_vector const, true>>::value);
    CHECK(std::is_same<decltype(cview_st_cref),
                       compressed_state_view<state_vector const, true>>::value);
    CHECK(
        std::is_same<decltype(cview_tmp),
                     compressed_state_view<state_vector const, false>>::value);
    CHECK(
        std::is_same<decltype(cview_rref),
                     compressed_state_view<state_vector const, false>>::value);
  }

  SECTION("loperator") {
    hs_type hs{make_space_fermion(0),
               make_space_fermion(1),
               make_space_boson(3, 1),
               make_space_boson(3, 2),
               make_space_boson(4, 3)};

    state_vector in(hs.dim());
    state_vector out(hs.dim());
    state_vector out_ref(hs.dim());

    auto view_in = compressed_state_view<state_vector const>(in, hs);
    auto view_out = compressed_state_view<state_vector>(out, hs);

    auto b1op = make_loperator(a_dag(1), hs);
    auto b2op = make_loperator(a_dag(2), hs);

    // in = |0>|0>|0>|0>
    std::fill(in.begin(), in.end(), 0);
    in[view_out.map_index(0)] = 1;

    // out = |0>|1>|0>|0>
    b1op(view_in, view_out);
    std::fill(out_ref.begin(), out_ref.end(), 0);
    out_ref[view_out.map_index(1 << 2)] = 1;
    CHECK_THAT(out, Catch::Matchers::Approx(out_ref));

    // out = |0>|1>|1>|0>
    in = out;
    b2op(view_in, view_out);
    std::fill(out_ref.begin(), out_ref.end(), 0);
    out_ref[view_out.map_index((1 << 2) + (1 << 4))] = 1;
    CHECK_THAT(out, Catch::Matchers::Approx(out_ref));

    // out = |0>|2>|1>|0>
    in = out;
    b1op(view_in, view_out);
    std::fill(out_ref.begin(), out_ref.end(), 0);
    out_ref[view_out.map_index((2 << 2) + (1 << 4))] = std::sqrt(2);
    CHECK_THAT(out, Catch::Matchers::Approx(out_ref));

    // out = |0>|2>|2>|0>
    in = out;
    b2op(view_in, view_out);
    std::fill(out_ref.begin(), out_ref.end(), 0);
    out_ref[view_out.map_index((2 << 2) + (2 << 4))] = 2;
    CHECK_THAT(out, Catch::Matchers::Approx(out_ref));

    // in = |0>|0>|1>|0>
    // out = |0>|0>|2>|0>
    std::fill(in.begin(), in.end(), 0);
    in[view_out.map_index(1 << 4)] = 1;
    b2op(view_in, view_out);
    std::fill(out_ref.begin(), out_ref.end(), 0);
    out_ref[view_out.map_index(2 << 4)] = std::sqrt(2);
    CHECK_THAT(out, Catch::Matchers::Approx(out_ref));
  }
}
