/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <catch.hpp>

#include "my_complex.hpp"
#include "print_matcher.hpp"

#include <libcommute/expression/factories.hpp>

#include <complex>
#include <type_traits>

using namespace libcommute;

using monomial_t = monomial<int, std::string>;

template <typename E, typename S, typename... Generators>
void check_monomial(E const& expr, S ref_coeff, Generators&&... generators) {
  monomial_t ref_monomial(std::forward<Generators>(generators)...);

  CHECK(expr.get_monomials().size() == 1);
  CHECK((expr.get_monomials().begin()->first) == ref_monomial);
  CHECK((expr.get_monomials().begin()->second) == ref_coeff);
}

TEST_CASE("Factory functions", "[factories]") {

  using namespace static_indices;

  SECTION("make_complex()") {
    auto expr = make_complex(4.0 * c_dag(1, "up") * c(2, "dn") + 1.0);
    CHECK(std::is_same<
          decltype(expr),
          expression<std::complex<double>, int, std::string>>::value);
    CHECK(expr == std::complex<double>(4, 0) * c_dag(1, "up") * c(2, "dn") +
                      std::complex<double>(1, 0));
  }

  SECTION("my_complex") {
    SECTION("fermion") {
      auto c_dag_1_up = c_dag<my_complex>(1, "up");
      check_monomial(c_dag_1_up, my_complex{1, 0}, make_fermion(true, 1, "up"));

      auto c_2_dn = c<my_complex>(2, "dn");
      check_monomial(c_2_dn, my_complex{1, 0}, make_fermion(false, 2, "dn"));

      auto n_1_dn = n<my_complex>(1, "dn");
      check_monomial(n_1_dn,
                     my_complex{1, 0},
                     make_fermion(true, 1, "dn"),
                     make_fermion(false, 1, "dn"));
    }
    SECTION("boson") {
      auto a_dag_x = a_dag<my_complex>(0, "x");
      check_monomial(a_dag_x, my_complex{1, 0}, make_boson(true, 0, "x"));

      auto a_y = a<my_complex>(0, "y");
      check_monomial(a_y, my_complex{1, 0}, make_boson(false, 0, "y"));
    }

    SECTION("spin-1/2") {
      auto S_p_0_x = S_p<my_complex>(0, "x");
      check_monomial(S_p_0_x,
                     my_complex{1, 0},
                     make_spin(spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<my_complex>(0, "x");
      check_monomial(S_m_0_x,
                     my_complex{1, 0},
                     make_spin(spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<my_complex>(0, "x");
      check_monomial(S_z_0_x,
                     my_complex{1, 0},
                     make_spin(spin_component::z, 0, "x"));
    }
    SECTION("spin-1") {
      auto S_p_0_x = S_p<3, my_complex>(0, "x");
      check_monomial(S_p_0_x,
                     my_complex{1, 0},
                     make_spin(1.0, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<3, my_complex>(0, "x");
      check_monomial(S_m_0_x,
                     my_complex{1, 0},
                     make_spin(1.0, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<3, my_complex>(0, "x");
      check_monomial(S_z_0_x,
                     my_complex{1, 0},
                     make_spin(1.0, spin_component::z, 0, "x"));
    }
    SECTION("spin-3/2") {
      auto S_p_0_x = S_p<4, my_complex>(0, "x");
      check_monomial(S_p_0_x,
                     my_complex{1, 0},
                     make_spin(3.0 / 2, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<4, my_complex>(0, "x");
      check_monomial(S_m_0_x,
                     my_complex{1, 0},
                     make_spin(3.0 / 2, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<4, my_complex>(0, "x");
      check_monomial(S_z_0_x,
                     my_complex{1, 0},
                     make_spin(3.0 / 2, spin_component::z, 0, "x"));
    }
  }

  SECTION("real") {
    SECTION("fermion") {
      auto c_dag_1_up = c_dag(1, "up");
      check_monomial(c_dag_1_up, 1, make_fermion(true, 1, "up"));

      auto c_2_dn = c(2, "dn");
      check_monomial(c_2_dn, 1, make_fermion(false, 2, "dn"));

      auto n_1_dn = n(1, "dn");
      check_monomial(n_1_dn,
                     1,
                     make_fermion(true, 1, "dn"),
                     make_fermion(false, 1, "dn"));
    }
    SECTION("boson") {
      auto a_dag_x = a_dag(0, "x");
      check_monomial(a_dag_x, 1, make_boson(true, 0, "x"));

      auto a_y = a(0, "y");
      check_monomial(a_y, 1, make_boson(false, 0, "y"));
    }
    SECTION("spin-1/2") {
      auto S_p_0_x = S_p(0, "x");
      check_monomial(S_p_0_x, 1, make_spin(spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m(0, "x");
      check_monomial(S_m_0_x, 1, make_spin(spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z(0, "x");
      check_monomial(S_z_0_x, 1, make_spin(spin_component::z, 0, "x"));
    }
    SECTION("spin-1") {
      auto S_p_0_x = S_p<3>(0, "x");
      check_monomial(S_p_0_x, 1, make_spin(1.0, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<3>(0, "x");
      check_monomial(S_m_0_x, 1, make_spin(1.0, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<3>(0, "x");
      check_monomial(S_z_0_x, 1, make_spin(1.0, spin_component::z, 0, "x"));
    }
    SECTION("spin-3/2") {
      auto S_p_0_x = S_p<4>(0, "x");
      check_monomial(S_p_0_x,
                     1,
                     make_spin(3.0 / 2, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<4>(0, "x");
      check_monomial(S_m_0_x,
                     1,
                     make_spin(3.0 / 2, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<4>(0, "x");
      check_monomial(S_z_0_x, 1, make_spin(3.0 / 2, spin_component::z, 0, "x"));
    }
  }

  SECTION("complex") {
    std::complex<double> const c1(1.0);
    std::complex<double> const I(0, 1.0);

    SECTION("fermion") {
      auto c_dag_1_up = make_complex(c_dag(1, "up"));
      check_monomial(c_dag_1_up, c1, make_fermion(true, 1, "up"));

      auto c_2_dn = make_complex(c(2, "dn"));
      check_monomial(c_2_dn, c1, make_fermion(false, 2, "dn"));

      auto n_1_dn = make_complex(n(1, "dn"));
      check_monomial(n_1_dn,
                     c1,
                     make_fermion(true, 1, "dn"),
                     make_fermion(false, 1, "dn"));
    }
    SECTION("boson") {
      auto a_dag_x = make_complex(a_dag(0, "x"));
      check_monomial(a_dag_x, c1, make_boson(true, 0, "x"));

      auto a_y = make_complex(a(0, "y"));
      check_monomial(a_y, c1, make_boson(false, 0, "y"));
    }
    SECTION("spin-1/2") {
      auto S_p_0_x = make_complex(S_p(0, "x"));
      check_monomial(S_p_0_x, c1, make_spin(spin_component::plus, 0, "x"));
      auto S_m_0_x = make_complex(S_m(0, "x"));
      check_monomial(S_m_0_x, c1, make_spin(spin_component::minus, 0, "x"));
      auto S_z_0_x = make_complex(S_z(0, "x"));
      check_monomial(S_z_0_x, c1, make_spin(spin_component::z, 0, "x"));

      CHECK(S_p_0_x == S_x(0, "x") + I * S_y(0, "x"));
      CHECK(S_m_0_x == S_x(0, "x") - I * S_y(0, "x"));
    }
    SECTION("spin-1") {
      auto S_p_0_x = make_complex(S_p<3>(0, "x"));
      check_monomial(S_p_0_x, c1, make_spin(1.0, spin_component::plus, 0, "x"));
      auto S_m_0_x = make_complex(S_m<3>(0, "x"));
      check_monomial(S_m_0_x,
                     c1,
                     make_spin(1.0, spin_component::minus, 0, "x"));
      auto S_z_0_x = make_complex(S_z<3>(0, "x"));
      check_monomial(S_z_0_x, c1, make_spin(1.0, spin_component::z, 0, "x"));

      CHECK(S_p_0_x == S_x<3>(0, "x") + I * S_y<3>(0, "x"));
      CHECK(S_m_0_x == S_x<3>(0, "x") - I * S_y<3>(0, "x"));
    }
    SECTION("spin-3/2") {
      auto S_p_0_x = make_complex(S_p<4>(0, "x"));
      check_monomial(S_p_0_x,
                     c1,
                     make_spin(3.0 / 2, spin_component::plus, 0, "x"));
      auto S_m_0_x = make_complex(S_m<4>(0, "x"));
      check_monomial(S_m_0_x,
                     c1,
                     make_spin(3.0 / 2, spin_component::minus, 0, "x"));
      auto S_z_0_x = make_complex(S_z<4>(0, "x"));
      check_monomial(S_z_0_x,
                     c1,
                     make_spin(3.0 / 2, spin_component::z, 0, "x"));

      CHECK(S_p_0_x == S_x<4>(0, "x") + I * S_y<4>(0, "x"));
      CHECK(S_m_0_x == S_x<4>(0, "x") - I * S_y<4>(0, "x"));
    }
  }
}
