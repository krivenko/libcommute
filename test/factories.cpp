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

#include "catch2/catch.hpp"

#include "my_complex.hpp"

#include <libcommute/expression/factories.hpp>

using namespace libcommute;

using monomial_t = monomial<int, std::string>;

template<typename E, typename S, typename... Generators>
void check_monomial(E const& expr, S ref_coeff, Generators&&... generators) {
  monomial_t ref_monomial(std::forward<Generators>(generators)...);

  CHECK(expr.get_monomials().size() == 1);
  CHECK((expr.get_monomials().begin()->first) == ref_monomial);
  CHECK((expr.get_monomials().begin()->second) == ref_coeff);
}

TEST_CASE("Factory functions", "[factories]") {

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
                    make_fermion(false, 1, "dn")
                    );
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
                     make_spin(3.0/2, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<4, my_complex>(0, "x");
      check_monomial(S_m_0_x,
                     my_complex{1, 0},
                     make_spin(3.0/2, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<4, my_complex>(0, "x");
      check_monomial(S_z_0_x,
                     my_complex{1, 0},
                     make_spin(3.0/2, spin_component::z, 0, "x"));
    }
  }

  SECTION("real") {
    SECTION("fermion") {
      using real::c_dag;
      using real::c;
      using real::n;

      auto c_dag_1_up = c_dag(1, "up");
      check_monomial(c_dag_1_up, 1, make_fermion(true, 1, "up"));

      auto c_2_dn = c(2, "dn");
      check_monomial(c_2_dn, 1, make_fermion(false, 2, "dn"));

      auto n_1_dn = n(1, "dn");
      check_monomial(n_1_dn,
                    1,
                    make_fermion(true, 1, "dn"),
                    make_fermion(false, 1, "dn")
                    );
    }
    SECTION("boson") {
      using real::a_dag;
      using real::a;

      auto a_dag_x = a_dag(0, "x");
      check_monomial(a_dag_x, 1, make_boson(true, 0, "x"));

      auto a_y = a(0, "y");
      check_monomial(a_y, 1, make_boson(false, 0, "y"));
    }
    SECTION("spin-1/2") {
      using real::S_p;
      using real::S_m;
      using real::S_z;

      auto S_p_0_x = S_p(0, "x");
      check_monomial(S_p_0_x, 1, make_spin(spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m(0, "x");
      check_monomial(S_m_0_x, 1, make_spin(spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z(0, "x");
      check_monomial(S_z_0_x, 1, make_spin(spin_component::z, 0, "x"));
    }
    SECTION("spin-1") {
      using real::S_p;
      using real::S_m;
      using real::S_z;

      auto S_p_0_x = S_p<3>(0, "x");
      check_monomial(S_p_0_x, 1, make_spin(1.0, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<3>(0, "x");
      check_monomial(S_m_0_x, 1, make_spin(1.0, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<3>(0, "x");
      check_monomial(S_z_0_x, 1, make_spin(1.0, spin_component::z, 0, "x"));
    }
    SECTION("spin-3/2") {
      using real::S_p;
      using real::S_m;
      using real::S_z;

      auto S_p_0_x = S_p<4>(0, "x");
      check_monomial(S_p_0_x,
                     1,
                     make_spin(3.0/2, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<4>(0, "x");
      check_monomial(S_m_0_x,
                     1,
                     make_spin(3.0/2, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<4>(0, "x");
      check_monomial(S_z_0_x,
                     1,
                     make_spin(3.0/2, spin_component::z, 0, "x"));
    }
  }

  SECTION("complex") {
    const std::complex<double> c1(1.0);

    SECTION("fermion") {
      using complex::c_dag;
      using complex::c;
      using complex::n;

      auto c_dag_1_up = c_dag(1, "up");
      check_monomial(c_dag_1_up, c1, make_fermion(true, 1, "up"));

      auto c_2_dn = c(2, "dn");
      check_monomial(c_2_dn, c1, make_fermion(false, 2, "dn"));

      auto n_1_dn = n(1, "dn");
      check_monomial(n_1_dn,
                    c1,
                    make_fermion(true, 1, "dn"),
                    make_fermion(false, 1, "dn")
                    );
    }
    SECTION("boson") {
      using complex::a_dag;
      using complex::a;

      auto a_dag_x = a_dag(0, "x");
      check_monomial(a_dag_x, c1, make_boson(true, 0, "x"));

      auto a_y = a(0, "y");
      check_monomial(a_y, c1, make_boson(false, 0, "y"));
    }
    SECTION("spin-1/2") {
      using complex::S_p;
      using complex::S_m;
      using complex::S_z;

      auto S_p_0_x = S_p(0, "x");
      check_monomial(S_p_0_x, c1, make_spin(spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m(0, "x");
      check_monomial(S_m_0_x, c1, make_spin(spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z(0, "x");
      check_monomial(S_z_0_x, c1, make_spin(spin_component::z, 0, "x"));
    }
    SECTION("spin-1") {
      using complex::S_p;
      using complex::S_m;
      using complex::S_z;

      auto S_p_0_x = S_p<3>(0, "x");
      check_monomial(S_p_0_x, c1, make_spin(1.0, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<3>(0, "x");
      check_monomial(S_m_0_x, c1, make_spin(1.0, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<3>(0, "x");
      check_monomial(S_z_0_x, c1, make_spin(1.0, spin_component::z, 0, "x"));
    }
    SECTION("spin-3/2") {
      using complex::S_p;
      using complex::S_m;
      using complex::S_z;

      auto S_p_0_x = S_p<4>(0, "x");
      check_monomial(S_p_0_x,
                     c1,
                     make_spin(3.0/2, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<4>(0, "x");
      check_monomial(S_m_0_x,
                     c1,
                     make_spin(3.0/2, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<4>(0, "x");
      check_monomial(S_z_0_x,
                     c1,
                     make_spin(3.0/2, spin_component::z, 0, "x"));
    }
  }
}
