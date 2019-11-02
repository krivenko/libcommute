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

#include "catch2/catch.hpp"

#include "int_complex.hpp"
#include "print_matcher.hpp"

#include <libcommute/expression/factories_dyn.hpp>

using namespace libcommute;
using namespace libcommute::dynamic_indices;

using monomial_t = monomial<dyn_indices>;

template<typename E, typename S, typename... Generators>
void check_monomial(E const& expr, S ref_coeff, Generators&&... generators) {
  monomial_t ref_monomial(std::forward<Generators>(generators)...);

  CHECK(expr.get_monomials().size() == 1);
  CHECK((expr.get_monomials().begin()->first) == ref_monomial);
  CHECK((expr.get_monomials().begin()->second) == ref_coeff);
}

TEST_CASE("Factory functions", "[factories]") {

  using namespace dynamic_indices;

  SECTION("int_complex") {
    SECTION("fermion") {
      auto c_dag_1_up = c_dag<int_complex>(1, "up");
      check_monomial(c_dag_1_up,
                     int_complex{1, 0},
                     make_fermion(true, 1, "up"));

      auto c_2_dn = c<int_complex>(2, "dn");
      check_monomial(c_2_dn,
                     int_complex{1, 0},
                     make_fermion(false, 2, "dn"));

      auto n_1_dn = n<int_complex>(1, "dn");
      check_monomial(n_1_dn,
                    int_complex{1, 0},
                    make_fermion(true, 1, "dn"),
                    make_fermion(false, 1, "dn")
                    );
    }
    SECTION("boson") {
      auto a_dag_x = a_dag<int_complex>(0, "x");
      check_monomial(a_dag_x, int_complex{1, 0}, make_boson(true, 0, "x"));

      auto a_y = a<int_complex>(0, "y");
      check_monomial(a_y, int_complex{1, 0}, make_boson(false, 0, "y"));
    }

    SECTION("spin-1/2") {
      auto S_p_0_x = S_p<int_complex>(0, "x");
      check_monomial(S_p_0_x,
                     int_complex{1, 0},
                     make_spin(spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<int_complex>(0, "x");
      check_monomial(S_m_0_x,
                     int_complex{1, 0},
                     make_spin(spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<int_complex>(0, "x");
      check_monomial(S_z_0_x,
                     int_complex{1, 0},
                     make_spin(spin_component::z, 0, "x"));
    }
    SECTION("spin-1") {
      auto S_p_0_x = S_p<3, int_complex>(0, "x");
      check_monomial(S_p_0_x,
                     int_complex{1, 0},
                     make_spin(1.0, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<3, int_complex>(0, "x");
      check_monomial(S_m_0_x,
                     int_complex{1, 0},
                     make_spin(1.0, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<3, int_complex>(0, "x");
      check_monomial(S_z_0_x,
                     int_complex{1, 0},
                     make_spin(1.0, spin_component::z, 0, "x"));
    }
    SECTION("spin-3/2") {
      auto S_p_0_x = S_p<4, int_complex>(0, "x");
      check_monomial(S_p_0_x,
                     int_complex{1, 0},
                     make_spin(3.0/2, spin_component::plus, 0, "x"));
      auto S_m_0_x = S_m<4, int_complex>(0, "x");
      check_monomial(S_m_0_x,
                     int_complex{1, 0},
                     make_spin(3.0/2, spin_component::minus, 0, "x"));
      auto S_z_0_x = S_z<4, int_complex>(0, "x");
      check_monomial(S_z_0_x,
                     int_complex{1, 0},
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

  SECTION("Mixed indices") {
    using namespace real;

    auto expr = c_dag(0, "up")*c(1, "dn") +
                a_dag("x") * n() + a("y") * n() + S_p()*S_m();
    CHECK_THAT(expr, Prints<decltype(expr)>(
      "1*C+(0,up)C(1,dn) + 1*S+()S-() + 1*C+()C()A+(x) + 1*C+()C()A(y)")
    );
  }
}
