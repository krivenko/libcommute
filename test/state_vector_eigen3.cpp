/*******************************************************************************
 *
 * This file is part of libcommute, a quantum operator algebra DSL and
 * exact diagonalization toolkit for C++11/14/17.
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include "catch2/catch.hpp"

#include <libcommute/loperator/state_vector_eigen3.hpp>

#include <type_traits>
#include <vector>

using namespace libcommute;

template<typename VectorType, typename StateVector>
void test_functions(StateVector sv, typename VectorType::Scalar ref_val) {
  CHECK(get_element(sv, 1) == ref_val);
  update_add_element(sv, 1, 4.0);
  CHECK(get_element(sv, 1) == ref_val + 4.0);
  CHECK(zeros_like(sv) == VectorType::Zero(3));
  set_zeros(sv);
  CHECK(sv == VectorType::Zero(3));
}

TEST_CASE("Implementation of StateVector interface for Eigen3 types",
          "[state_vector]") {

  SECTION("VectorXd") {
    Eigen::VectorXd v(3); v << 1, 2, 3;

    CHECK(std::is_same<element_type_t<Eigen::VectorXd>, double>::value);
    CHECK(std::is_same<element_type_t<const Eigen::VectorXd>, double>::value);
    test_functions<Eigen::VectorXd>(v, 2.0);

    SECTION("foreach()") {
      Eigen::VectorXd v(4); v << 1, 2, 0, 4;
      foreach(v, [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("VectorXcd") {
    Eigen::VectorXcd v(3); v << 1, 2, 3;

    CHECK(std::is_same<element_type_t<decltype(v)>,
          std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<const decltype(v)>,
          std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(v, 2.0);

    SECTION("foreach()") {
      Eigen::VectorXcd v(4); v << 1, 2, 0, 4;
      foreach(v, [](int i, std::complex<double> a) {
        CHECK(double(i + 1) == a);
      });
    }
  }

  SECTION("Segment of VectorXd") {
    Eigen::VectorXd v(5); v << 1, 2, 3, 4, 5;
    auto vs = v.segment(1, 3);

    CHECK(std::is_same<element_type_t<decltype(vs)>, double>::value);
    CHECK(std::is_same<element_type_t<const decltype(vs)>, double>::value);
    test_functions<Eigen::VectorXd>(vs, 3.0);

    SECTION("foreach()") {
      Eigen::VectorXd v(6); v << 10, 1, 2, 0, 4, 5;
      auto vs = v.segment(1, 4);
      foreach(vs, [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("Segment of VectorXd") {
    Eigen::VectorXcd v(5); v << 1, 2, 3, 4, 5;
    auto vs = v.segment(1, 3);

    CHECK(std::is_same<element_type_t<decltype(vs)>,
                       std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<const decltype(vs)>,
                       std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(vs, 3.0);

    SECTION("foreach()") {
      Eigen::VectorXcd v(6); v << 10, 1, 2, 0, 4, 5;
      auto vs = v.segment(1, 4);
      foreach(vs, [](int i, std::complex<double> a) {
        CHECK(double(i + 1) == a);
      });
    }
  }

  SECTION("One column of MatrixXd") {
    Eigen::MatrixXd m(3, 3);
    m << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    auto v = m.col(1);

    CHECK(std::is_same<element_type_t<decltype(v)>, double>::value);
    CHECK(std::is_same<element_type_t<const decltype(v)>,
          double>::value);
    test_functions<Eigen::VectorXd>(v, 5.0);

    SECTION("foreach()") {
      Eigen::MatrixXd m(4, 4);
      m << 9, 1, 9, 9,
           9, 2, 9, 9,
           9, 0, 9, 9,
           9, 4, 9, 9;
      foreach(m.col(1), [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("One column of MatrixXcd") {
    Eigen::MatrixXcd m(3, 3);
    m << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    auto v = m.col(1);

    CHECK(std::is_same<element_type_t<decltype(v)>,
                       std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<const decltype(v)>,
                       std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(v, 5.0);

    SECTION("foreach()") {
      Eigen::MatrixXcd m(4, 4);
      m << 9, 1, 9, 9,
           9, 2, 9, 9,
           9, 0, 9, 9,
           9, 4, 9, 9;
      foreach(m.col(1), [](int i, std::complex<double> a) {
        CHECK(double(i + 1) == a);
      });
    }
  }

  SECTION("Part of a column of MatrixXd") {
    Eigen::MatrixXd m(5, 5);
    m << 1,   2,  3,  4,  5,
         6,   7,  8,  9, 10,
         11, 12, 13, 14, 15,
         16, 17, 18, 19, 20,
         21, 22, 23, 24, 25;
    auto v = m.block(1, 2, 3, 1);

    CHECK(std::is_same<element_type_t<decltype(v)>, double>::value);
    CHECK(std::is_same<element_type_t<const decltype(v)>, double>::value);
    test_functions<Eigen::VectorXd>(v, 13.0);

    SECTION("foreach()") {
      Eigen::MatrixXd m(5, 5);
      m << 1,   2, 5,  4,  5,
           6,   7, 1,  9, 10,
           11, 12, 2, 14, 15,
           16, 17, 0, 19, 20,
           21, 22, 4, 24, 25;
      auto v = m.block(1, 2, 4, 1);
      foreach(v, [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("Part of a column of MatrixXcd") {
    Eigen::MatrixXcd m(5, 5);
    m << 1,   2,  3,  4,  5,
         6,   7,  8,  9, 10,
         11, 12, 13, 14, 15,
         16, 17, 18, 19, 20,
         21, 22, 23, 24, 25;
    auto v = m.block(1, 2, 3, 1);

    CHECK(std::is_same<element_type_t<decltype(v)>,
                       std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<const decltype(v)>,
                       std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(v, 13.0);

    SECTION("foreach()") {
      Eigen::MatrixXcd m(5, 5);
      m << 1,   2, 5,  4,  5,
           6,   7, 1,  9, 10,
           11, 12, 2, 14, 15,
           16, 17, 0, 19, 20,
           21, 22, 4, 24, 25;
      auto v = m.block(1, 2, 4, 1);
      foreach(v, [](int i, std::complex<double> a) {
        CHECK(double(i + 1) == a);
      });
    }
  }

  SECTION("Map<VectorXd>") {
    std::vector<double> v{0, 1, 2, 3, 4};
    auto m = Eigen::Map<Eigen::VectorXd>(v.data() + 1, 3, 1);

    CHECK(std::is_same<element_type_t<decltype(m)>, double>::value);
    CHECK(std::is_same<element_type_t<const decltype(m)>, double>::value);
    test_functions<Eigen::VectorXd>(m, 2.0);

    SECTION("foreach()") {
      std::vector<double> v{5, 1, 2, 0, 4, 5};
      auto m = Eigen::Map<Eigen::VectorXd>(v.data() + 1, 4, 1);
      foreach(m, [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("Map<VectorXcd>") {
    std::vector<std::complex<double>> v{0, 1, 2, 3, 4};
    auto m = Eigen::Map<Eigen::VectorXcd>(v.data() + 1, 3, 1);

    CHECK(std::is_same<element_type_t<decltype(m)>,
                       std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<const decltype(m)>,
                       std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(m, 2.0);

    SECTION("foreach()") {
      std::vector<std::complex<double>> v{5, 1, 2, 0, 4, 5};
      auto m = Eigen::Map<Eigen::VectorXcd>(v.data() + 1, 4, 1);
      foreach(m, [](int i, std::complex<double> a) {
        CHECK(double(i + 1) == a);
      });
    }
  }
}
