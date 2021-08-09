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

#include <catch.hpp>

#include <libcommute/loperator/state_vector_eigen3.hpp>

#include <type_traits>
#include <vector>

using namespace libcommute;

template <typename VectorType, typename StateVector>
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
    // clang-format off
    // cppcheck-suppress constStatement
    Eigen::VectorXd v(3); v << 1, 2, 3;
    // clang-format on

    CHECK(std::is_same<element_type_t<Eigen::VectorXd>, double>::value);
    CHECK(std::is_same<element_type_t<Eigen::VectorXd const>, double>::value);
    test_functions<Eigen::VectorXd>(v, 2.0);

    SECTION("foreach()") {
      // clang-format off
      // cppcheck-suppress constStatement
      Eigen::VectorXd v_fe(4); v_fe << 1, 2, 0, 4;
      // clang-format on
      foreach(v_fe, [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("VectorXcd") {
    // clang-format off
    // cppcheck-suppress constStatement
    Eigen::VectorXcd v(3); v << 1, 2, 3;
    // clang-format on

    CHECK(
        std::is_same<element_type_t<decltype(v)>, std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<decltype(v) const>,
                       std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(v, 2.0);

    SECTION("foreach()") {
      // clang-format off
      // cppcheck-suppress constStatement
      Eigen::VectorXcd v_fe(4); v_fe << 1, 2, 0, 4;
      // clang-format on
      foreach(v_fe,
              [](int i, std::complex<double> a) { CHECK(double(i + 1) == a); });
    }
  }

  SECTION("Segment of VectorXd") {
    // clang-format off
    // cppcheck-suppress constStatement
    Eigen::VectorXd v(5); v << 1, 2, 3, 4, 5;
    // clang-format on
    auto vs = v.segment(1, 3);

    CHECK(std::is_same<element_type_t<decltype(vs)>, double>::value);
    CHECK(std::is_same<element_type_t<decltype(vs) const>, double>::value);
    test_functions<Eigen::VectorXd>(vs, 3.0);

    SECTION("foreach()") {
      // clang-format off
      // cppcheck-suppress constStatement
      Eigen::VectorXd v_fe(6); v_fe << 10, 1, 2, 0, 4, 5;
      // clang-format on
      auto vs_fe = v_fe.segment(1, 4);
      foreach(vs_fe, [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("Segment of VectorXd") {
    // clang-format off
    // cppcheck-suppress constStatement
    Eigen::VectorXcd v(5); v << 1, 2, 3, 4, 5;
    // clang-format on
    auto vs = v.segment(1, 3);

    CHECK(std::is_same<element_type_t<decltype(vs)>,
                       std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<decltype(vs) const>,
                       std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(vs, 3.0);

    SECTION("foreach()") {
      // clang-format off
      // cppcheck-suppress constStatement
      Eigen::VectorXcd v_fe(6); v_fe << 10, 1, 2, 0, 4, 5;
      // clang-format on
      auto vs_fe = v_fe.segment(1, 4);
      foreach(vs_fe,
              [](int i, std::complex<double> a) { CHECK(double(i + 1) == a); });
    }
  }

  SECTION("One column of MatrixXd") {
    Eigen::MatrixXd m(3, 3);
    // clang-format off
    m << 1, 2, 3,
         4, 5, 6,
         // cppcheck-suppress constStatement
         7, 8, 9;
    // clang-format on
    auto v = m.col(1);

    CHECK(std::is_same<element_type_t<decltype(v)>, double>::value);
    CHECK(std::is_same<element_type_t<decltype(v) const>, double>::value);
    test_functions<Eigen::VectorXd>(v, 5.0);

    SECTION("foreach()") {
      Eigen::MatrixXd m_fe(4, 4);
      // clang-format off
      m_fe << 9, 1, 9, 9,
              9, 2, 9, 9,
              9, 0, 9, 9,
              // cppcheck-suppress constStatement
              9, 4, 9, 9;
      // clang-format on
      foreach(m_fe.col(1), [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("One column of MatrixXcd") {
    Eigen::MatrixXcd m(3, 3);
    // clang-format off
    m << 1, 2, 3,
         4, 5, 6,
         // cppcheck-suppress constStatement
         7, 8, 9;
    // clang-format on
    auto v = m.col(1);

    CHECK(
        std::is_same<element_type_t<decltype(v)>, std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<decltype(v) const>,
                       std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(v, 5.0);

    SECTION("foreach()") {
      Eigen::MatrixXcd m_fe(4, 4);
      // clang-format off
      m_fe << 9, 1, 9, 9,
              9, 2, 9, 9,
              9, 0, 9, 9,
              // cppcheck-suppress constStatement
              9, 4, 9, 9;
      // clang-format on
      foreach(m_fe.col(1),
              [](int i, std::complex<double> a) { CHECK(double(i + 1) == a); });
    }
  }

  SECTION("Part of a column of MatrixXd") {
    Eigen::MatrixXd m(5, 5);
    // clang-format off
    m << 1,   2,  3,  4,  5,
         6,   7,  8,  9, 10,
         11, 12, 13, 14, 15,
         16, 17, 18, 19, 20,
         // cppcheck-suppress constStatement
         21, 22, 23, 24, 25;
    // clang-format on
    auto v = m.block(1, 2, 3, 1);

    CHECK(std::is_same<element_type_t<decltype(v)>, double>::value);
    CHECK(std::is_same<element_type_t<decltype(v) const>, double>::value);
    test_functions<Eigen::VectorXd>(v, 13.0);

    SECTION("foreach()") {
      Eigen::MatrixXd m_fe(5, 5);
      // clang-format off
      m_fe << 1,   2, 5,  4,  5,
              6,   7, 1,  9, 10,
              11, 12, 2, 14, 15,
              16, 17, 0, 19, 20,
              // cppcheck-suppress constStatement
              21, 22, 4, 24, 25;
      // clang-format on
      auto v_fe = m_fe.block(1, 2, 4, 1);
      foreach(v_fe, [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("Part of a column of MatrixXcd") {
    Eigen::MatrixXcd m(5, 5);
    // clang-format off
    m << 1,   2,  3,  4,  5,
         6,   7,  8,  9, 10,
         11, 12, 13, 14, 15,
         16, 17, 18, 19, 20,
         // cppcheck-suppress constStatement
         21, 22, 23, 24, 25;
    // clang-format on
    auto v = m.block(1, 2, 3, 1);

    CHECK(
        std::is_same<element_type_t<decltype(v)>, std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<decltype(v) const>,
                       std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(v, 13.0);

    SECTION("foreach()") {
      Eigen::MatrixXcd m_fe(5, 5);
      // clang-format off
      m_fe << 1,   2, 5,  4,  5,
              6,   7, 1,  9, 10,
              11, 12, 2, 14, 15,
              16, 17, 0, 19, 20,
              // cppcheck-suppress constStatement
              21, 22, 4, 24, 25;
      // clang-format on
      auto v_fe = m_fe.block(1, 2, 4, 1);
      foreach(v_fe,
              [](int i, std::complex<double> a) { CHECK(double(i + 1) == a); });
    }
  }

  SECTION("Map<VectorXd>") {
    std::vector<double> v{0, 1, 2, 3, 4};
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    auto m = Eigen::Map<Eigen::VectorXd>(v.data() + 1, 3, 1);

    CHECK(std::is_same<element_type_t<decltype(m)>, double>::value);
    CHECK(std::is_same<element_type_t<decltype(m) const>, double>::value);
    test_functions<Eigen::VectorXd>(m, 2.0);

    SECTION("foreach()") {
      std::vector<double> v_fe{5, 1, 2, 0, 4, 5};
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      auto m_fe = Eigen::Map<Eigen::VectorXd>(v_fe.data() + 1, 4, 1);
      foreach(m_fe, [](int i, double a) { CHECK(i + 1 == a); });
    }
  }

  SECTION("Map<VectorXcd>") {
    std::vector<std::complex<double>> v{0, 1, 2, 3, 4};
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    auto m = Eigen::Map<Eigen::VectorXcd>(v.data() + 1, 3, 1);

    CHECK(
        std::is_same<element_type_t<decltype(m)>, std::complex<double>>::value);
    CHECK(std::is_same<element_type_t<decltype(m) const>,
                       std::complex<double>>::value);
    test_functions<Eigen::VectorXcd>(m, 2.0);

    SECTION("foreach()") {
      std::vector<std::complex<double>> v_fe{5, 1, 2, 0, 4, 5};
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      auto m_fe = Eigen::Map<Eigen::VectorXcd>(v_fe.data() + 1, 4, 1);
      foreach(m_fe,
              [](int i, std::complex<double> a) { CHECK(double(i + 1) == a); });
    }
  }
}
