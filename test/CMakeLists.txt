#
# This file is part of libcommute, a quantum operator algebra DSL and
# exact diagonalization toolkit for C++11/14/17.
#
# Copyright (C) 2016-2025 Igor Krivenko
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Build Catch2 main object file
add_library(catch2 OBJECT catch2/catch2-main.cpp)
# Exclude Catch2 from static analysis
if(STATIC_ANALYSIS)
  set_property(TARGET catch2 PROPERTY CXX_CLANG_TIDY)
  set_property(TARGET catch2 PROPERTY CXX_CPPCHECK)
endif(STATIC_ANALYSIS)

# Make catch.hpp visible as a system header
include_directories(BEFORE SYSTEM catch2)

set(TESTS
  utility
  metafunctions
  generator
  monomial
  scalar_traits
  factories
  expression
  expression.addition
  expression.subtraction
  expression.multiplication
  expression.const.addition
  expression.const.subtraction
  expression.const.multiplication
  expression.mixed_arithmetics
  hc
  commutators
  new_algebra
  state_vector
  elementary_space
  hilbert_space
  monomial_action
  monomial_action_fermion
  monomial_action_boson
  monomial_action_spin
  loperator
  new_algebra.loperator
  disjoint_sets
  sparse_state_vector
  space_partition
  mapped_basis_view
  n_fermion_sector_view
  n_fermion_multisector_view
)

# Build C++ unit tests
foreach(t ${TESTS})
  set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
  add_executable(${t} ${s})
  target_link_libraries(${t} PRIVATE libcommute catch2)
  add_test(NAME ${t} COMMAND ${t})
endforeach()

set(CXX17_TESTS
  dyn_indices
  generator_dyn
  factories_dyn
  elementary_space_dyn
)

# C++17 tests
if("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
  message(STATUS "Enabling C++17 tests")
  foreach(t ${CXX17_TESTS})
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    add_executable(${t} ${s})
    target_compile_features(${t} PRIVATE cxx_std_17)
    target_link_libraries(${t} PRIVATE libcommute catch2)
    add_test(NAME ${t} COMMAND ${t})
  endforeach()
else()
  message(STATUS "Disabling C++17 tests")
endif()

# Tests using Eigen 3
if(Eigen3_FOUND)
  message(STATUS "Enabling Eigen 3 tests")
  set(EIGEN3_TESTS state_vector_eigen3)
  foreach(t ${EIGEN3_TESTS})
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    add_executable(${t} ${s})
    target_link_libraries(${t} PRIVATE libcommute catch2)
    link_to_eigen3(${t})
    add_test(NAME ${t} COMMAND ${t})
  endforeach()
endif(Eigen3_FOUND)
