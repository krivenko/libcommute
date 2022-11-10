#
# This file is part of libcommute, a quantum operator algebra DSL and
# exact diagonalization toolkit for C++11/14/17.
#
# Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

set(EXAMPLES
  holstein
  gamma.expression
  gamma.loperator
  partition
  hubbard_holstein
  exchange_matrix
  kanamori
  jaynes_cummings
  heisenberg_chain
  parametric_loperator
  virasoro_algebra
)

# Build C++11 examples
foreach(e ${EXAMPLES})
  set(s ${CMAKE_CURRENT_SOURCE_DIR}/${e}.cpp)
  add_executable(${e} ${s})
  target_link_libraries(${e} PRIVATE libcommute)
endforeach()

# Examples using Eigen 3
if(Eigen3_FOUND)
  message(STATUS "Enabling Eigen 3 examples")
  set(EIGEN3_EXAMPLES n_fermion_sectors hubbard_holstein_1d)
  foreach(e ${EIGEN3_EXAMPLES})
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/${e}.cpp)
    add_executable(${e} ${s})
    target_link_libraries(${e} PRIVATE libcommute)
    link_to_eigen3(${e})
  endforeach()
endif(Eigen3_FOUND)
