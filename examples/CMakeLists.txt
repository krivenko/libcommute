#
# This file is part of libcommute, a quantum operator algebra DSL and
# exact diagonalization toolkit for C++11/14/17.
#
# Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

set(EXAMPLES
  hubbard_holstein
  exchange_matrix
  heisenberg_chain
  parametric_qoperator
  virasoro_algebra
)

# Build C++11 examples
foreach(e ${EXAMPLES})
  set(s ${CMAKE_CURRENT_SOURCE_DIR}/${e}.cpp)
  add_executable(${e} ${s})
  target_link_libraries(${e} PRIVATE libcommute)
endforeach()
