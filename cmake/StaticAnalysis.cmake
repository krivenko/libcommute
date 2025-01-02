#
# This file is part of libcommute, a quantum operator algebra DSL and
# exact diagonalization toolkit for C++11/14/17.
#
# Copyright (C) 2016-2025 Igor Krivenko
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#
# clang-tidy
#

find_program(CLANG_TIDY_EXECUTABLE
             NAMES "clang-tidy" REQUIRED
             DOC "Path to the clang-tidy executable")
mark_as_advanced(CLANG_TIDY_EXECUTABLE)

message(STATUS "Using clang-tidy: ${CLANG_TIDY_EXECUTABLE}")
set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")

#
# cppcheck
#

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.10)
  find_program(CPPCHECK_EXECUTABLE
               NAMES "cppcheck" REQUIRED
               DOC "Path to the cppcheck executable")
  set(CPPCHECK_EXTRA_FLAGS "" CACHE STRING
      "Semicolon-separated list of additional flags to be passed to cppcheck")

  mark_as_advanced(CPPCHECK_EXECUTABLE)
  mark_as_advanced(CPPCHECK_EXTRA_FLAGS)

  message(STATUS "Using cppcheck: ${CPPCHECK_EXECUTABLE}")
  set(CMAKE_CXX_CPPCHECK
      "${CPPCHECK_EXECUTABLE}"
      "--enable=warning,style,performance,portability"
      "--std=c++17"
      "--template=gcc"
      "--inline-suppr"
      "--verbose"
      "--force"
      "--quiet"
  )

  if(CPPCHECK_EXTRA_FLAGS)
    message(STATUS "Extra flags passed to cppcheck: ${CPPCHECK_EXTRA_FLAGS}")
    list(APPEND CMAKE_CXX_CPPCHECK ${CPPCHECK_EXTRA_FLAGS})
  endif(CPPCHECK_EXTRA_FLAGS)
else(CMAKE_VERSION VERSION_GREATER_EQUAL 3.10)
  message(WARNING "CMake >= 3.10 is required to run cppcheck")
endif(CMAKE_VERSION VERSION_GREATER_EQUAL 3.10)
