#
# This file is part of libcommute, a C++11/14/17 header-only library
# to manipulate polynomial expressions with quantum-mechanical operators.
#
# Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Start configuration
cmake_minimum_required(VERSION 3.8.0 FATAL_ERROR)

# Project name and version
set(LIBCOMMUTE_VERSION 0.1)
project(libcommute VERSION ${LIBCOMMUTE_VERSION} LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

message(STATUS "Configuring ${PROJECT_NAME} version ${PROJECT_VERSION}")

# Prohibit in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. Please make a separate build directory.")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

# CMake options
option(Tests "Build unit tests" ON)
option(Examples "Build examples" ON)

# Install C++ headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include
        DESTINATION ${CMAKE_INSTALL_PREFIX})

# Define interface library target
add_library(libcommute INTERFACE)
target_include_directories(libcommute INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Build unit tests
if(Tests)
  enable_testing()
  message(STATUS "Building unit tests")
  add_subdirectory(test)
endif(Tests)

# Build examples
if(Examples)
  message(STATUS "Building examples")
  add_subdirectory(example)
endif(Examples)
