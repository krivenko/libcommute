#
# This file is part of libcommute, a quantum operator algebra DSL and
# exact diagonalization toolkit for C++11/14/17.
#
# Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Start configuration
cmake_minimum_required(VERSION 3.8.0 FATAL_ERROR)

# Project name and version
set(LIBCOMMUTE_VERSION 0.3)
project(libcommute VERSION ${LIBCOMMUTE_VERSION} LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

message(STATUS "Configuring ${PROJECT_NAME} version ${PROJECT_VERSION}")

# Prohibit in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. "
                      "Please make a separate build directory.")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

# CMake options
option(TESTS "Build unit tests" ON)
option(EXAMPLES "Build examples" ON)

# Install C++ headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include
        DESTINATION ${CMAKE_INSTALL_PREFIX})

# Define interface library target
add_library(libcommute INTERFACE)
target_include_directories(libcommute INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Write config version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/libcommuteConfigVersion.cmake"
  VERSION ${LIBCOMMUTE_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Install CMake files
install(TARGETS libcommute
  EXPORT libcommuteTargets
  PUBLIC_HEADER DESTINATION include COMPONENT Development
)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/libcommuteConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/libcommuteConfig.cmake"
  INSTALL_DESTINATION lib/cmake
)
install(EXPORT libcommuteTargets DESTINATION lib/cmake)
install(FILES "${PROJECT_BINARY_DIR}/libcommuteConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/libcommuteConfig.cmake"
        DESTINATION lib/cmake)

# Build unit tests
if(TESTS)
  enable_testing()
  message(STATUS "Building unit tests")
  add_subdirectory(test)
endif(TESTS)

# Build examples
if(EXAMPLES)
  message(STATUS "Building examples")
  add_subdirectory(example)
endif(EXAMPLES)
