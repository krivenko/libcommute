#
# This file is part of libcommute, a quantum operator algebra DSL and
# exact diagonalization toolkit for C++11/14/17.
#
# Copyright (C) 2016-2025 Igor Krivenko
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Start configuration
cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)

# Project name and version
set(LIBCOMMUTE_VERSION 0.8.0)
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
option(DOCUMENTATION "Build documentation" OFF)

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
install(EXPORT libcommuteTargets
        NAMESPACE libcommute::
        DESTINATION lib/cmake/libcommute)
install(FILES "${PROJECT_BINARY_DIR}/libcommuteConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/libcommuteConfig.cmake"
        DESTINATION lib/cmake/libcommute)

# Install pkg-config file
configure_file(libcommute.pc.in libcommute.pc @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/libcommute.pc"
        DESTINATION share/pkgconfig
)

# Optionally find Boost
find_package(Boost 1.70 CONFIG)
if(Boost_FOUND)
  message(STATUS "Found Boost version ${Boost_VERSION}")
endif(Boost_FOUND)

# Optionally find Eigen 3
if(NOT POLICY CMP0074)
  set(Eigen3_DIR ${Eigen3_ROOT}/share/eigen3/cmake)
endif(NOT POLICY CMP0074)
find_package(Eigen3 CONFIG)
if(Eigen3_FOUND)
  if(NOT Eigen3_VERSION)
    set(Eigen3_VERSION "${EIGEN3_VERSION_STRING}")
  endif(NOT Eigen3_VERSION)
  message(STATUS "Found Eigen3 version ${Eigen3_VERSION}")

  macro(link_to_eigen3 target)
    if(TARGET Eigen3::Eigen)
      target_link_libraries(${target} PRIVATE Eigen3::Eigen)
    else(TARGET Eigen3::Eigen)
      target_compile_definitions(${target} PRIVATE ${EIGEN3_DEFINITIONS})
      target_include_directories(${target}
                                 SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIRS})
    endif(TARGET Eigen3::Eigen)
    if(NOT ${Eigen3_VERSION} VERSION_LESS "3.3.0")
      target_compile_definitions(${target}
                                 PRIVATE EIGEN_CAN_MIX_REAL_COMPLEX_EXPR)
    endif(NOT ${Eigen3_VERSION} VERSION_LESS "3.3.0")
  endmacro(link_to_eigen3)
endif(Eigen3_FOUND)

# Build documentation
if(DOCUMENTATION)
  message(STATUS "Building documentation")

  # FindPythonInterp has been deprecated since 3.12,
  # but we still require only 3.10.
  cmake_policy(SET CMP0148 OLD)

  # Detect Sphinx, sphinx_rtd_theme and myst-parser
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  find_package(PythonInterp 3 REQUIRED)
  find_package(Sphinx 2.1.0 REQUIRED)
  find_package(SphinxRTDTheme 0.3.0 REQUIRED)
  find_package(MySTParser REQUIRED)
  add_subdirectory(doc)
endif(DOCUMENTATION)

# Enable static analysis
option(STATIC_ANALYSIS "Run static analysis tools on C++ sources" OFF)
if(STATIC_ANALYSIS)
  message(STATUS "Static analysis enabled")
  include(cmake/StaticAnalysis.cmake)
endif(STATIC_ANALYSIS)

# Build unit tests
if(TESTS)
  enable_testing()
  message(STATUS "Building unit tests")
  add_subdirectory(test)
endif(TESTS)

# Build examples
if(EXAMPLES)
  message(STATUS "Building examples")
  add_subdirectory(examples)
endif(EXAMPLES)
