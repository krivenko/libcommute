#
# This file is part of libcommute, a quantum operator algebra DSL and
# exact diagonalization toolkit for C++11/14/17.
#
# Copyright (C) 2016-2025 Igor Krivenko
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

macro(FindPythonModule Name ImportName)
  execute_process(COMMAND ${PYTHON_EXECUTABLE}
                  "-c" "import ${ImportName}; print(${ImportName}.__version__)"
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  RESULT_VARIABLE ${Name}_FOUND
                  OUTPUT_VARIABLE ${Name}_VERSION)

  if(${${Name}_FOUND} EQUAL 0)
    set(${Name}_FOUND TRUE)
    string(REPLACE "." ";" ${Name}_VERSION_LIST ${${Name}_VERSION})
    list(GET ${Name}_VERSION_LIST 0 ${Name}_VERSION_MAJOR)
    list(GET ${Name}_VERSION_LIST 1 ${Name}_VERSION_MINOR)
    list(GET ${Name}_VERSION_LIST 2 ${Name}_VERSION_PATCH)
    unset(${Name}_VERSION_LIST)
  else(${${Name}_FOUND} EQUAL 0)
    unset(${Name}_FOUND)
    unset(${Name}_VERSION)
  endif(${${Name}_FOUND} EQUAL 0)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(${Name}
    FOUND_VAR ${Name}_FOUND
    VERSION_VAR ${Name}_VERSION
    REQUIRED_VARS ${Name}_FOUND
    FAIL_MESSAGE "Failed to find ${Name} Python module"
  )
endmacro(FindPythonModule Name ImportName)
