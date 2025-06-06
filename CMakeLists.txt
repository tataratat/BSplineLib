# Copyright (c) 2018–2021 SplineLib
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# Project
project(
  BSplineLib
  VERSION 1.0.0
  LANGUAGES CXX)
set(namespace "${PROJECT_NAME}::")
set(targets_export_name "${PROJECT_NAME}Targets")
option(BSPLINELIB_SHARED "Build shared library" OFF)
option(BSPLINELIB_BUILD_TOOLS "Build tools" OFF)

# Overwrite some options if this is for splinepy
if(SPLINEPY_BUILD_BSPLINELIB)
  set(BSPLINELIB_SHARED OFF)
  set(BSPLINELIB_BUILD_TOOLS OFF)
endif()

# Setup
set(include_install_directory "include")
set(include_install_directory_headers
    "${include_install_directory}/${PROJECT_NAME}")
set(INCLUDE_DIRECTORIES
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Sources>
    $<INSTALL_INTERFACE:${include_install_directory_headers}>)
set(library_install_directory "lib")
set(executable_install_directory "bin")
set(configuration_install_directory "lib/cmake/${PROJECT_NAME}")

# Compilation
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(OPTIMIZATION_FLAGS $<$<NOT:$<CONFIG:Debug>>:-O3>)
  set(PARALLELIZATION_FLAGS)
  set(RUNTIME_CHECKS_DEBUG)
  set(WARNING_FLAGS
      -Wall
      -Wextra
      -Wmost
      -Wextra
      -Wpedantic
      -Wunreachable-code
      -Wshadow
      -Wfloat-equal
      -Weffc++
      $<$<NOT:$<CONFIG:Debug>>:-Wno-unused-parameter>
      $<$<NOT:$<CONFIG:Debug>>:-Wno-unused-variable>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(OPTIMIZATION_FLAGS $<$<NOT:$<CONFIG:Debug>>:-O3>)
  set(PARALLELIZATION_FLAGS)
  set(RUNTIME_CHECKS_DEBUG)
  set(WARNING_FLAGS -Wall -Wextra -Wpedantic -Wzero-as-null-pointer-constant
                    $<$<NOT:$<CONFIG:Debug>>:-Wno-unused>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  set(OPTIMIZATION_FLAGS -fast)
  set(PARALLELIZATION_FLAGS -vec)
  set(RUNTIME_CHECKS_DEBUG -ftrapv)
  set(WARNING_FLAGS -Wall -Wextra)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(OPTIMIZATION_FLAGS -fast)
  set(PARALLELIZATION_FLAGS -vec -parallel)
  set(RUNTIME_CHECKS_DEBUG -fp-trap=all -fp-trap-all=all -ftrapuv)
  set(WARNING_FLAGS -Wall -Wextra -Wcheck -Weffc++)
else()
  message(WARNING "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()
set(COMPILE_DEFINITIONS $<$<BOOL:${BSPLINELIB_SHARED}>:BSPLINELIB_SHARED>)
set(COMPILE_OPTIONS
    ${OPTIMIZATION_FLAGS} $<IF:$<CONFIG:Release>,${PARALLELIZATION_FLAGS},
    ${RUNTIME_CHECKS_DEBUG}> ${WARNING_FLAGS} ${SPLINEPY_FLAGS})
message(STATUS "Compile using ${CMAKE_CXX_COMPILER_ID} compiler!")
set(CMAKE_DEBUG_POSTFIX "d")

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "Build debug!")
  set(BSPLINELIB_COMPILE_FEATURES cxx_std_20)
else()
  message(STATUS "Build release!")
  # MSVC has issues evaluating constexpr @ string_operation.inc l.125 with c++17
  if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(BSPLINELIB_COMPILE_FEATURES cxx_std_20)
  else()
    set(BSPLINELIB_COMPILE_FEATURES cxx_std_17)
  endif()
endif()
if(BSPLINELIB_SHARED)
  message(STATUS "Build shared!")
  set(BUILD_SHARED_LIBS ON)
else()
  message(STATUS "Build static!")
  set(BUILD_SHARED_LIBS OFF)
endif()

# Content
add_subdirectory(Sources/BSplineLib)

if(BSPLINELIB_BUILD_TOOLS)
  add_subdirectory(Tools)
endif()

# Packaging
set(generated_directory ${CMAKE_CURRENT_BINARY_DIR}/Generated)
set(version_configuration
    ${generated_directory}/${PROJECT_NAME}ConfigVersion.cmake)
set(project_configuration ${generated_directory}/${PROJECT_NAME}Config.cmake)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  CMake/config.cmake.in ${project_configuration}
  INSTALL_DESTINATION ${configuration_install_directory})
write_basic_package_version_file(${version_configuration}
                                 COMPATIBILITY SameMajorVersion)

install(FILES ${project_configuration} ${version_configuration}
        DESTINATION ${configuration_install_directory})
install(
  EXPORT ${targets_export_name}
  NAMESPACE ${namespace}
  DESTINATION ${configuration_install_directory})
