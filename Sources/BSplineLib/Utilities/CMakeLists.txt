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

set(HEADERS
    error_handling.inl
    error_handling.hpp
    index.inl
    index.hpp
    math_operations.hpp
    named_type.inl
    named_type.hpp
    numeric_operations.inl
    numeric_operations.hpp
    containers.inl
    containers.hpp
    string_operations.inl
    string_operations.hpp
    system_operations.hpp
    system_operations.inl)

set(SOURCES
    math_operations.cpp string_operations.cpp system_operations.cpp
    #
    ${HEADERS})
set_source_files_properties(${HEADERS} PROPERTIES LANGUAGE CXX HEADER_FILE_ONLY
                                                               TRUE)

add_library(utilities ${SOURCES})
add_library(BSplineLib::utilities ALIAS utilities)

target_include_directories(utilities PUBLIC ${INCLUDE_DIRECTORIES})
target_compile_definitions(utilities PUBLIC ${COMPILE_DEFINITIONS})
target_compile_options(utilities PRIVATE ${COMPILE_OPTIONS})
target_compile_features(utilities PUBLIC ${BSPLINELIB_COMPILE_FEATURES})

install(
  TARGETS utilities
  EXPORT ${targets_export_name}
  LIBRARY DESTINATION ${library_install_directory}
  ARCHIVE DESTINATION ${library_install_directory}/static
  RUNTIME DESTINATION ${executable_install_directory}
  INCLUDES
  DESTINATION ${include_install_directory})

install(FILES ${HEADERS}
        DESTINATION ${include_install_directory_headers}/Utilities)
