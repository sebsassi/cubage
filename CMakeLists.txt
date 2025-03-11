# Copyright (c) 2024 Sebastian Sassi

# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal # in the Software without restriction, including without limitation the rights # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.
cmake_minimum_required(VERSION 3.12)
project(cubage VERSION 0.1.0 LANGUAGES CXX)

add_library(cubage INTERFACE)

target_include_directories(cubage
    INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION include/cubage)

install(TARGETS cubage EXPORT cubage-targets)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "cmake/config.cmake.in"
    "cubage/cubage-config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cubage
    PATH_VARS CMAKE_INSTALL_LIBDIR
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cubage/cubage-config-version.cmake"
    VERSION ${CMAKE_PROJECT_VERSION}
    COMPATIBILITY ExactVersion
)

install(EXPORT cubage-targets
    FILE cubage-targets.cmake
    NAMESPACE cubage::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cubage
)

install(
    FILES 
        "${CMAKE_CURRENT_BINARY_DIR}/cubage/cubage-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/cubage/cubage-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cubage
    COMPONENT Devel
)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()