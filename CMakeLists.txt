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