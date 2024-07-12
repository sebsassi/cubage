cmake_minimum_required(VERSION 3.12)
project(cubage VERSION 0.1.0 LANGUAGES CXX)

add_library(cubage INTERFACE)

target_include_directories(cubage
    INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)