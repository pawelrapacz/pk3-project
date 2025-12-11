cmake_minimum_required(VERSION 3.28.3)
project(crs-matrix VERSION 1.0.0
                   DESCRIPTION "pk3-project"
                   LANGUAGES CXX)

set(PROJECT_AUTHOR "Paweł Rapacz")
string(TIMESTAMP PROJECT_PUBLICATION_YEAR "%Y")

configure_file(config.h.in config.h @ONLY)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

add_subdirectory(Matrix/)

add_executable(${PROJECT_NAME}
    src/main.cpp
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 23
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${PROJECT_BINARY_DIR}/
    src/
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    Matrix
)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic -Werror --pedantic-errors)

target_precompile_headers(${PROJECT_NAME} PRIVATE
    <cstddef>
    <cstdint>
    <memory>
    <functional>
    <utility>
    <cstring>
)
