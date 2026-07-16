cmake_minimum_required(VERSION 3.20)


project(
    dynamont
    VERSION 0.0.0
    LANGUAGES CXX
)


# ------------------------------------------------------------
# Options
# ------------------------------------------------------------

option(BUILD_TESTING "Build tests" OFF)


# ------------------------------------------------------------
# Compiler configuration
# ------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
        -Wno-unused-parameter
    )
endif()


# ------------------------------------------------------------
# Version
# ------------------------------------------------------------

if(PROJECT_VERSION)
    set(DYNAMONT_VERSION "${PROJECT_VERSION}")
else()
    set(DYNAMONT_VERSION "unknown")
endif()


configure_file(
    ${PROJECT_SOURCE_DIR}/include/version.hpp.in
    ${PROJECT_BINARY_DIR}/version.hpp
    @ONLY
)



# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------

find_package(OpenMP REQUIRED)



# ------------------------------------------------------------
# Common include paths
# ------------------------------------------------------------

set(PROJECT_INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/include
)


# ------------------------------------------------------------
# Utils
# ------------------------------------------------------------

add_library(utils
    src/cpp/utils.cpp
)

target_link_libraries(
    utils
    PUBLIC
        OpenMP::OpenMP_CXX
)

target_include_directories(
    utils
    PUBLIC
        ${PROJECT_INCLUDE_DIR}
)


# ------------------------------------------------------------
# NT
# ------------------------------------------------------------

add_library(NT_library
    src/cpp/NT.cpp
)


target_link_libraries(
    NT_library
    PUBLIC
        utils
        OpenMP::OpenMP_CXX
)


target_include_directories(
    NT_library
    PUBLIC
        ${PROJECT_INCLUDE_DIR}
        ${PROJECT_BINARY_DIR}
)


add_executable(
    dynamont-NT
    src/cpp/NT_main.cpp
)


target_link_libraries(
    dynamont-NT
    PRIVATE
        NT_library
)



# ------------------------------------------------------------
# NTC
# ------------------------------------------------------------

add_library(
    NTC_library
    src/cpp/NTC.cpp
)


target_link_libraries(
    NTC_library
    PUBLIC
        utils
        OpenMP::OpenMP_CXX
)


target_include_directories(
    NTC_library
    PUBLIC
        ${PROJECT_INCLUDE_DIR}
        ${PROJECT_BINARY_DIR}
)



add_executable(
    dynamont-NTC
    src/cpp/NTC_main.cpp
)


target_link_libraries(
    dynamont-NTC
    PRIVATE
        NTC_library
)



# ------------------------------------------------------------
# NT banded
# ------------------------------------------------------------

add_library(
    NT_banded_library
    src/cpp/NT_banded.cpp
)


target_link_libraries(
    NT_banded_library
    PUBLIC
        utils
        OpenMP::OpenMP_CXX
)


target_include_directories(
    NT_banded_library
    PUBLIC
        ${PROJECT_INCLUDE_DIR}
        ${PROJECT_BINARY_DIR}
)



add_executable(
    dynamont-NT-banded
    src/cpp/NT_banded_main.cpp
)


target_link_libraries(
    dynamont-NT-banded
    PRIVATE
        NT_banded_library
)



# ------------------------------------------------------------
# Tests
# ------------------------------------------------------------

if(BUILD_TESTING)

    include(FetchContent)

    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
    )

    FetchContent_MakeAvailable(googletest)


    enable_testing()


    add_executable(
        test_dynamont
        test/test.cpp
    )


    target_compile_definitions(
        test_dynamont
        PRIVATE
            UNIT_TESTING
    )


    target_link_libraries(
        test_dynamont
        PRIVATE
            NTC_library
            GTest::gtest
            GTest::gtest_main
    )

endif()



# ------------------------------------------------------------
# Installation
# ------------------------------------------------------------

install(
    TARGETS
        dynamont-NTC
        dynamont-NT-banded
        dynamont-NT

    DESTINATION bin
)