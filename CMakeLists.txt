cmake_minimum_required(VERSION 3.20)


project(
    ${SKBUILD_PROJECT_NAME}
    VERSION ${SKBUILD_PROJECT_VERSION}
    LANGUAGES CXX
)


# ------------------------------------------------------------
# Options
# ------------------------------------------------------------

option(BUILD_TESTING "Build tests" OFF)
option(BUILD_PYTHON_BINDING "Build the pybind11 Python extension" OFF)


# ------------------------------------------------------------
# Compiler configuration
# ------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)

if(IPO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()


if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
        -Wno-unused-parameter
        $<$<CONFIG:Release>:-O3>
    )
endif()


# ------------------------------------------------------------
# Version
# ------------------------------------------------------------

configure_file(
    ${PROJECT_SOURCE_DIR}/include/dynamont/version.hpp.in
    ${PROJECT_BINARY_DIR}/generated/dynamont/version.hpp
    @ONLY
)
# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------

# find_package(OpenMP REQUIRED)

if(BUILD_PYTHON_BINDING)
    find_package(pybind11 CONFIG REQUIRED)
endif()

# ------------------------------------------------------------
# Include directories
# ------------------------------------------------------------

set(DYNAMONT_INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/include
)

set(DYNAMONT_GENERATED_INCLUDE_DIR
    ${PROJECT_BINARY_DIR}/generated
)

# ------------------------------------------------------------
# Core C++ library
# ------------------------------------------------------------

add_library(
    dynamont_core
    src/cpp/aligner.cpp
    src/cpp/NT_aligner_api.cpp
    src/cpp/NTK_aligner_api.cpp
)

set_target_properties(
    dynamont_core
    PROPERTIES
        POSITION_INDEPENDENT_CODE ON
)

# target_link_libraries(
#     dynamont_core
#     PUBLIC
#         OpenMP::OpenMP_CXX
# )

target_include_directories(
    dynamont_core
    PUBLIC
        ${DYNAMONT_INCLUDE_DIR}
        ${DYNAMONT_GENERATED_INCLUDE_DIR}
)

# ------------------------------------------------------------
# Python binding
# ------------------------------------------------------------

if(BUILD_PYTHON_BINDING)

    pybind11_add_module(
        _dynamont
        src/cpp/aligner_bindings.cpp
    )

    target_link_libraries(
        _dynamont
        PRIVATE
            dynamont_core
    )

    target_include_directories(
        _dynamont
        PRIVATE
            ${DYNAMONT_INCLUDE_DIR}
            ${DYNAMONT_GENERATED_INCLUDE_DIR}
    )

    install(
        TARGETS _dynamont
        LIBRARY DESTINATION dynamont
        RUNTIME DESTINATION dynamont
    )

endif()