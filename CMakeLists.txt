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
        $<$<CONFIG:Release>:-march=native>
        $<$<CONFIG:Release>:-mtune=native>
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

if(BUILD_PYTHON_BINDING)
    find_package(pybind11 CONFIG REQUIRED)
endif()

# ------------------------------------------------------------
# Common include paths
# ------------------------------------------------------------

set(PROJECT_INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/include
)


add_library(
    aligner_api_library
    src/cpp/aligner.cpp
    src/cpp/NT_aligner_api.cpp
    src/cpp/NTK_aligner_api.cpp
)

set_target_properties(
    aligner_api_library
    PROPERTIES
        POSITION_INDEPENDENT_CODE ON
)

target_link_libraries(
    aligner_api_library
    PUBLIC
        OpenMP::OpenMP_CXX
)

target_include_directories(
    aligner_api_library
    PUBLIC
        ${PROJECT_INCLUDE_DIR}
        ${PROJECT_BINARY_DIR}
)

if(BUILD_PYTHON_BINDING)
    pybind11_add_module(
        _dynamont
        src/cpp/aligner_bindings.cpp
    )

    target_link_libraries(
        _dynamont
        PRIVATE
            aligner_api_library
    )

    install(
        TARGETS
            _dynamont
        LIBRARY DESTINATION python
        RUNTIME DESTINATION python
    )
endif()