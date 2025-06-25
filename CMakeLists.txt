cmake_minimum_required(VERSION 3.14)
project(dynamont)

# Enable OpenMP
find_package(OpenMP REQUIRED)

# GoogleTest requires at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Link libstdc++ and libgcc statically, but avoid fully static linking of glibc
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")

# Ensure compatibility with older glibc versions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

# Set RPATH to include the Conda environment lib directory
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# fetch GTest version 1.15.2 from github
include(FetchContent)
FetchContent_Declare(
  GTest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE  # Add this line
)
FetchContent_MakeAvailable(GTest)

# Disable installation of GTest system-wide
set(INSTALL_GTEST OFF CACHE BOOL "Disable GTest installation" FORCE)

# Add the include directory to the search path
include_directories(${PROJECT_SOURCE_DIR}/include)

# generate version.hpp
# CMakeLists.txt
if(DEFINED ENV{PROJECT_VERSION})
    set(PROJECT_VERSION $ENV{PROJECT_VERSION})
endif()

# Write the version to a header file
configure_file(
    ${CMAKE_SOURCE_DIR}/include/version.hpp.in
    ${CMAKE_BINARY_DIR}/version.hpp
    @ONLY
)

include_directories(${CMAKE_BINARY_DIR})


# Create the first executable dynamont-NT
add_library(NT_library src/cpp/utils.cpp src/cpp/NT.cpp)
add_executable(dynamont-NT src/cpp/NT_main.cpp)
target_link_libraries(dynamont-NT PRIVATE NT_library)
if(OpenMP_CXX_FOUND)
  target_link_libraries(dynamont-NT PUBLIC OpenMP::OpenMP_CXX)
endif()
# Include header files for NT_library
target_include_directories(dynamont-NT PRIVATE include)

# Create the second executable dynamont-NTC
add_library(NTC_library src/cpp/utils.cpp src/cpp/NTC.cpp)
add_executable(dynamont-NTC src/cpp/NTC_main.cpp)
# Link NTC_library to dynamont-NTC
target_link_libraries(dynamont-NTC PRIVATE NTC_library)
if(OpenMP_CXX_FOUND)
  target_link_libraries(dynamont-NTC PUBLIC OpenMP::OpenMP_CXX)
endif()
# Include header files for NTC_library
target_include_directories(dynamont-NTC PRIVATE include)

# Create the second executable dynamont-NT-banded
add_library(NT_banded_library src/cpp/utils.cpp src/cpp/NT_banded.cpp)
add_executable(dynamont-NT-banded src/cpp/NT_banded_main.cpp)
# Link NT_banded_library to dynamont-NTbanded
target_link_libraries(dynamont-NT-banded PRIVATE NT_banded_library)
if(OpenMP_CXX_FOUND)
  target_link_libraries(dynamont-NT-banded PUBLIC OpenMP::OpenMP_CXX)
endif()
# Include header files for NT_banded_library
target_include_directories(dynamont-NT-banded PRIVATE include)

# Create the test executable
# need global variables from dynamont-NTC.cpp
add_executable(test_dynamont test/test.cpp src/cpp/NTC.cpp)
# exclude main from dynamont-NTC.cpp with UNIT_TESTING definition
target_compile_definitions(test_dynamont PRIVATE UNIT_TESTING)
# add GTest
find_package(GTest REQUIRED)
# Link GTest
target_link_libraries(test_dynamont PRIVATE GTest::GTest GTest::Main)
# Link dynamont_lib to the test executable
target_link_libraries(test_dynamont PRIVATE NTC_library)
if(OpenMP_CXX_FOUND)
  target_link_libraries(test_dynamont PUBLIC OpenMP::OpenMP_CXX)
endif()
# Include header files for the test target
target_include_directories(test_dynamont PRIVATE include)
# Add the test to CTest
add_test(NAME DynamontTests COMMAND test_dynamont)

# Specify installation location
install(TARGETS dynamont-NTC DESTINATION bin)
install(TARGETS dynamont-NT-banded DESTINATION bin)
install(TARGETS dynamont-NT DESTINATION bin)
install(TARGETS test_dynamont DESTINATION bin)