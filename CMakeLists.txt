cmake_minimum_required(VERSION 3.14)
project(dynamont)

# GoogleTest requires at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set RPATH to include the Conda environment lib directory
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# fetch GTest version 1.15.2 from github
include(FetchContent)
FetchContent_Declare(
  GTest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
)
FetchContent_MakeAvailable(GTest)

# Disable installation of GTest system-wide
set(INSTALL_GTEST OFF CACHE BOOL "Disable GTest installation" FORCE)

# Add the include directory to the search path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add library files
add_library(dynamont_lib src/cpp/utils.cpp)
# target_include_directories(dynamont_lib PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Create the first executable dynamont-NT
add_executable(dynamont-NT src/cpp/dynamont_NT.cpp)
# Link dynamont_lib to dynamont-NT
target_link_libraries(dynamont-NT PRIVATE dynamont_lib)
# Include header files for dynamont_lib
target_include_directories(dynamont-NT PRIVATE include)

# Create the second executable dynamont-NTC
add_executable(dynamont-NTC src/cpp/dynamont_NTC.cpp)
# Link dynamont_lib to dynamont-NTC
target_link_libraries(dynamont-NTC PRIVATE dynamont_lib)
# Include header files for dynamont_lib
target_include_directories(dynamont-NTC PRIVATE include)


# Create the test executable
# need global variables from dynamont-NTC.cpp
add_executable(test_dynamont test/test.cpp src/cpp/dynamont_NTC.cpp)
# exclude main from dynamont-NTC.cpp with UNIT_TESTING definition
target_compile_definitions(test_dynamont PRIVATE UNIT_TESTING)
# add GTest
find_package(GTest REQUIRED)
# Link GTest
target_link_libraries(test_dynamont PRIVATE GTest::GTest GTest::Main)
# Link dynamont_lib to the test executable
target_link_libraries(test_dynamont PRIVATE dynamont_lib)
# Include header files for the test target
target_include_directories(test_dynamont PRIVATE include)
# Add the test to CTest
add_test(NAME DynamontTests COMMAND test_dynamont)

# Specify installation location
install(TARGETS dynamont-NT DESTINATION bin)
install(TARGETS dynamont-NTC DESTINATION bin)
install(TARGETS test_dynamont DESTINATION bin)