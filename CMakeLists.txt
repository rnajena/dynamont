cmake_minimum_required(VERSION 3.14)
project(Dynamont)

# GoogleTest requires at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# fetch GTest version 1.15.2 from github
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
)
FetchContent_MakeAvailable(googletest)

# Add the include directory to the search path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add library files
add_library(dynamont_lib src/cpp/utils.cpp)
# target_include_directories(dynamont_lib PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Create the first executable dynamont_NT
add_executable(dynamont_NT src/cpp/dynamont_NT.cpp)
# Link dynamont_lib to dynamont_NT
target_link_libraries(dynamont_NT PRIVATE dynamont_lib)
# Include header files for dynamont_lib
target_include_directories(dynamont_NT PRIVATE include)

# Create the second executable dynamont_NTC
add_executable(dynamont_NTC src/cpp/dynamont_NTC.cpp)
# Link dynamont_lib to dynamont_NTC
target_link_libraries(dynamont_NTC PRIVATE dynamont_lib)
# Include header files for dynamont_lib
target_include_directories(dynamont_NTC PRIVATE include)


# Create the test executable
# need global variables from dynamont_NTC.cpp
add_executable(test_dynamont test/test.cpp src/cpp/dynamont_NTC.cpp)
# exclude main from dynamont_NTC.cpp with UNIT_TESTING definition
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