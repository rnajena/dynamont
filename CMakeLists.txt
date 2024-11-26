cmake_minimum_required(VERSION 3.10)
project(Dynamont)

# Enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the include directory to the search path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add library files
add_library(dynamont_lib src/dynamont/utils.cpp)
# target_include_directories(dynamont_lib PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Create the first executable dynamont_NT
add_executable(dynamont_NT src/dynamont/dynamont_NT.cpp)
# Link dynamont_lib to dynamont_NT
target_link_libraries(dynamont_NT PRIVATE dynamont_lib)
# Include header files for dynamont_lib
target_include_directories(dynamont_NT PRIVATE include)

# Create the second executable dynamont_NTC
add_executable(dynamont_NTC src/dynamont/dynamont_NTC.cpp)
# Link dynamont_lib to dynamont_NTC
target_link_libraries(dynamont_NTC PRIVATE dynamont_lib)
# Include header files for dynamont_lib
target_include_directories(dynamont_NTC PRIVATE include)


# Create the test executable
add_executable(test_dynamont test/test.cpp src/dynamont/dynamont_NTC.cpp src/dynamont/utils.cpp)
# add GTest
find_package(GTest REQUIRED)
target_link_libraries(test_dynamont PRIVATE GTest::GTest GTest::Main)
# Link dynamont_lib to the test executable
target_link_libraries(test_dynamont PRIVATE dynamont_lib)
# Include header files for the test target
target_include_directories(test_dynamont PRIVATE include)
# Add the test to CTest
add_test(NAME DynamontTests COMMAND test_dynamont)