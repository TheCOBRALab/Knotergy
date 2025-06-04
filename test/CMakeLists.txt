cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GoogleTest)

set(TEST_TARGET ComputeEnergyTests)

set(SOURCES
    tests/test_RNAEntry.cpp
)

add_executable(${TEST_TARGET} ${SOURCES})
target_link_libraries(${TEST_TARGET} PUBLIC
    gtest_main
    ComputeEnergyLib
)
target_include_directories(${TEST_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/src)

gtest_discover_tests(${TEST_TARGET}
    DISCOVERY_MODE PRE_TEST
    EXTRA_ARGS --gtest_color=yes
)
