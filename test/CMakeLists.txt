cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GoogleTest)

set(TEST_TARGET KnotergyTests)

file(GLOB SOURCES
    tests/*.cpp
)

add_executable(${TEST_TARGET} ${SOURCES})

target_link_libraries(${TEST_TARGET} PUBLIC
    gtest_main
    KnotergyLib
)
target_include_directories(${TEST_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Conditionally apply sanitizers if in Debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug build detected in test target, enabling sanitizers")
    target_compile_options(${TEST_TARGET} PRIVATE        
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
        -fno-omit-frame-pointer
    )
    target_link_options(${TEST_TARGET} PRIVATE         
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
    )
endif()

gtest_discover_tests(${TEST_TARGET}
    DISCOVERY_MODE PRE_TEST
    EXTRA_ARGS --gtest_color=yes
)
