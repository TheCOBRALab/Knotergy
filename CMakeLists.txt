cmake_minimum_required(VERSION 3.15)

project(Knotergy
    LANGUAGES CXX
    DESCRIPTION "Knotergy: A C++ library for computing energy values of pseudoknotten structures"
)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --------------------------------- Google Test ---------------------------------------
option(ENABLE_TESTS "Build unit tests" ON)

if (NOT ENABLE_TESTS)
    message(STATUS "Unit tests are disabled. To enable them, set -DBUILD_TESTING=ON")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/googletest/CMakeLists.txt")
    enable_testing()
    add_subdirectory(test)
    add_subdirectory(test/googletest)
else()
    message(WARNING "
============================================================
 GoogleTest submodule not found. Unit tests will be disabled.
------------------------------------------------------------
To enable unit testing, please do one of the following:

  1. If you already cloned the repository:
     Run the following command to initialize submodules:
         git submodule update --init --recursive

  2. If you haven't cloned the repository yet:
     Clone it with submodules using:
         git clone --recursive <repo-url>

After initializing the submodule, re-run CMake to enable tests.
============================================================")
endif()

# ------------------------- Source Files -------------------------
file(GLOB KnotergySources
    ${CMAKE_CURRENT_SOURCE_DIR}/src/helpers/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/rna_regions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/loops/*.cpp
)

add_library(KnotergyLib ${KnotergySources})

add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE KnotergyLib)

# ------------------------- Configurations -------------------------

# Default build type if none is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Debug configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in Debug mode")
    add_compile_options(
        -g
        -O0
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wconversion
        -Wsign-conversion
        -fsanitize=address
        -fno-omit-frame-pointer
    )
    add_link_options(-fsanitize=address)
else()
    # Normal (Release or other non-Debug builds)
    message(STATUS "Building in Release/Normal mode")
    add_compile_options(
        -O3
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

# ------------------------- Install -------------------------
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
