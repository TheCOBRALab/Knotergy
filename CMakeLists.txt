cmake_minimum_required(VERSION 3.15)

project(ComputeEnergy
    LANGUAGES CXX
    DESCRIPTION "ComputeEnergy: A C++ library for computing energy values of pseudoknotten structures"
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


# ------------------------- build -------------------------
set(ComputEnergySources
    src/helpers/main_helpers.cpp
    src/rna_regions/Bands.cpp
    src/loop/Loop.cpp
)

# Enable debug symbols and extra warnings
set(CMAKE_BUILD_TYPE Debug)

# Enable maximum warnings and debugging info for GCC, Clang, and MSVC
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

add_link_options(
    -fsanitize=address
)

add_library(ComputeEnergyLib ${ComputEnergySources})
target_include_directories(ComputeEnergyLib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/rna_regions
        ${CMAKE_CURRENT_SOURCE_DIR}/src/helpers
        ${CMAKE_CURRENT_SOURCE_DIR}/src/loop
)

add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE ComputeEnergyLib)


install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
