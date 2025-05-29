cmake_minimum_required(VERSION 3.15)

project(ComputeEnergy
    LANGUAGES CXX
    DESCRIPTION "ComputeEnergy: A C++ library for computing energy values of pseudoknotten structures"
)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(ComputEnergySources
    src/ComputeEnergy.cpp
    src/RNARegions/Bands.cpp
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
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE ComputeEnergyLib)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
