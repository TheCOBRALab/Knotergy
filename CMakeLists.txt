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
)

add_library(ComputeEnergyLib ${ComputEnergySources})
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE ComputeEnergyLib)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
