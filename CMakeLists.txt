cmake_minimum_required(VERSION 3.15)

project(Knotergy
    VERSION 0.3.0
    LANGUAGES CXX
    DESCRIPTION "Knotergy: A C++ library for computing energy values of pseudoknotted structures"
)


# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Define a global compile definition for the source directory,
# so it can be used in both the library and tests
add_compile_definitions(
    KNOTERGY_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

# ------------------------- External Dependencies ------------------------- 
include(FetchContent)

message(STATUS "Fetching external dependencies... this may take a while on the first run")

# ----------------- nlohmann/json -----------------
set(JSON_Install OFF CACHE BOOL "" FORCE)

FetchContent_Declare(nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.12.0
)
FetchContent_MakeAvailable(nlohmann_json)

# ----------------- uni-algo -----------------
set(UNI_ALGO_HEADER_ONLY ON CACHE BOOL "" FORCE)

FetchContent_Declare(uni_algo
  GIT_REPOSITORY https://github.com/uni-algo/uni-algo.git
  GIT_TAG        v1.2.0
)
FetchContent_MakeAvailable(uni_algo)


#------------------ ThermoRNA (ViennaRNA) -----------------
set(THERMORNA_KNOTERGY OFF CACHE BOOL "" FORCE) # Knotergy implements this already
set(THERMORNA_VIENNARNA ON CACHE BOOL "" FORCE)
# set(FETCHCONTENT_SOURCE_DIR_THERMORNA
#     "path/to/ThermoRNA" # Replace with the actual path to ThermoRNA for local development
# )
FetchContent_Declare(ThermoRNA
    GIT_REPOSITORY https://github.com/TheCOBRALab/ThermoRNA.git
    GIT_TAG v0.1.2
)
FetchContent_MakeAvailable(ThermoRNA)
# ------------------------- Configurations -------------------------

# Default build type if none is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# -------------------------------------------------------------------
#                             C++ Flags
# -------------------------------------------------------------------

# Debug configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "Strict")
    message(STATUS "Building in Debug mode")
    add_compile_options(
        -g
        -O0
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wsign-compare
        -Wconversion
        -Wsign-conversion
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
        -fno-omit-frame-pointer
    )

    if(CMAKE_BUILD_TYPE STREQUAL "Strict")
        add_compile_options(-Werror)
    endif()

    add_link_options(
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
    )
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


# -------------------------------------------------------------------
#                      GoogleTest Configuration
# -------------------------------------------------------------------
option(BUILD_TESTING "Build unit tests" ON)

if (NOT BUILD_TESTING)
    message(STATUS "Unit tests are disabled. To enable them, set -DBUILD_TESTING=ON")
elseif(EXISTS "${CMAKE_SOURCE_DIR}/test/googletest/CMakeLists.txt")
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

# -------------------------------------------------------------------
#                     Knotergy Configuration
# -------------------------------------------------------------------

# ------------------------- Source Files -------------------------
file(GLOB_RECURSE KnotergySources
    ${CMAKE_SOURCE_DIR}/src/io/*.cpp
    ${CMAKE_SOURCE_DIR}/src/preprocessing/*.cpp
    ${CMAKE_SOURCE_DIR}/src/loop_tree/*.cpp
    ${CMAKE_SOURCE_DIR}/src/energy/*.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/*.cpp
)

add_library(KnotergyLib ${KnotergySources})

target_include_directories(KnotergyLib PUBLIC ${CMAKE_SOURCE_DIR}/src)

# ------------------------- Link Libraries -------------------------
target_link_libraries(KnotergyLib
  PUBLIC
    nlohmann_json::nlohmann_json
    uni-algo::uni-algo
    ThermoRNA::ThermoRNA
)

# ------------------------- Executable -------------------------
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE KnotergyLib)


# ------------------------ Platform-Specific Link Options -------------------------
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_link_options(${PROJECT_NAME} PRIVATE
            "$<$<CONFIG:Release>:-static-libstdc++>"
            "$<$<CONFIG:Release>:-static-libgcc>"
        )
    endif()

elseif(APPLE)
    target_link_options(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Release>:-Wl,-dead_strip>"
    )
endif()


# ------------------------- Install -------------------------
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/params DESTINATION share/knotergy)

# -------------------------------------------------------------------
#                         clang-format
# -------------------------------------------------------------------

find_program(CLANG_FORMAT_EXE NAMES clang-format)

file(GLOB_RECURSE KNOTERGY_FORMAT_FILES
    CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_SOURCE_DIR}/test/*.cpp"
    "${CMAKE_SOURCE_DIR}/test/*.hpp"
)

if(CLANG_FORMAT_EXE)
    add_custom_target(format
        COMMAND ${CLANG_FORMAT_EXE} -i ${KNOTERGY_FORMAT_FILES}
        COMMENT "Formatting C++ source files with clang-format"
        VERBATIM
    )

    add_custom_target(format-check
        COMMAND ${CLANG_FORMAT_EXE} --dry-run --Werror ${KNOTERGY_FORMAT_FILES}
        COMMENT "Checking C++ source formatting with clang-format"
        VERBATIM
    )
else()
    add_custom_target(format
        COMMAND ${CMAKE_COMMAND} -E echo "clang-format not found"
        COMMAND ${CMAKE_COMMAND} -E false
        VERBATIM
    )

    add_custom_target(format-check
        COMMAND ${CMAKE_COMMAND} -E echo "clang-format not found"
        COMMAND ${CMAKE_COMMAND} -E false
        VERBATIM
    )
endif()