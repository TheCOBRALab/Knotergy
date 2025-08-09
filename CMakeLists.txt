cmake_minimum_required(VERSION 3.15)

project(Knotergy
    LANGUAGES CXX
    DESCRIPTION "Knotergy: A C++ library for computing energy values of pseudoknotten structures"
)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------- Configurations -------------------------

# Default build type if none is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# -------------------------------------------------------------------
#                             C++ Flags
# -------------------------------------------------------------------

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
        -Wsign-compare
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


# -------------------------------------------------------------------
#                      GoogleTest Configuration
# -------------------------------------------------------------------
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




# -------------------------------------------------------------------
#                  ViennaRNA Library Configuration
# -------------------------------------------------------------------

# Get Raw user input for library path
set(_raw_user_input_lib "${VIENNARNA_LIB}") # user input using -D flag
if (NOT _raw_user_input_lib)
    set(_raw_user_input_lib "$ENV{VIENNARNA_LIB}") # user input using environment variable
endif()

# Set path to the cache variable VIENNARNA_LIB
if(_raw_user_input_lib)
    set(VIENNARNA_LIB
        "${_raw_user_input_lib}"
        CACHE PATH "Path to ViennaRNA headers" FORCE)
endif()

# Locate ViennaRNA library
find_library(VIENNARNA_LIBRARY
    NAMES vrna RNA
    HINTS ${VIENNARNA_LIB}
    PATHS
        /usr/local/lib/
        /usr/lib
        /opt/homebrew/lib
        /opt/vienna/lib
        /usr/local/lib64
        /usr/lib64
    DOC "Path to ViennaRNA library"
)

# ---------------- Clean include DIR path ----------------

# Get user input for include directory
set(_raw_user_input_dir "${VIENNARNA_INCLUDE_DIR}") # user input using -D flag
if (NOT _raw_user_input_dir)
    set(_raw_user_input_dir "$ENV{VIENNARNA_INCLUDE_DIR}") # user input using environment variable
endif()

set(_clean_user_input_dir "${_raw_user_input_dir}")

if(_clean_user_input_dir MATCHES "(.*/ViennaRNA)/pair_mat.h$")
    get_filename_component(_clean_user_input_dir "${_clean_user_input_dir}" DIRECTORY)  # strip pair_mat.h
    get_filename_component(_clean_user_input_dir "${_clean_user_input_dir}" DIRECTORY)  # strip ViennaRNA
elseif(_clean_user_input_dir MATCHES ".*/ViennaRNA/?$")
    get_filename_component(_clean_user_input_dir "${_clean_user_input_dir}" DIRECTORY) # strip ViennaRNA
endif()

# Set the cleaned path to the cache variable VIENNARNA_INCLUDE_DIR
if(_clean_user_input_dir)
    set(VIENNARNA_INCLUDE_DIR
        "${_clean_user_input_dir}"
        CACHE PATH "Path to ViennaRNA headers" FORCE)
endif()


# ------------------ Find the include directory -----------------

# Search for the pair_mat.h file
find_path(VIENNARNA_INCLUDE_DIRECTORY
          NAMES  ViennaRNA/pair_mat.h
          HINTS "${VIENNARNA_INCLUDE_DIR}" # first try the cleaned path
          PATHS
                /usr/local/include
                /usr/include
                /opt/homebrew/include
                /opt/vienna/include
          DOC   "Path to ViennaRNA headers")

message(STATUS "VIENNARNA_LIBRARY = ${VIENNARNA_LIBRARY}")
message(STATUS "VIENNARNA_INCLUDE_DIRECTORY = ${VIENNARNA_INCLUDE_DIRECTORY}")

# ------------------ Check if the library and include directory were found -----------------

if(NOT VIENNARNA_LIBRARY OR NOT VIENNARNA_INCLUDE_DIRECTORY)
message(FATAL_ERROR "
❌ ViennaRNA library not found!

To fix this, either:

1. Set environment variables:
   export VIENNARNA_LIB=/path/to/libRNA.a
   export VIENNARNA_INCLUDE_DIR=/path/to/include

2. Or pass them to CMake directly:
   cmake -S . -B build \\
     -DVIENNARNA_LIB=/path/to/libRNA.a \\
     -DVIENNARNA_INCLUDE_DIR=/path/to/include

Tip: To find the files, try:
   find /usr/local /home -name libRNA.a
   find /usr/local /home -name pair_mat.h
")

endif()




# -------------------------------------------------------------------
#                     Knotergy Configuration
# -------------------------------------------------------------------

# ------------------------- Source Files -------------------------
file(GLOB KnotergySources
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pipeline/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/preprocessing/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/loop_tree/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/energy/*.cpp
)

add_library(KnotergyLib ${KnotergySources})

target_include_directories(KnotergyLib SYSTEM PUBLIC ${VIENNARNA_INCLUDE_DIRECTORY})
target_link_libraries(KnotergyLib PUBLIC ${VIENNARNA_LIBRARY})
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE KnotergyLib)

# ------------------------- Install -------------------------
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
