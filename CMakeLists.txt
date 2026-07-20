cmake_minimum_required(VERSION 3.15)

project(Knotergy
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
#                  ViennaRNA Library Configuration
# This section was drafted with the help of AI. Contributions and improvements are welcome.
# -------------------------------------------------------------------

# Users may pass any of these:
#  -DVIENNARNA_LIBRARY=/path/to/libRNA.a        (file)
#  -DVIENNARNA_LIB=/path/to/lib_or_dir          (file OR directory)
#  -DVIENNARNA_INCLUDE_DIRECTORY=/path/to/include
#  -DVIENNARNA_INCLUDE_DIR=/path/to/include OR /path/to/include/ViennaRNA/pair_mat.h

# ---------- Normalize user inputs ----------
set(_raw_lib "${VIENNARNA_LIBRARY}")
if (NOT _raw_lib)
  set(_raw_lib "${VIENNARNA_LIB}")
endif()
if (NOT _raw_lib)
  set(_raw_lib "$ENV{VIENNARNA_LIBRARY}")
endif()
if (NOT _raw_lib)
  set(_raw_lib "$ENV{VIENNARNA_LIB}")
endif()

# If user pointed to the static file, accept as-is. If it's a directory, search inside.
if (_raw_lib)
  if (EXISTS "${_raw_lib}" AND NOT IS_DIRECTORY "${_raw_lib}")
    set(VIENNARNA_LIBRARY "${_raw_lib}" CACHE FILEPATH "Path to ViennaRNA library" FORCE)
  elseif (IS_DIRECTORY "${_raw_lib}")
    list(APPEND _LIB_HINTS "${_raw_lib}")
  endif()
endif()

# Include dir normalization (user may pass the header file path or the ViennaRNA subdir).
set(_raw_inc "${VIENNARNA_INCLUDE_DIRECTORY}")
if (NOT _raw_inc)
  set(_raw_inc "${VIENNARNA_INCLUDE_DIR}")
endif()
if (NOT _raw_inc)
  set(_raw_inc "$ENV{VIENNARNA_INCLUDE_DIRECTORY}")
endif()
if (NOT _raw_inc)
  set(_raw_inc "$ENV{VIENNARNA_INCLUDE_DIR}")
endif()

set(_inc_hint "")
if (_raw_inc)
  set(_inc_hint "${_raw_inc}")
  # If a file was passed (…/ViennaRNA/pair_mat.h), strip file and ViennaRNA folder.
  if (_inc_hint MATCHES "(.*/ViennaRNA)/pair_mat.h$")
    get_filename_component(_inc_hint "${_inc_hint}" DIRECTORY)  # -> …/ViennaRNA
    get_filename_component(_inc_hint "${_inc_hint}" DIRECTORY)  # -> …/include
  elseif (_inc_hint MATCHES ".*/ViennaRNA/?$")
    get_filename_component(_inc_hint "${_inc_hint}" DIRECTORY)  # -> …/include
  endif()
  if (IS_DIRECTORY "${_inc_hint}")
    set(VIENNARNA_INCLUDE_DIRECTORY "${_inc_hint}" CACHE PATH "Path to ViennaRNA headers" FORCE)
  endif()
endif()

# ---------- Try pkg-config first (if available) ----------
include(CheckIncludeFileCXX)
find_package(PkgConfig QUIET)
if (PkgConfig_FOUND AND (NOT VIENNARNA_LIBRARY OR NOT VIENNARNA_INCLUDE_DIRECTORY))
  # Try several common pc names; the first that resolves will fill VIENNARNA_* vars
  foreach(pcname ViennaRNA RNAlib2 RNAlib vrna rna)
    if (NOT VIENNARNA_LIBRARY OR NOT VIENNARNA_INCLUDE_DIRECTORY)
      pkg_check_modules(PC_VIENNARNA QUIET ${pcname})
      if (PC_VIENNARNA_FOUND)
        # Prefer the first library; fall back to -l flags if needed
        if (PC_VIENNARNA_STATIC_LIBRARIES)
          list(GET PC_VIENNARNA_STATIC_LIBRARIES 0 _pc_lib)
        elseif (PC_VIENNARNA_LIBRARIES)
          list(GET PC_VIENNARNA_LIBRARIES 0 _pc_lib)
        endif()
        if (_pc_lib)
          # If pkg-config gave "-lRNA" style, let find_library resolve it
          if (_pc_lib MATCHES "^-l(.+)$")
            set(_libname "${CMAKE_MATCH_1}")
            find_library(VIENNARNA_LIBRARY NAMES ${_libname} vrna RNA)
          else()
            # Absolute path or full archive path
            set(VIENNARNA_LIBRARY "${_pc_lib}")
          endif()
        endif()
        if (PC_VIENNARNA_INCLUDE_DIRS AND NOT VIENNARNA_INCLUDE_DIRECTORY)
          list(GET PC_VIENNARNA_INCLUDE_DIRS 0 _pc_inc)
          set(VIENNARNA_INCLUDE_DIRECTORY "${_pc_inc}")
        endif()
      endif()
    endif()
  endforeach()
endif()

# ---------- Manual probing if still not found ----------
# Common prefixes
set(_COMMON_PREFIXES
  "$ENV{CONDA_PREFIX}"
  "$ENV{HOMEBREW_PREFIX}"
  "/opt/homebrew"
  "/usr/local"
  "/opt"
  "$ENV{HOME}/local"
)

# Library: if not already a full path from above, try to locate it
if (NOT VIENNARNA_LIBRARY)
  # Hints: user-provided dir plus common prefixes
  if (_LIB_HINTS)
    list(APPEND _COMMON_PREFIXES ${_LIB_HINTS})
  endif()
  list(REMOVE_DUPLICATES _COMMON_PREFIXES)

  find_library(VIENNARNA_LIBRARY
    NAMES vrna RNA
    HINTS ${_COMMON_PREFIXES}
    PATH_SUFFIXES lib lib64
    DOC "Path to ViennaRNA library")
endif()

# Include directory: must contain ViennaRNA/pair_mat.h
if (NOT VIENNARNA_INCLUDE_DIRECTORY)
  find_path(VIENNARNA_INCLUDE_DIRECTORY
    NAMES ViennaRNA/pair_mat.h
    HINTS ${_inc_hint} ${_COMMON_PREFIXES}
    PATH_SUFFIXES include
    DOC "Path to ViennaRNA headers")
endif()

# If VIENNARNA_LIBRARY is a bare name (e.g. "RNA"), try to resolve to an absolute path.
# If that fails, record likely link directories to use later.
set(VIENNARNA_LINK_DIRS "")
if (VIENNARNA_LIBRARY AND NOT IS_ABSOLUTE "${VIENNARNA_LIBRARY}")
  find_library(_VIENNA_ABS
    NAMES ${VIENNARNA_LIBRARY} vrna RNA
    HINTS
      $ENV{CONDA_PREFIX}
      $ENV{HOMEBREW_PREFIX}
      /opt/homebrew
      /usr/local
      /opt
    PATH_SUFFIXES lib lib64
  )

  if (_VIENNA_ABS)
    set(VIENNARNA_LIBRARY "${_VIENNA_ABS}" CACHE FILEPATH "Path to ViennaRNA library" FORCE)
  else()
    set(VIENNARNA_LINK_DIRS
      $ENV{CONDA_PREFIX}/lib
      $ENV{HOMEBREW_PREFIX}/lib
      /opt/homebrew/lib
      /usr/local/lib
      /opt/lib)
    list(REMOVE_DUPLICATES VIENNARNA_LINK_DIRS)
  endif()
endif()

message(STATUS "VIENNARNA_LIBRARY = ${VIENNARNA_LIBRARY}")
message(STATUS "VIENNARNA_INCLUDE_DIRECTORY = ${VIENNARNA_INCLUDE_DIRECTORY}")


# ---------- Friendly failure message ----------
if (NOT VIENNARNA_LIBRARY OR NOT VIENNARNA_INCLUDE_DIRECTORY)
  string(JOIN "\n" _help
    "❌ ViennaRNA not found."
    ""
    "This project needs the ViennaRNA library and headers to build."
    "If you already have ViennaRNA installed, CMake just couldn't find it with the default settings."
    ""
    "Quick fix:"
    "  Step 1: Locate the ViennaRNA library and header on your system."
    "  They are often named 'libRNA.a' (or similar) and 'pair_mat.h'."
    ""
    " 🔎 To locate these files, try:"
    "    find \$HOME /usr/local /opt -type f -name 'libRNA*' 2>/dev/null"
    "    find \$HOME /usr/local /opt -type f -name 'pair_mat.h' 2>/dev/null"
    ""
    "  Step 2: Set environment variables with the paths you found:"
    "    export VIENNARNA_LIBRARY=/path/to/libRNA.a"
    "    export VIENNARNA_INCLUDE_DIRECTORY=/path/to/include"
    ""
    "  For example, if you found:"
    "    /usr/local/lib/libRNA.a"
    "    /usr/local/include/ViennaRNA/pair_mat.h"
    "  then use:"
    "    export VIENNARNA_LIBRARY=/usr/local/lib/libRNA.a"
    "    export VIENNARNA_INCLUDE_DIRECTORY=/usr/local/include"
    ""
    "  Step 3: Re-run CMake:"
    "    cmake -S . -B build"
    ""
    "Optional: pass the paths directly to CMake instead of environment variables:"
    "    cmake -S . -B build \\"
    "      -DVIENNARNA_LIBRARY=/path/to/libRNA.a \\"
    "      -DVIENNARNA_INCLUDE_DIRECTORY=/path/to/include"
  )
  message(FATAL_ERROR "${_help}")
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

target_include_directories(KnotergyLib SYSTEM PUBLIC ${VIENNARNA_INCLUDE_DIRECTORY})

# Add link directories if we couldn't resolve an absolute path earlier.
if (VIENNARNA_LINK_DIRS)
  target_link_directories(KnotergyLib PUBLIC ${VIENNARNA_LINK_DIRS})
endif()

# ------------------------- Link Libraries -------------------------
target_link_libraries(KnotergyLib
  PUBLIC
    ${VIENNARNA_LIBRARY}
    nlohmann_json::nlohmann_json
    uni-algo::uni-algo
)

# ------------------------- Executable -------------------------
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE KnotergyLib)

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