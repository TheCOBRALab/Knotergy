#!/usr/bin/env bash

# ──────────────────────────────────────────────────────────────
# 0 · Compiler flags
# ──────────────────────────────────────────────────────────────
export CFLAGS="${CFLAGS:-}  -std=c11"
export CXXFLAGS="${CXXFLAGS:-} -std=c++17 -include cstdint"

if [[ "$OSTYPE" == "darwin"* ]]; then
  export CXXFLAGS="${CXXFLAGS}"
else
  export CPPFLAGS="${CPPFLAGS:-} -D_POSIX_C_SOURCE=200809L -D_GNU_SOURCE"
fi 


# ──────────────────────────────────────────────────────────────
# 1 · Fresh build directory
# ──────────────────────────────────────────────────────────────a
rm -rf build
mkdir  build
cd     build


# ──────────────────────────────────────────────────────────────
# 2 · CMake phase for your own project that links ViennaRNA
# ──────────────────────────────────────────────────────────────
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="${CC}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_C_FLAGS="-DHAVE_STRDUP=1" \
  -DVIENNARNA_LIB="${PREFIX}" \
  -DVIENNARNA_INCLUDE_DIR="${PREFIX}/include" \
  -DBUILD_TESTING=OFF

cmake --build . --parallel -- VERBOSE=1
cmake --install . --prefix "${PREFIX}"