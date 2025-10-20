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
# 2 · Unpack ViennaRNA, patch and build ViennaRNA
# ──────────────────────────────────────────────────────────────
curl -L -O https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.0/ViennaRNA-2.7.0.tar.gz
tar -xf  ./ViennaRNA-2.7.0.tar.gz
cd   ViennaRNA-2.7.0

# Bug fixes in dlib until ViennaRNA updates to a newer version
if [[ "$OSTYPE" == "darwin"* ]]; then
    # “template” keyword.  Add the empty <> that upstream already has.
    sed -i.bak 's/::template go(/::template go<>(/' \
        src/dlib-19.24/dlib/global_optimization/find_max_global.h
fi

./configure --without-perl --prefix="${PREFIX}"
make  -j"${CPU_COUNT}"

cd .. # back to build dir

# ──────────────────────────────────────────────────────────────
# 4 · CMake phase for your own project that links ViennaRNA
# ──────────────────────────────────────────────────────────────
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="${CC}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_C_FLAGS="-DHAVE_STRDUP=1" \
  -DVIENNARNA_LIB="${PREFIX}" \
  -DVIENNARNA_INCLUDE_DIR="${PREFIX}/include" \
  -DENABLE_TESTS=OFF

cmake --build . --parallel -- VERBOSE=1
cmake --install . --prefix "${PREFIX}"
