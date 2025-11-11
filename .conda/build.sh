#!/usr/bin/env bash

mkdir build
cd build

cmake ..                                                                       \
  -DCMAKE_CXX_COMPILER=${BUILD_PREFIX}/bin/$(basename ${CXX})                  \
  -DCMAKE_INSTALL_PREFIX=${PREFIX}                                             \
  -DCMAKE_BUILD_TYPE=Release                                                   \
  -DTESTS=ON                                                                   \
  -DEXAMPLES=OFF                                                               \
  -DDocumentation=OFF

make -j2 VERBOSE=1
ctest --output-on-failure
make install
