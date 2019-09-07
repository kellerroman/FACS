#!/bin/bash

set -euo pipefail

rm -rf build
mkdir -p build && cd build

# Configure
cmake -DCODE_COVERAGE=ON -DPFUNIT=ON -DCMAKE_BUILD_TYPE=Debug ..
# Build (for Make on Unix equivalent to `make -j $(nproc)`)
cmake --build . --config Debug -- -j $(nproc) VERBOSE=1
# Test
ctest -j $(nproc) --output-on-failure
