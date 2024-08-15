#!/bin/bash

# This file builds, runs, executse test.
# Usages
# 1. `./do build`
#    Builds all files
# 2. `./do run <target>`
#    Builds the target and run the binary.
# 3. `./do test`
#    Executes all test.
# 4. `./do test <target>`
#    Executes all test in the `<target>_test` file.
# 5. `./do draw`
#    Draws super singular isogeny graph

set -eux

# Buils executable, it has optional 1 following argument
# 1. target: if this argment is set, only build this executable
#
# e.g. build CGLHash
function build() {
    local target_option=""
    if [ $# -eq 1 ]; then
        target_option="--target $1"
    fi

    cmake -S . -B build
    cmake --build build $target_option
}

cmd=$1

if [ $cmd = "build" ]; then
    build
elif [ $cmd = "run" ]; then
    target=$2
    build $target
    ./build/src/$target
elif [ $cmd = "test" ]; then
    if [ $# -eq 1 ]; then
        build 
        cd build/test && ctest
    else 
        target=$2
        build "${target}_test"
        cd build/test && ./${target}_test
    fi
elif [ $cmd = "python" ]; then
    target=$2
    cd operation-calc
    poetry run python3 exec/$target.py
elif [ $cmd = "draw" ]; then
    build drawSSIG
    ./build/src/drawSSIG > graph.dot
    dot -Tpng -o graph.png graph.dot
fi