#!/bin/sh
# Build and run the native bit-exactness test without make:  CC=gcc sh test.sh   (default cc)
set -e
cd "$(dirname "$0")"
CC=${CC:-cc}
mkdir -p build
$CC -std=c99 -O2 -Wall -Wextra -Ialg -Ivectors -o build/test_native_madow alg/*.c tests/test_native_madow_main.c
./build/test_native_madow | tail -n 6
python3 tools/check_madow_vectors.py
