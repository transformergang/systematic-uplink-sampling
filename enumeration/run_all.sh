#!/bin/sh
# Runs all enumeration checks (a few seconds each). Needs Python 3 with NumPy and SciPy.
# usage: sh run_all.sh [python-executable]
PY=${1:-python3}
cd "$(dirname "$0")" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
set -e
$PY enum_s.py
$PY test_madow_exact.py | tail -n 5
$PY ../code/python/test_reference.py
