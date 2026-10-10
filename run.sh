#!/bin/bash
# Run a calculation script with lunasrim (pure Python, no SRIM.exe needed).
#
#   ./run.sh                 runs scripts/H_in_C.py
#   ./run.sh scripts/F_in_Fe.py
#
# Results are written to output/<name>/ ("SRIM Outputs" holds the TRIM files).
# To run the original SRIM-2013 TRIM.exe on a TRIM.IN instead, see
# tools/run_real_trim.sh.
set -e
cd "$(dirname "$0")"
mkdir -p output
python3 "${1:-scripts/H_in_C.py}"
