#!/bin/sh
# Run the original SRIM-2013 TRIM.exe (wine, costrouc/pysrim docker image)
# on a TRIM.IN file and copy all output files back.
#
# usage: tools/run_real_trim.sh <directory containing TRIM.IN>
#
# The output files (RANGE.txt, ..., TRANSMIT.txt, ...) are written into the
# same directory; TRIM.IN.used is the input as TRIM saved it.
set -e
if [ -z "$1" ] || [ ! -f "$1/TRIM.IN" ]; then
    echo "usage: $0 <directory containing TRIM.IN>" >&2
    exit 1
fi
d=$(cd "$1" && pwd)
docker run --rm -e TRIM_TIMEOUT="${TRIM_TIMEOUT:-3000}" -v "$d":/case costrouc/pysrim sh -c '
cd /tmp/srim && cp /case/TRIM.IN TRIM.IN && printf 1 > TRIMAUTO && rm -f "SRIM Outputs"/*.txt *.txt;
timeout ${TRIM_TIMEOUT:-3000} xvfb-run -a wine TRIM.exe >/dev/null 2>&1;
cp "SRIM Outputs"/*.txt /case/ 2>/dev/null; cp *.txt /case/ 2>/dev/null; cp TRIM.IN /case/TRIM.IN.used'
