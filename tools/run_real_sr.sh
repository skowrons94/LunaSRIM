#!/bin/sh
# Run the original SRIM-2013 "SR Module" (SRModule.exe, wine, costrouc/pysrim
# docker image) on an SR.IN file; the output file it names is copied back.
#
# usage: tools/run_real_sr.sh <directory containing SR.IN>
set -e
if [ -z "$1" ] || [ ! -f "$1/SR.IN" ]; then
    echo "usage: $0 <directory containing SR.IN>" >&2
    exit 1
fi
d=$(cd "$1" && pwd)
docker run --rm -e SR_TIMEOUT="${SR_TIMEOUT:-600}" -v "$d":/case costrouc/pysrim sh -c '
cd "/tmp/srim/SR Module" && cp /case/SR.IN SR.IN && touch /tmp/stamp &&
timeout ${SR_TIMEOUT:-600} xvfb-run -a wine SRModule.exe >/dev/null 2>&1;
find . -maxdepth 1 -type f -newer /tmp/stamp ! -name SR.IN -exec cp {} /case/ \;'
