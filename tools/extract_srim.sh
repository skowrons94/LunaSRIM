#!/bin/sh
# Copy the SRIM-2013 installation (TRIM.exe, SRModule.exe, data files) out
# of the costrouc/pysrim docker image, e.g. to look at the binaries or to
# refresh lunasrim/data.
#
# usage: tools/extract_srim.sh [destination directory, default ./srim]
set -e
dest=${1:-srim}
cid=$(docker create costrouc/pysrim)
trap 'docker rm "$cid" >/dev/null' EXIT
docker cp "$cid":/tmp/srim "$dest"
echo "SRIM copied to $dest"
ls "$dest"
