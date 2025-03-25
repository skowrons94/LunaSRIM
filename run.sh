#!/bin/bash

# If output directory does not exist, create it
if [ ! -d "output" ]; then
  mkdir output
fi

docker run -v $PWD/scripts:/opt/pysrim/ \
           -v $PWD/output/:/tmp/output \
           -it costrouc/pysrim sh -c "xvfb-run -a python3.6 /opt/pysrim/H_in_C.py"