#!/bin/bash

# Compress output files, excluding model state files and components. (i.e. Only include logs)

LABEL=$1

# Set root directories
if [ -d "/fastdata/ac1drb" ]; then
	FD="/fastdata/ac1drb"
else
	FD="/fastdata-sharc/ac1drb"
fi

NOW=$(date +%y.%m.%d_%H.%M)
INPUT_DIR=${FD}/output/$LABEL/
OUTPUT_DIR=${FD}

tar -zvcf $OUTPUT_DIR/output_$LABEL\_$NOW.tar.gz $INPUT_DIR \
--exclude 'model/*' \
--exclude 'model' \
--exclude 'run/*' \
--exclude 'run' \
--exclude 'output.script' \
--exclude 'run_brahms*'\

echo "Done!"
