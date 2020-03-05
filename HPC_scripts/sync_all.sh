#!/bin/bash

# Set paths
DATA_DIR="${HOME}/Dropbox/University/2019 (Sheffield) - Striatal oscillations"

echo "Syncing models from ${DATA_DIR}/Models…"
rsync -ahirstuvz "${DATA_DIR}/Models/" "ac1drb@iceberg.sheffield.ac.uk:/data/ac1drb/models"

echo "Syncing striatums from ${DATA_DIR}/Striatums…"
rsync -ahirstuvz "${DATA_DIR}/Striatums/" "ac1drb@iceberg.sheffield.ac.uk:/data/ac1drb/striatums"

echo "Syncing bash scripts from ${DATA_DIR}/Scripts…"
rsync -ahirstuvz "${DATA_DIR}/Scripts/" "ac1drb@iceberg.sheffield.ac.uk:/home/ac1drb/scripts"
rsync -ahirstuvz "ac1drb@iceberg.sheffield.ac.uk:/home/ac1drb/scripts/" "${DATA_DIR}/Scripts/"

echo "Syncing MatLab scripts from ${DATA_DIR}/MatLab…"
rsync -ahirstuvz "${DATA_DIR}/MatLab/" "ac1drb@iceberg.sheffield.ac.uk:/home/ac1drb/MatLab"
