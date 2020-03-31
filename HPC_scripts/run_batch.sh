#!/bin/bash

# Set HPC parameters
export USER="ac1drb"
export RMEM="50G"
export TIME="6:00:00"

# Set model name and experiment number
export MODEL="Physical_1CH"
export STRIATUM="20.03.27_12.09_84900+864_1CH_0sep"
export EXP_NO="0"

####################

# Set root directories
if [ -d "/fastdata/${USER}" ]; then
	export FD="/fastdata/${USER}"
else
	export FD="/fastdata-sharc/${USER}"
fi
export LISTS_ROOT="/data/${USER}/striatums"
export LOGS_DIR="${HOME}/logs"
export MODEL_ROOT="/data/${USER}/models"
export OUTPUT_ROOT="${FD}/output"
export S2B_DIR="${HOME}/SpineML_2_BRAHMS"
export SCRIPTS_DIR="${HOME}/scripts"
export WORK_ROOT="${FD}/temp"


# Clear previous data
echo "Clearing old logs…"
rm ${LOGS_DIR}/*
echo "Clearing old temp files…"
rm -rf ${WORK_ROOT}
echo "Clearing old output…"
rm -rf ${OUTPUT_ROOT}

# Set model and connection list directories
export MODEL_DIR="${MODEL_ROOT}/${MODEL}"
export LISTS_DIR="${LISTS_ROOT}/${STRIATUM}/connection_lists"

# Send jobs to SGE
qsub -V -rmem=${RMEM} -h_rt=${TIME} batch_submit.sge

