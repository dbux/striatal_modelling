#!/bin/bash

# Set environment variables
# Recommended settings:
# Physical: RMEM=60G, TIME=4:00:00
# Statistical: RMEM=60G, TIME=3:00:00
RMEM="60G"
TIME="4:00:00"
MAIL="d.r.buxton@sheffield.ac.uk"
USER="ac1drb"

# Set model name and experiment number
export MODEL="Physical_1CH"
export STRIATUM="20.02.25_10.51_84900+834_1CH_0sep_0overlap"
# export STRIATUM="19.03.27_11.05_680+3_1CH_phys_0sep"
export EXP_NO="0"

####################

# Set root directories
if [ -d "/fastdata/${USER}" ]; then
	FD="/fastdata/${USER}"
else
	FD="/fastdata-sharc/${USER}"
fi

export LOGS_DIR="${HOME}/logs"
export MODEL_ROOT="/data/${USER}/models"
export OUTPUT_ROOT="${FD}/output"
export SCRIPTS_DIR="${HOME}/scripts"
export WORK_ROOT="${FD}/temp"

# Clear previous data
echo "Clearing old logs…"
rm ${LOGS_DIR}/*
echo "Clearing old output…"
rm -rf ${OUTPUT_ROOT}

# Set model and connection list directories
export MODEL_DIR="${MODEL_ROOT}/${MODEL}"
export OUTPUT_DIR="${OUTPUT_ROOT}/${MODEL}"
export WORK_DIR="${WORK_ROOT}/${MODEL}"

# Execute model on SGE
echo "Sending model '${MODEL}' experiment #${EXP_NO} to Sun Grid Engine…"
qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -N ${MODEL} \
	-e ${LOGS_DIR}/${MODEL}.err -o ${LOGS_DIR}/${MODEL}.out \
	${SCRIPTS_DIR}/s2b.sh ${MODEL_DIR} ${EXP_NO} ${WORK_DIR} ${OUTPUT_DIR} 
