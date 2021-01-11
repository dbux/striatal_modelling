#!/bin/bash

# Set HPC parameters
export USER="ac1drb"
export RMEM="60G"
export TIME="1:30:00"

# Set model name and experiment number
export MODEL="Physical_2CH"
# export MODEL_NEW="Physical_2CH"
# export MODEL_OLD="Physical_2CH_old"
# export STRIATUM="20.04.10_17.00_84900+849"
export STRIATUM_HI-R="21.01.08_15.30_84900+811(r0.95)"
export STRIATUM_LO-R="21.01.08_16.16_84900+811(r0.05)"
export CHANNELS=1
export EXP_NO=1

# Set parallel job stride values (modified later in script)
export T_START="1"
export T_STOP="1"
export T_STRIDE="1"

# Number of variations in FSI / MSN / WCH variables
export NUM_FSI="6"
export NUM_MSN="5"
export NUM_WCH="5"

####################

# Set root directories
export FD="/fastdata/${USER}"
export LISTS_ROOT="/data/${USER}/striatums"
export LOGS_DIR="${HOME}/logs"
export MODEL_ROOT="/data/${USER}/models"
export OUTPUT_ROOT="${FD}/output"
export S2B_DIR="${HOME}/SpineML_2_BRAHMS"
export SCRIPTS_DIR="${HOME}/scripts"
export WORK_ROOT="${FD}/temp"

# Clear previous data
echo -n	"Clearing logs…"
rm -rf	${LOGS_DIR}/*
echo -n	" temp files…"
rm -rf	${WORK_ROOT}
echo -n	" output…"
rm -rf	${OUTPUT_ROOT}
echo	" done!"

# Set model and connection list directories
export MODEL_DIR="${MODEL_ROOT}/${MODEL}"

# export STRIATUM=${STRIATUM_HI-R}
# export LISTS_DIR="${LISTS_ROOT}/${STRIATUM}/connection_lists"

# Send jobs to SGE
echo "Executing model '${MODEL}' experiment #${EXP_NO} on striatum ${STRIATUM_HI-R} (${CHANNELS} channels)"
export STRIATUM=${STRIATUM_HI-R}
export LISTS_DIR="${LISTS_ROOT}/${STRIATUM}/connection_lists"
qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -t ${T_START}:${T_STOP}:${T_STRIDE} -N "Physical_NEW" pair_submit.sge

echo "Executing model '${MODEL}' experiment #${EXP_NO} on striatum ${STRIATUM_LO-R} (${CHANNELS} channels)"
export STRIATUM=${STRIATUM_LO-R}
export LISTS_DIR="${LISTS_ROOT}/${STRIATUM}/connection_lists"
qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -t ${T_START}:${T_STOP}:${T_STRIDE} -N "Physical_OLD" pair_submit.sge