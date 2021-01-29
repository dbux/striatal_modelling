#!/bin/bash

# Set HPC parameters
export USER="ac1drb"
export RMEM="30G"
# export TIME="1:30:00"
export TIME="8:00:00"

# Set model name and experiment number
export MODEL="Universal"
# export MODEL_NEW="Physical_2CH"
# export MODEL_OLD="Physical_2CH_old"
# export STRIATUM="20.04.10_17.00_84900+849"
export LRG="21.01.25_20.06_84900+845"
export MED="21.01.29_11.24_10612+106"
export SML="21.01.28_12.53_1326+15"
export STAT="21.01.26_09.51_6000+60"
export DELAY="1.0"
export CHANNELS=1
export EXP_NO=1

# Set parallel job stride values (modified later in script)
export T_START="1"
export T_STOP="61"
export T_STRIDE="60"

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


# Send jobs to SGE
# Large phys model
export STRIATUM=${LRG}
echo "Executing experiment #${EXP_NO} on striatum ${STRIATUM} (${CHANNELS} channels)"
export LISTS_DIR="${LISTS_ROOT}/${STRIATUM}/connection_lists"
qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -t ${T_START}:${T_STOP}:${T_STRIDE} -N "Striatum_LRG" pair_submit.sge

# Medium phys model
export STRIATUM=${MED}
echo "Executing experiment #${EXP_NO} on striatum ${STRIATUM} (${CHANNELS} channels)"
export LISTS_DIR="${LISTS_ROOT}/${STRIATUM}/connection_lists"
qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -t ${T_START}:${T_STOP}:${T_STRIDE} -N "Striatum_MED" pair_submit.sge

# Small phys model
export STRIATUM=${SML}
echo "Executing experiment #${EXP_NO} on striatum ${STRIATUM} (${CHANNELS} channels)"
export LISTS_DIR="${LISTS_ROOT}/${STRIATUM}/connection_lists"
qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -t ${T_START}:${T_STOP}:${T_STRIDE} -N "Striatum_SML" pair_submit.sge

# Statistical model
export STRIATUM=${STAT}
echo "Executing experiment #${EXP_NO} on striatum ${STRIATUM} (${CHANNELS} channels)"
export LISTS_DIR="${LISTS_ROOT}/${STRIATUM}/connection_lists"
qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -t ${T_START}:${T_STOP}:${T_STRIDE} -N "Striatum_STAT" pair_submit.sge