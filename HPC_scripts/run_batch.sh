#!/bin/bash

# Set HPC parameters
export USER="ac1drb"
export RMEM="3G"
# export RMEM="30G"
export TIME="0:30:00"
# export TIME="1:30:00"

# Set model name and experiment number
export MODEL="Physical_2CH"
export STRIATUM="20.04.10_17.00_84900+849"
export EXP_NO="0"

# Set parallel job stride values
export T_START="1"
export T_STOP="1"
export T_STRIDE="1"

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

# TODO: Sanity check on batch variables
# TODO: Cleaner printout of batch variables
# TODO: Update batch variable code to be less fragile, iterate over any two varaibles

# Set number of channels, variation flags, and number of parallel jobs
while getopts ":c:f:m:w:" opt; do
	case $opt in
		c) CHANNELS=${OPTARG}
		;;
		f) 
			VAR_FSI=${OPTARG}
		  	if [ ! -z "${VAR_FSI}" ]; then
		  		export T_STOP=`echo ${T_STOP} \* 5 | bc`
				printf "TSTOP now ${TSTOP}"
		  	fi
		;;
		m)
			VAR_MSN=${OPTARG}
		  	if [ ! -z "${VAR_MSN}" ]; then
		  		export T_STOP=`echo ${T_STOP} \* 5 | bc`
				printf "TSTOP now ${TSTOP}"
		  	fi
		;;
		w) 
			VAR_WCH=${OPTARG}
		  	if [ ! -z "${VAR_WCH}" ]; then
		  		export T_STOP=`echo ${T_STOP} \* 5 | bc`
				printf "TSTOP now ${TSTOP}"
		  	fi
		;;
		\?) echo "Invalid option -${OPTARG}" >&2
		;;
	esac
done
	
# Array Jobs on ShARC and Iceberg can have a maximum of 75000 tasks
# https://docs.hpc.shef.ac.uk/en/latest/parallel/JobArray.html

# number of parallel jobs should be 5 * number of values to vary over
# uid = ceiling (current job / 5*(?)) - i think?

# Send jobs to SGE
echo "Executing model '${MODEL}' experiment #${EXP_NO} on striatum ${STRIATUM} with ${CHANNELS} channels"
qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -t ${T_START}:${T_STOP}:${T_STRIDE} -N ${MODEL} batch_submit.sge \
	-c${CHANNELS} -f${VAR_FSI} -m${VAR_MSN} -w${VAR_WCH}

