#!/bin/bash

# Set HPC parameters
export USER="ac1drb"
export RMEM="30G"
export TIME="1:30:00"

# Set model name and experiment number
export MODEL="Physical_2CH_old"
export STRIATUM="20.04.10_17.00_84900+849"
export CHANNELS=1

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
echo -n	"Clearing logs…"
rm -rf	${LOGS_DIR}/*
echo -n	" temp files…"
rm -rf	${WORK_ROOT}
echo -n	" output…"
rm -rf	${OUTPUT_ROOT}
echo	" done!"

# Set model and connection list directories
export MODEL_DIR="${MODEL_ROOT}/${MODEL}"
export LISTS_DIR="${LISTS_ROOT}/${STRIATUM}/connection_lists"

# Set number of channels, variation flags, and number of parallel jobs
while getopts ":c:f:m:w:" opt; do
	case $opt in
		c) 
			export CHANNELS=${OPTARG}
			if [ "${CHANNELS}" -eq "1" ]; then
				export EXP_NO=1
		  	elif [ "${CHANNELS}" -eq "2" ]; then
				export EXP_NO=0
			fi
		;;
		f) 
			export VAR_FSI=${OPTARG}
		  	if [ ! -z "${VAR_FSI}" ]; then
		  		export T_STOP=`echo ${T_STOP} \* ${NUM_FSI} | bc`
		  	fi
		;;
		m)
			export VAR_MSN=${OPTARG}
		  	if [ ! -z "${VAR_MSN}" ]; then
		  		export T_STOP=`echo ${T_STOP} \* ${NUM_MSN} | bc`
		  	fi
		;;
		w) 
			export VAR_WCH=${OPTARG}
		  	if [ ! -z "${VAR_WCH}" ]; then
		  		export T_STOP=`echo ${T_STOP} \* ${NUM_WCH} | bc`
		  	fi
		;;
		\?) echo "Invalid option -${OPTARG}" >&2
		;;
	esac
done

# Send jobs to SGE
echo "Executing model '${MODEL}' experiment #${EXP_NO} on striatum ${STRIATUM} (${CHANNELS} channels)"
qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -t ${T_START}:${T_STOP}:${T_STRIDE} -N ${MODEL} batch_submit.sge