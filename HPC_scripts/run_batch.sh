#!/bin/bash

# Set environment variables
RMEM="60G"
TIME="4:00:00"
# MAIL="d.r.buxton@sheffield.ac.uk"
USER="ac1drb"

# Set model name and experiment number
export MODEL="Physical_1CH"
export STRIATUM="19.07.29_11.52_84900+827_1CH_0sep"
export EXP_NO="0"

####################

# Set root directories
if [ -d "/fastdata/${USER}" ]; then
	FD="/fastdata/${USER}"
else
	FD="/fastdata-sharc/${USER}"
fi

export LISTS_ROOT="/data/${USER}/striatums"
export LOGS_DIR="${HOME}/logs"
export MODEL_ROOT="/data/${USER}/models"
export OUTPUT_ROOT="${FD}/output"
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

# Model construction and execution loop
cd $LISTS_DIR
# For each descriptor variant
for i in {0..90..10}; do 
	for j in {0..90..10}; do 

		# Create UID for each simulation (must match the naming pattern of connection lists)         
		if [ -z "$i" ]; then
			MSN_UID="bkMSN0"
		else                                    
			MSN_UID="bkMSN${i}" 
		fi

		if [ -z "$j" ]; then
			FSI_UID="bkFSI0"
		else
			FSI_UID="bkFSI${j}"
		fi 
		SIM_UID="${MSN_UID}_${FSI_UID}"                                             

		# Get connection probability
		# CON_PROB=$(bc <<< "scale=1;1 - ${i} / 100")							

		# Set unique per-simulation work and output directories
		export WORK_DIR="${WORK_ROOT}/${MODEL}/${SIM_UID}"
		export OUTPUT_DIR="${OUTPUT_ROOT}/${MODEL}/${SIM_UID}"
		echo ""
		echo "Working directory is ${WORK_DIR}"
		echo "Output directory is ${OUTPUT_DIR}"

		# Check for directory existence
		if [ ! -d ${WORK_DIR} ]; then
		    #echo "Making work dir"
		    mkdir -p ${WORK_DIR}
			else echo "WARNING! Working directory already exists, proceeding anyway"
		fi

		# Create output directory
		if [ ! -d ${OUTPUT_DIR} ]; then
		    mkdir -p ${OUTPUT_DIR}
			else echo "WARNING! Output directory already exists, proceeding anyway"
		fi

		# For each connection list file containing experiment UID
		echo "Copying connection lists…"
		for file in $(find . -maxdepth 1 -name "*${SIM_UID}*"); do 
			# Remove initial ./             
		    FILE_FULL=$(echo $file | sed "s|^\./||")   
			# Remove MSN/FSI-specific descriptors                         
		    FILE_TRIM=$(echo $FILE_FULL | sed "s/_bkMSN[0-9]\+_bkFSI[0-9]\+//") 
			# Copy connection list to work directory with unique descriptors removed			
		    cp $file ${WORK_DIR}/${FILE_TRIM}									
		done
		
		# Copy other model files to work directory
		cp -n ${MODEL_DIR}/* ${WORK_DIR}
		# Rewrite model connection lists	
		echo "Modifying model.xml…"						            
		${SCRIPTS_DIR}/modify_connections.sh								      

		# Execute model on SGE
		echo "Sending model '${MODEL}' experiment #${EXP_NO} version ${SIM_UID} to Sun Grid Engine…"
		qsub -V -l rmem=${RMEM} -l h_rt=${TIME} -N ${MODEL}_${SIM_UID} \
			-e ${LOGS_DIR}/${MODEL}_${SIM_UID}.err -o ${LOGS_DIR}/${MODEL}_${SIM_UID}.out \
			${SCRIPTS_DIR}/s2b.sh ${WORK_DIR} ${EXP_NO} ${WORK_DIR} ${OUTPUT_DIR} ${CON_PROB}
	done
done
