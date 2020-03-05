#!/bin/bash
#$ -m beas
#$ -M d.r.buxton@sheffield.ac.uk

# Reassign experiment variables
MODEL_DIR=$1
EXP_NO=$2
WORK_DIR=$3
OUTPUT_DIR=$4
CON_PROB=$5

# Set paths
S2B_DIR="${HOME}/SpineML_2_BRAHMS"

# Run experiment; -g prevents GUI, -s prevents extra xsltproc calls
if [ -z "$5" ]
	then
		${S2B_DIR}/convert_script_s2b -gs \
			-m ${MODEL_DIR} -e ${EXP_NO} -w ${WORK_DIR} -o ${OUTPUT_DIR} 
	else
		${S2B_DIR}/convert_script_s2b -gs \
			-m ${MODEL_DIR} -e ${EXP_NO} -w ${WORK_DIR} -o ${OUTPUT_DIR}  \
			# Make these populations imported from previous script as well
			-f "'CH1 input':'Striatum D1':0:${CON_PROB}" \
			-f "'CH1 input':'Striatum D1':1:${CON_PROB}" \
			-f "'CH1 input':'Striatum D2':0:${CON_PROB}" \
			-f "'CH1 input':'Striatum D2':1:${CON_PROB}"
fi
