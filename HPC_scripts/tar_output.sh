#!/bin/bash

# Compress output files, excluding model state files and components. (i.e. Only include logs)

# LABEL=$1

# Set root directories
#if [ -d "/fastdata/ac1drb" ]; then
#	FD="/fastdata/ac1drb"
	#else
#	FD="/fastdata-sharc/ac1drb"
	#fi

FD="/fastdata/ac1drb"
NOW=$(date +%y.%m.%d_%H.%M)
INPUT_DIR=${FD}/output
OUTPUT_DIR=${FD}

# From https://stackoverflow.com/a/2108296
for dir in ${INPUT_DIR}/*/		# List directories in the form "/tmp/dirname/"
do
	dir=${dir%*/}     		 	# Remove the trailing "/"
	echo "Dir is"
	#echo ${dir##*/}				# Print everything after the final "/"
	exp=${dir##*/}
	echo ${exp}
	#echo "Dir is"
	#echo ${exp}
	
	#echo "other dir is"
	#echo ${dir##*/}
	
	tar -C ${INPUT_DIR} -zvcf ${OUTPUT_DIR}/output_${exp}\_${NOW}.tar.gz ${exp} \
	--exclude 'model/*' \
	--exclude 'model' \
	--exclude 'run/*' \
	--exclude 'run' \
	--exclude 'output.script' \
	--exclude 'run_brahms*'
	
done

#tar -C ${INPUT_DIR} -zvcf ${OUTPUT_DIR}/output_${LABEL}\_${NOW}.tar.gz ${LABEL} \
#--exclude 'model/*' \
#--exclude 'model' \
#--exclude 'run/*' \
#--exclude 'run' \
#--exclude 'output.script' \
#--exclude 'run_brahms*'

echo "Done!"