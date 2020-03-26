#!/bin/bash

# Go to simulation working directory
cd ${WORK_DIR}

# Check for existence of model.xml                                                            
if [ ! -e model.xml ]; then                                                
	echo "No model.xml found in $(pwd); aborting."
else
	echo "In ${WORK_DIR}:"
	
	# For every .bin file
	for file in *.bin; do

		# (Break out of the loop if no .bin files exist)                                                
		[ -f "$file" ] || break                                             
		FILE_CSV=$(basename $file .bin).csv

		# If a corresponding CSV exists
		if [ -f ${FILE_CSV} ]; then  
			# Get the number of CSV rows                                                
			CSV_LINES=$(awk 'END{ print NR }' ${FILE_CSV})                  
			echo "(✓) ${FILE_CSV} found (${CSV_LINES} connections)"

			# Construct the string occurring immediately before the number of connections
			STR_MATCH=${file}'" num_connections="'

			# Correct the number of connections in model.xml                         
			sed -i "s/\(${STR_MATCH}\)[0-9]\+/\1${CSV_LINES}/" model.xml    
		else
			# Report if no matching CSV found
			echo "(!) ${FILE_CSV} not found"                                
		fi
	done
fi
