string_to_search_for="run_tempering"
$(cd '.\Tempering Inputs\')
current_dir=$(pwd)
for FILE in $(ls $current_dir); do
	if [[ ${FILE} =~ ${string_to_search_for} ]]; then
		echo ${FILE}
		# sbatch -p ccq -C rome ${FILE}
	sleep 1 # pause to be kind to the scheduler
	fi

done


#!/bin/bash

# Get the current directory.
current_dir=$(pwd)

