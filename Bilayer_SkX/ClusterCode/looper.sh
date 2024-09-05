#!/bin/bash
string_to_search_for="run_MPI_Aminus"
current_dir=$(pwd)
for FILE in $(ls $current_dir) 
do
	if [[ ${FILE} =~ ${string_to_search_for} ]] 
 	then
  		dos2unix ${FILE}
		echo ${FILE}
		sbatch ${FILE}
		sleep 1 # pause to be kind to the scheduler
	fi
done

#!/bin/bash
