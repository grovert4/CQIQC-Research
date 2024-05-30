#!/bin/bash
string_to_search_for="run_tempering"
current_dir=$(pwd)
for FILE in $(ls $current_dir); do
	if [[ ${FILE} =~ ${string_to_search_for} ]]; then
		echo ${FILE}
		# sbatch -p ccq -C rome ${FILE}
	sleep 1 # pause to be kind to the scheduler
	fi
done
#!/bin/bas
