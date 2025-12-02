#!/bin/bash

# Usage:
# 
#     cd peakpatch
#     ./scripts/pksc2npz.sh /path/to/peak/patch/run/

# Load in the run directory (get absoltue path if relative path given)
run_dir="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"

# Load in time limit
if [ "$2" = "" ]; then
    time_limit=1:00:00
else
    time_limit=$2
fi

# Load in the catalogue file (get absoltue path if relative path given)
if [ "$3" = "" ]; then
    catalogue_file="$(ls ${run_dir}/output/*merge*)"
else
    catalogue_file="$(cd "$(dirname "$2")";pwd)/$(basename "$2")"
fi

# Load in the number of parallel cores for the run
if [ "$4" = "" ]; then
    parallel_cores=40
else
    parallel_cores=$3
fi

# Niagara batch job parameters
job_name=pksc2npz
bash_output=./pksc2npz.out

# For jobs less than 1 hour and using 4 nodes or less, run in debug queue
h=$(echo "$time_limit" | cut -d':' -f1)
m=$(echo "$time_limit" | cut -d':' -f2)
s=$(echo "$time_limit" | cut -d':' -f3)
h=$(( h + m/60 + s/3600 ))
if [ "$h" -le 1 ]; then
    debug="#SBATCH -p debug"
else
    debug=""
fi

# Change directories to the halo catalogue's directory
cd $run_dir/output/.

# Make a copy the Niagara batch job template script
cp $PP_DIR/templates/pksc2npz_niagara_submit_script.sh ./pksc2npz.sh

# Set job script parameters
sed -i "s!DEBUG!${debug}!g" ./pksc2npz.sh
sed -i "s!RUN_DIR!${run_dir}!g" ./pksc2npz.sh
sed -i "s!CATALOGUE_FILE!${catalogue_file}!g" ./pksc2npz.sh
sed -i "s!PARALLEL_CORES!${parallel_cores}!g" ./pksc2npz.sh
sed -i "s!TIME_LIMIT!${time_limit}!g" ./pksc2npz.sh
sed -i "s!JOB_NAME!${job_name}!g" ./pksc2npz.sh
sed -i "s!BASH_OUTPUT!${bash_output}!g" ./pksc2npz.sh
sed -i "s!PP_DIR!${PP_DIR}!g" ./pksc2npz.sh

# Submit job to Niagara scheduler
sbatch ./pksc2npz.sh
