#!/bin/bash

# This script is for running multiple debug jobs without having to call
# `setup.sh <run> debug` for each run, waiting for each run to finish
# before you start the next. 
# 
# USAGE:
#     cd peakpatch
#     ./debugsetup.sh <run1> <run2> ... <runN>
# where the runs you want to submit are in direcotries `<run1>` etc.
 
# Make sure runs are passed at command line
if [ "$#" == 0 ]; then
    echo "You must pass runs at command line"
    exit
fi

# If hold flag passed, then jobs aren't submitted
if [ "$1" == "-h" ]; then
    shift
    for run in "$@"; do
        ./setup.sh $run hold
    done
    exit
fi

# Run debug jobs sequentially
for run in "$@"; do

    # Set up job
    echo "Preparing run $run."
    ./setup.sh $run debug

    # Read parameter file to get run name
    filename="$run/param/param.params"
    while read line; do
        line=${line%'#'*} # Trim any trailing comments beginning with '#'
        variable=`echo ${line%'='*}`
        value=`echo ${line#*'='}`

        # Get seed number
        if [ "$variable" = "seed" ]; then
            seed=$value

        # Get file name for run
        elif [ "$variable" = "short_name" ]; then
            short_name=${value%"'"*}       # trim trailing single quote
            short_name=${short_name#*"'"}  # trim leading single quote
            short_name=${short_name%"\""*} # trim trailing double quote
            short_name=${short_name#*"\""} # trim leading double quote
        fi
    done < $filename

    # Get the job ID by reading stdout file
    flag=false
    while [ "$flag" != true ]; do
        if [ -f "$run/${short_name}_${seed}_"*".txt" ]; then
            flag=true
        else
            sleep 15s
        fi
    done
    jobid=$(ls $run/${short_name}_${seed}_*.txt)
    jobid=${jobid%".txt"*}                        # trim trailing .txt
    jobid=${jobid#*"$run/${short_name}_${seed}_"} # trim leading stuff

    # Don't continue until the current job has finished running
    flag=false
    verbose=false
    while [ "$flag" != true ]; do

        # Read state of specified job from scheduler, check if the job has completed
        jobstate=$(sacct -j $jobid --format=State)

        if [[ "$run" == "${!#}" ]]; then
            flag=true

        # If job is still waiting in scheduler, wait 30 seconds and then try again
        elif [[ "$jobstate" == *"PENDING"* ]]; then
            if [ "$verbose" == true ]; then echo "job $jobid is pending"; fi
            sleep 15s

        # If job is running, wait 30 seconds and then try again
        elif [[ "$jobstate" == *"RUNNING"* ]]; then
            if [ "$verbose" == true ]; then echo "job $jobid is running"; fi
            sleep 15s

        # If there was a nodefail, warn the user and continue
        elif [[ "$jobstate" == *"NODE_FAIL"* ]]; then
            echo "job $run experienced a node fail."
            flag=true

        # If job was cancelled, warn the user and exit
        elif [[ "$jobstate" == *"CANCELLED"* ]]; then
            echo "job $run was cancelled, sorry."
            exit

        # If the job is completed, set flag to true to exit the loop
        elif [[ "$jobstate" == *"COMPLETED"* ]]; then
            echo "The run $run (jobid=$jobid) has completed."
            flag=true
        fi
    done

    if [[ "$run" != "${!#}" ]]; then
        echo
    else
        echo "Last job submitted."
    fi

done
