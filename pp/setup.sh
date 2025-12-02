#!/bin/bash

# To run:
#     cd peakpatch
#     ./setup.sh <run-name>
# Where the run you want to submit is in directory `<run-name>`
# 
# For small test runs, you can submit to the Niagara debug queue with
#     ./setup.sh <run-name> debug

# Configure Peak Patch
cd $SCRATCH/peakpatch
./configure
source ~/.bashrc

if [ "$1" == '' ]; then
    exit 0 # Exits if no directory passed at command line
fi

# Load modules on Niagara supercomputer
module load NiaEnv/2019b intel/2019u4 intelmpi/2019u4 fftw/3.3.8 gsl/2.5 cfitsio/3.450 python/3.9.8 #intelpython3/2019u5

# Setting up the run
echo "Compiling Peak Patch run"
cd $1
python3 $SCRATCH/peakpatch/python/peak-patch.py ./param/param.params

# Making collapse tables
echo "Making collapse tables"
./bin/hpkvd 1

# Get the name of the job script from parameter file
filename="param/param.params"
while read line; do # read each line in param/param.params

    line=${line%'#'*} # Trim any trailing comments beginning with '#'
    variable=`echo ${line%'='*}`
    value=`echo ${line#*'='}`

    # Get seed number
    if [ "$variable" = "seed" ]; then
        seed=$value

    # Get number of nodes requested
    elif [ "$variable" = "nnodes" ]; then
        nnodes=$value

    # Get time limit requested
    elif [ "$variable" = "tlimit" ]; then
        tlimit=${value%"'"*}   # trim trailing single quote
        tlimit=${tlimit#*"'"}  # trim leading single quote
        tlimit=${tlimit%"\""*} # trim trailing double quote
        tlimit=${tlimit#*"\""} # trim leading double quote
        h=${tlimit%":"*};h=${h%":"*}
        m=${tlimit%":"*};m=${m#*":"}
        s=${tlimit#*":"};s=${s#*":"}
        t=$(( 3600*$h + 60*$m + $s ))

    # Get file name for run
    elif [ "$variable" = "short_name" ]; then
        short_name=${value%"'"*}       # trim trailing single quote
        short_name=${short_name#*"'"}  # trim leading single quote
        short_name=${short_name%"\""*} # trim trailing double quote
        short_name=${short_name#*"\""} # trim leading double quote
    fi
done < $filename
job_script="${short_name}_${seed}.sh" # name of the job script

# Check job node and time limits
if ((nnodes>1000)) || ((t>86400)); then
    echo "ABORTING: Niagara jobs must use 1000 or fewer nodes and take no more than 24 hours."
    exit
fi

# Submit job to the Niagara scheduler
if [ "$2" = "" ]; then
    echo "Submitting $job_script to Niagara scheduler"
    sbatch $job_script

# Prepare job script without submitting to scheduler
elif [ "$2" = "hold" ]; then
    echo "Job $job_script ready to submit"

# Submit debug job to the Niagara scheduler
elif [ "$2" = "debug" ]; then
    if ((nnodes>4)) || ((t>3600)); then
        echo "ABORTING: debug jobs must use 4 nodes or fewer and take no more than 1 hour."
    else 
        echo "Submitting $job_script to Niagara debug queue"
        sed -i '1 a\#SBATCH -p debug' $job_script
        sbatch $job_script
    fi

else
    echo "Command line argument \"$2\" not recognized."
    echo "Job $job_script ready to submit"
fi
