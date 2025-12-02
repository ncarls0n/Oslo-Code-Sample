#!/bin/sh

# This script sets up map-making jobs for the Niagara supercomputer. For a given Peak Patch
# lightcone run `<run_name>` (i.e. a run with ievol=1 in the parameter file
# `<run_name>/param/param.params`), the script reads the parameter file and makes jobscripts for MPI
# runs of pks2cmb.f90. The final maps are output as .fits files to the directory `<run_name>/maps`
# and can be plotted using the script `mapview.py`.
#
# USAGE: first run peak patch:
#
#     cd peak-patch
#     ./setup.sh <run_name>
#
# Then run this mapping setup script
#
#     ./mapsetup.sh <run_name>
# 
# For smaller Peak Patch halo catalogues, making maps is usually very unintensive. For such runs,
# you can sumbit to the debug queue by changing the final command to
# 
#     ./mapsetup.sh <run_name> debug
# 
# Sometimes when debugging, you may want to make edits to the run directory after running Peak Patch
# and before trying to make maps. For this reason, I've introduced the flag `recompile', which
# recompiles pks2map.f90 before submitting map making scripts to the Niagara scheduler.
# 
#     ./mapsetup.sh <run_name> recompile

# Configure Peak Patch
cd $SCRATCH/peakpatch
./configure
source ~/.bashrc

# Load modules on Niagara supercomputer
module load NiaEnv/2019b intel/2019u4 intelmpi/2019u4
echo

# Go to maps directory
cd $1/maps

# Delete existing standard output and error files
rm -f ./*.stdout ./*.stderr

# Get the name of the job script from parameter file
echo "Reading parameter file"
while read line; do # read each line in <run-name>/param/param.params

    line=${line%"#"*}            # Trim trailing comments beginning with #
    variable=`echo ${line%"="*}` # Trim end of `line` beginning with `=`
    value=`echo ${line#*"="}`    # Trim begining of `line` ending in `=`

    # Get seed number
    if [ "$variable" = "seed" ]; then
        seed=$value

    # Get "short_name" used in job submission scripts
    elif [ "$variable" = "short_name" ]; then
        short_name=${value%"'"*}      # trim trailing single quote
        short_name=${short_name#*"'"} # trim leading single quote
        short_name=${short_name%'"'*} # trim trailing double quote
        short_name=${short_name#*'"'} # trim leading double quote

    # Get number of nodes requested
    elif [ "$variable" = "nnodes" ]; then
        nnodes=$value

    # Get list of maps to make
    elif [ "$variable" = "maps" ]; then
        maps=${value%"'"*} # trim trailing single quote
        maps=${maps#*"'"}  # trim leading single quote
        maps=${maps%'"'*}  # trim trailing double quote
        maps=${maps#*'"'}  # trim leading double quote

    # Get map-making time limit
    elif [ "$variable" = "tlimit_map" ]; then
        tlimit_map=${value%"'"*}       # trim trailing single quote
        tlimit_map=${tlimit_map#*"'"}  # trim leading single quote
        tlimit_map=${tlimit_map%"\""*} # trim trailing double quote
        tlimit_map=${tlimit_map#*"\""} # trim leading double quote
        h=${tlimit_map%":"*};h=${h%":"*}
        m=${tlimit_map%":"*};m=${m#*":"}
        s=${tlimit_map#*":"};s=${s#*":"}
        t=$(( 3600*$h + 60*$m + $s ))
    fi
done < ../param/param.params

# Check job node and time limits
if ((nnodes>1000)) || ((t>86400)); then
    echo "ABORTING: Niagara jobs must use 1000 or fewer nodes and take no more than 24 hours."
    exit
fi 

# Loop over map types
i=0 ; maptypes=("tsz" "ksz" "tau" "kap" "tco" "thi")
while [ $i -le 5 ]; do

    # Check if parameter file calls for map type i
    if [[ $maps == *"${maptypes[$i]}"* ]]; then
        #script="make_${maptypes[$i]}_map_${short_name}_${seed}.sh"
        script="make_${maptypes[$i]}_map_${short_name}.sh"

        # Submit job script to Niagara scheduler
        if [ "$2" = '' ]; then
            echo "Submitting WebSky ${maptypes[$i]} map-making run ${script}."
            sbatch $script

        # Make job scripts but don't submit them
        elif [ "$2" = "hold" ]; then
            echo "${maptypes[$i]} run ready to submit: `pwd`//${script}"

        # Prepare debug job script
        elif [ "$2" = "debug" ]; then
            if ((nnodes>4)) || ((t>3600)); then
                echo "ABORTING: debug jobs must use 4 nodes or fewer and take no more than 1 hour."
            else
                sed -i '1 a\#SBATCH -p debug' $script
                echo "${maptypes[$i]} debug run ready to submit: `pwd`/${script}"
            fi

        # Re-compile pks2map if the recompile flag is passed
        elif [ "$2" = "recompile" ]; then
            echo "Recompiling WebSky pks2map.f90."
            cd ..
            pwd
            ml python
            python3 $SCRATCH/peakpatch/python/pks2map.py .
            cd maps
            echo "Submitting WebSky ${maptypes[$i]} map-making run ${script}."
            #sbatch $script
        fi
    fi

    i=$(( $i + 1 ))
done
