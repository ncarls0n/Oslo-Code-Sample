#!/bin/sh

# This script cleans up everything in a Peak Patch run directory <run-dir>
# except its parameter file <run-dir>/param/param.params.
# 
# To run:
#
#     cd peakpatch
#     ./cleanup.sh <run-dir>
#
# Where the run you want to erase everyting in directory `<run-dir>` except
# `<run-dir>/param/param.params`.
# 
# Often one only wants to erase certain directories. Various subsets of the
# Peak Patch run directory can be cleaned up using an additional command
# line argument `<arg>`
# 
#     ./cleanup.sh <run-dir> <arg>
# 
# where the optional command line argument `<arg>` can be:
# 
# <arg> = archive
#     Remove everythin in <run-dir> except the parameter file and run data,
#     move run data to a new subdirecotry archive/<time & date of run>
# 
# <arg> = maps
#     To return the run to the state after finishing the Peak Patch
#     calculation but before the WebSky calculation (i.e. erase everything
#     in the maps subdirectory except for shell scripts).
# 
# <arg> = fail
#     For runs that failed, it is often useful to save space by deleting
#     results, but keeping the log files so you can see what went wrong.

# Use extended globbing to allow patterning with except (!) and match (*)
shopt -s extglob

# Confirm the directory <run-dir> exists
if [ -d "$1" ]; then

    # Change directory to <run-dir>
    cd $1

    # Check that <run-dir> is formatted as a Peak Patch run
    if [ -d "param" ] && [ -f "param/param.params" ]; then

        # If no additional arguments, delete everything in <run-dir>
        if [ "$2" = '' ]; then

            rm -rf !("param")        # Erase everything except params
            cd param                 # Change directory to <run-dir>/params
            # Erase everything except parameter files
            find . -maxdepth 1 ! -name 'param.params' ! -name '*.ini' ! -name '*.conf' ! -name . -print0 | \
					xargs -0 rm -rf

        # If second command-line argument is "archive"
        elif [ "$2" = 'archive' ]; then

            # Erase everything in <run-dir> except params and run data
            rm -rf !("param"|"output"|"maps"|"logfiles"|"archive")

            # Move old run data to archive/<time & date of run>
            olddir=`date -r output "%H:%M:%S_+%m-%d-%Y"`
            mkdir -p archive/$olddir
            mv !("param") archive/$olddir/.

            cd param                 # Change directory to <run-dir>/params
            rm -rf !("param.params") # Erase everything except param.params

        # Return to state after finishing run with setup.sh
        elif [ "$2" = 'maps' ]; then

            # Check that maps directory exists
            if [ -d "maps" ]; then
                cd maps
                rm -rf !(*.sh) # Erase all but shell scripts
            else
                echo "Directory $1/maps/ does not exist"
            fi

        # Clear output, sourcefiles, etc. for failed runs
        elif [ "$2" = 'fail' ]; then
            rm -rf !("param"|"logfiles")

        else

            echo "Command line argument $2 not recognized, aborting."

        fi

    else

        echo "Direcotry $1 not formatted as a Peak Patch run."
        echo "$1/param/param.params must be present."

    fi

else

    echo "Directory $1 does not exist."

fi
