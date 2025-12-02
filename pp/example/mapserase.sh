#!/bin/sh

# To run:
#
#     cd peak-patch
#     ./example/mapserase.sh <run-dir>
#
# Where the run you want to erase everyting in directory `<run-dir>/maps`
# except the three shell scripts made when peak-patch.py was run to
# configure the directory.

# Use extended globbing to allow patterning with except (!) and match (*)
shopt -s extglob

if [ -d "$1/maps" ]; then # Confirm the directory <run-dir>/maps exists

        cd $1/maps       # Change directory to <run-dir>/maps/
        rm -r -f !(*.sh) # Erase everything but shell scripts

fi
