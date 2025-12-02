#!/bin/bash

# USAGE / help string to be output of promted
helpstring="\n
 Run this script BEFORE RUNNING ssh-local.sh ON YOUR LOCAL MACHINE.                 \n
                                                                                    \n
 Run this on script on CITA machine in the directory that houses your notebooks:    \n
                                                                                    \n
     cd ~/Repositories/peakpatch/python/notebooks                                   \n
     ./ssh-remote.sh                                                                \n
                                                                                    \n
 Or, to specify a remote host port, e.g. 8889 instead of the default value of 8898, \n
                                                                                    \n
     ./ssh-remote.sh 8889                                                           \n
                                                                                    \n
 For help, you can run                                                              \n
                                                                                    \n
     ./ssh-remote.sh -h                                                             \n
                                                                                    \n
 or                                                                                 \n
                                                                                    \n
     ./ssh-remote.sh --help                                                         \n
                                                                                    \n
 It may be useful to run this script in a screen session so that you can disconnect \n
 from the server without killing your notebook (this is particularly useful if you  \n
 have bad internet connection and don't want to lose work!). To do this, you can run\n
 a screen session:                                                                  \n
                                                                                    \n
     cd ~/Repositories/peakpatch/python/notebooks                                   \n
     screen -S session_name                                                         \n
     ./ssh-remote.sh                                                                \n
                                                                                    \n
 To detach from the screen session (and leave it running in the background), type   \n
 [Ctrl]+[a] to tell bash you want to enter a screen command, then type [d]. You can \n
 see a list of all active screen sessions using \`screen -ls\` and re-attach to a     \n
 screen in the list with \`screen -r session_name\`. You can kill an attached screen  \n
 with the \`exit\` command.                                                           \n"

# If no command line arguments, assume default remote host port
if [ "$#" == 0 ]; then
    remotehost=8898

# If one command line argument is encountered:
elif [ "$#" == 1 ]; then
    case "$1" in
        "-h"    ) echo -e $helpstring ; exit ;; # Return helpstring and exit if flag -h
        "--help") echo -e $helpstring ; exit ;; # Return helpstring and exit if flag --help
        *       ) remotehost=$1              ;; # Else assume cmd line arg is port for remote host
    esac

# Return error message if more than one command line argument is passed
else
    echo "Incorrect number of command line arguments, see docs for more details."
    exit
fi

# Load python module
ml python/3.10.2

# Create SSH tunnel to jupyter notebooks
jupyter notebook --no-browser --port=${remotehost}
