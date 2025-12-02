#!/bin/bash

# USAGE / help string to be output of promted
helpstring="\n
 Run this script AFTER RUNNING ssh-remote.sh ON REMOTE SERVER.                                  \n
                                                                                                \n
 Once you've run ssh-remote.sh on your remote machine, run this script to create an ssh tunnel  \n
 to your remote Jupyter notebooks:                                                              \n
                                                                                                \n\t 
     cd <...>/peakpatch/python/notebooks                                                        \n\t
     ./ssh-local.sh -u username -m remote -d domain -l localhost -r remotehost                  \n
                                                                                                \n
 I've also made a convenience feature so that I can quickly ssh into CITA machines for myself.  \n
 If just one argument is passed at command line, it is assumed to be a cita machine \`remote\`, \n
 and all other values are used as defaults.                                                     \n\t
     ./ssh-local.sh remote                                                                      \n
                                                                                                \n
 If the port 'localhost' is in use, a connection will not be made, and an error message will    \n
 show that explains how to check and kill processes running at a given port. You can also use   \n
 flag '-f' to force any process at port 'localhost' to be terminated, but you should obviously  \n
 be careful about not killing things that are important, so probably just don't use this unless \n
 you definitely know what you're doing. Note that -f must be the first flag for this to work.   \n
                                                                                                \n
 For help, you can run                                                                          \n\t
     ./ssh-remote.sh -h                                                                         \n
 or                                                                                             \n\t
     ./ssh-remote.sh --help                                                                     \n"

# Default parameters
user=njcarlson
remote=mussel
domain=cita.utoronto.ca
localhost=8888
remotehost=8898

# Add default parameters to helpstring
helpstring+="\nDefault parameters:\n\tuser=${user}\n\tremote=${remote}\n\tdomain=${domain}\n
\tlocalhost=${localhost}\n\tremotehost=${remotehost}\n" 

# If odd number of arguments passed at command line
force=false
if [ $(($#%2)) -eq 1 ]; then

    # If prompted for help, print the help string to screen and exit
    if [ "$1" == "-h" -o "$1" == "--help" ]; then
        echo -e $helpstring
	exit

    # If flag -f is passed, then set force to true to kill any process running at port localhost
    elif [ "$1" == "-f" -o "$1" == "--force" ]; then
        force=true

        # Shift command line arguments one to the left
	shift

    else
        # Assume first command line argument is remote if odd number passed
        remote=$1

        # Shift command line arguments one to the left
        shift
    fi
fi

# Loop over even number of remaining command line arguments
while [ $(($#)) -gt 0 ]; do

    # Read command line arguments
    case "$1" in
        "-u") user=$2       ;; # username
        "-m") remote=$2     ;; # remote machine
        "-d") domain=$2     ;; # remote machine's domain
        "-l") localhost=$2  ;; # local host number
        "-r") remotehost=$2 ;; # remote host number
        *)
            echo "Error: flag $1 not understood."
            exit
            ;;
    esac

    # Shift two arguments to the left to read the next flag and its value
    shift
    shift
done

# If the remote server is Niagara, suggest they use SciNet's dedicated JupyterHub
if [ "$remote" == "niagara" ]; then
    echo "You'll probably find it simpler to use SciNet's JupyterHub:"
    echo "https://jupyter.scinet.utoronto.ca/"
fi

# Check that the requested port localhost is available
if lsof -i :$localhost &> /dev/null; then
    
    # Get ID of process running at port localhost	
    PID=$(sudo lsof -i :${localhost} | awk 'NR>1 && !seen[$2]++ {print $2; exit}')

    # If flag -f passed, the process at port localhost will be killed automatically
    if [ "${force}" == true ]; then
        kill -9 $PID

    # Otherwise just display a message telling the user how to kill that process
    else
        echo "The port you've selected ${localhost} is in use."
        echo "To check if a port is in use e.g.:"
        echo "    lsof -i :${localhost}"
        echo "To terminate the process (PID=${PID}) running at port 8888 run:"
        echo "    kill -9 ${PID}"
        echo "Aborting."
        exit
    fi
fi

# Connect ssh tunnel for Jupyter Notebooks
ssh -N -f -L localhost:${localhost}:localhost:${remotehost} ${user}@${remote}.${domain}

# Print the address to connect to Jupyter notebooks in a web browser
echo "You're all set, paste the following into a web browser:"
echo
echo localhost:${localhost}
echo
