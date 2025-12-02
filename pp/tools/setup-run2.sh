#!/bin/bash

# USAGE:
# module load intel/2019u4
# bash setup-run.sh <n> <s> [ <cluster> <n_proc> <s_buff> <iLPT> <cat_len> ]
# 
# <n>^3  the simulation resolution (including buffers)
# <s>^3  the Lagrangian volume in Mpc of the simulation cube (excluding
#     buffers).
#
# Optional parameters:
# <cluster> the computer to run on (Peak Patch is optimized for UofT's
#               Niagara supercomputer, so 'Niagara' is default)
# <n_proc>  self-explanatory (Niagara has 2048 nodes each with 40 CPUs)
# <s_buff>  the Lagrangian space thickness of the buffers in Mpc
# <iLPT>    use linear or quadratic Lagrangian perturbation theory
#               (accepts values 1 or 2)
# <cat_len> the number of 32-bit floating point values per halo in the
#              catalogue files (11 generally, or 33 if outputting shear). 

compiler=ifort
executable=setuprun

# Read in Peak Patch run requirements (or set to defaults)
n=$1
s=$2
cluster=${3-"Niagara"} # Default to Niagara Supercomputer
n_proc=${4-40}         # Default to 40 CPUs/node on Niagara
s_buff=${5-36.0}       # Default to 36.0 Mpc
iLPT=${6-2}            # Default to 2nd order LPT
cat_len=${7-11}        # Default to 11 for no shear output

# Link and compile
${compiler} setup-run.f90 -o ${executable}

# Execute
./${executable} << EOF
${n}
${s}
${cluster}
${n_proc}
${s_buff}
${iLPT}
${cat_len}
EOF

# Clean up
rm -rf setup-run.mod setup-run.o
exit
