#!/bin/bash
# MOAB/Torque submission script for multiple serial jobs on
# SciNet GPC
#
#PBS -l nodes=1:m32g:ppn=8,walltime=1:00:00
#PBS -N halopyno

module purge
module load gcc intel/13.1.1 python/2.7.5                                            
module load openmpi/intel/1.6.4  
cd $PBS_O_WORKDIR

ntasks=1

stdout=logfiles/trial.stdout
stderr=logfiles/trial.stderr
    
mpirun -np $ntasks ./hod2map.py 2> $stderr 1> $stdout
