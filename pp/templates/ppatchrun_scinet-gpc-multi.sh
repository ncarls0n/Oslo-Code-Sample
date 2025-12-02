#!/bin/bash
# MOAB/Torque submission script for multiple serial jobs on
# SciNet GPC
#
#PBS -l nodes=NNODES_REPLACE:ppn=TPNODE_REPLACE,walltime=TLIMIT_REPLACE
#PBS -N SNAME_REPLACE

cd $PBS_O_WORKDIR

ntasks=NTASKS_REPLACE
ntasksmerge=NTASKSMERGE_REPLACE
lname=LNAME_REPLACE

seed=SEED_REPLACE
rm output/MNAME_REPLACE_merge.pksc.$seed

export OMP_NUM_THREADS=OMPNT_REPLACE

stdout=logfiles/tabmake.stdout
stderr=logfiles/tabmake.stderr
mpirun -np $ntasks --bynode ./bin/hpkvd 1 $seed 2> $stderr 1> $stdout

while : 
do

    stdout=logfiles/$lname\_$seed.stdout
    stderr=logfiles/$lname\_$seed.stderr
    
    mpirun -np $ntasks --bynode ./bin/hpkvd 0 $seed 2> $stderr 1> $stdout
    mpirun -np $ntasksmerge --bynode ./bin/merge_pkvd $seed 2>> $stderr 1>> $stdout
    
    mv Fvec_$lname $seed.delta

    ((seed+=2))

done
