#!/bin/bash
#SBATCH --nodes=NNODES_REPLACE
#SBATCH --ntasks=NTASKS_REPLACE
#SBATCH --ntasks-per-node=TPNODE_REPLACE
#SBATCH --cpus-per-task=OMPNT_REPLACE
#SBATCH --time=TLIMIT_REPLACE
#SBATCH --job-name SNAME_REPLACE
#SBATCH --output=SNAME_REPLACE_%j.txt
cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

ntasks=NTASKS_REPLACE
ntasksmerge=NTASKSMERGE_REPLACE
lname=LNAME_REPLACE

seed=SEED_REPLACE
rm output/MNAME_REPLACE_merge.pksc.$seed

stdout=logfiles/tabmake.stdout
stderr=logfiles/tabmake.stderr
mpirun -np $ntasks --bynode ./bin/hpkvd 1 $seed 2> $stderr 1> $stdout

while : 
do

    stdout=logfiles/$lname\_$seed.stdout
    stderr=logfiles/$lname\_$seed.stderr
    
    mpirun -np $ntasks --map-by node ./bin/hpkvd 0 $seed 2> $stderr 1> $stdout
    mpirun -np $ntasksmerge --map-by node ./bin/merge_pkvd $seed 2>> $stderr 1>> $stdout
    
    mv Fvec_$lname $seed.delta

    ((seed+=2))

done
