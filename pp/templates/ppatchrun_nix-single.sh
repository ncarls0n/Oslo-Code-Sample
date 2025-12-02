#!/usr/bin/env bash 

## Load Niagara modules
#module load NiaEnv/2019b intel/2019u4 intelmpi/2019u4 fftw/3.3.8 gsl/2.5 \
#    cfitsio/3.450 intelpython2/2019u5
#
## Go to SLURM submit directory
#cd $SLURM_SUBMIT_DIR

# For naming of output
lname=LNAME_REPLACE
seed=SEED_REPLACE

# Remove existing merged catalogues from output
old_catalogue="output/MNAME_REPLACE_merge.pksc.$seed"
if [ -f "$old_catalogue" ]; then
    rm -f $old_catalogue
fi

## Set OpenMP threading
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Set error and output logfile names
stdout=logfiles/$lname\_$seed.stdout
stderr=logfiles/$lname\_$seed.stderr

# MPI run of hierarchical peak/void finding script hpkvd.f90
mpirun -np NCPUS_REPLACE ./bin/hpkvd 0 $seed 2> $stderr 1> $stdout

# MPI run of merging & exclusion script merge_pkvd.f90
mpirun -np NCPUS_REPLACE ./bin/merge_pkvd $seed 2>> $stderr 1>> $stdout
