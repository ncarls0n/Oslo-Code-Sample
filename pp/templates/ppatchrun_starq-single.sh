#!/bin/bash 
#PBS -l nodes=NNODES_REPLACE:ppn=128
#PBS -l walltime=TLIMIT_REPLACE
#PBS -r n
#PBS -j oe
#PBS -q starq

# Load starq modules
module load openmpi/4.1.6-gcc-ucx fftw/3.3.10-openmpi-ucx gsl/2.7.1 cfitsio/4.0.0 python/3.10.2

# Go to SLURM submit directory
cd $PBS_O_WORKDIR

# For naming of output
lname=LNAME_REPLACE
seed=SEED_REPLACE

# Remove existing merged catalogues from output
old_catalogue="output/MNAME_REPLACE_merge.pksc.$seed"
if [ -f "$old_catalogue" ]; then
    rm -f $old_catalogue
fi

# Set OpenMP threading
export OMP_NUM_THREADS=2

# Set error and output logfile names
stdout=logfiles/$lname\_$seed.stdout
stderr=logfiles/$lname\_$seed.stderr

COMPILATION_REPLACE

./bin/hpkvd 1

# MPI run of hierarchical peak/void finding script hpkvd.f90
mpirun -np NCPUS_REPLACE ./bin/hpkvd 0 $seed 2> $stderr 1> $stdout

# MPI run of merging & exclusion script merge_pkvd.f90
mpirun -np NCPUS_REPLACE ./bin/merge_pkvd $seed 2>> $stderr 1>> $stdout
