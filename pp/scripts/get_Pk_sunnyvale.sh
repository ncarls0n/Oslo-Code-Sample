#!/bin/bash -l
#PBS -l nodes=1:ppn=128
#PBS -l walltime=6:00:00
#PBS -r n
#PBS -j oe
#PBS -q starq

# go to your working directory
cd $PBS_O_WORKDIR

# Intel compiler and MPI 
module load python/3.10.2

# Peak Patch dir
peakpatch_dir=/fs/lustre/project/act/njcarlson/peakpatch

# Python runner
runner=$peakpatch_dir/python/scripts/get_Pk_sunnyvale.py

# Peak Patch run directory
run_dir=/fs/lustre/project/act/njcarlson/LIMPING/24.10.21_LIMPING_CCAT_EoR-Spec/ng10_cenz7500.0_seed13579_mchi35.925/

# Which field to measure power spectrum from, allowed fields: rhog, rho, zetag, zeta, chi
field_type=chi

# Run script
cd $run_dir
python $runner $run_dir $field_type
