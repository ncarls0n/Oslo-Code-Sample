#!/bin/bash
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --job-name pksc2npz
#SBATCH --output ./pksc2npz.out

# Load virtual environemnt so we can use astropy
module load NiaEnv/2019b
module load python/3.8.5 gcc/9.2.0
# source ~/.virtualenvs/myenv/bin/activate
source $SCRATCH/.virtualenvs/astropy_env/bin/activate

export $XDG_CONFIG_HOME=/scratch/r/rbond/nate/
# export ASTROPY_CONFIGDIR=/scratch/r/rbond/nate/.astropy/config
# export ASTROPY_CACHE_DIR=/scratch/r/rbond/nate/.astropy/cache

# Change directory to the Niagara submit directory
cd $SLURM_SUBMIT_DIR

# Command line arguments for python script
run_dir=/scratch/r/rbond/nate/peakpatch/
catalogue_file=
parallel_cores=40

# Run the pksc2npz python script
python3 /scratch/r/rbond/nate/peakpatch/python/scripts/pksc2npz.py $run_dir $catalogue_file $parallel_cores

# Deactivate the virtual environment
deactivate
