#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=1:00:00
#SBATCH --job-name histograms
#SBATCH --output=histograms_%j.txt
#SBATCH --mail-type=END

# Load Niagara modules
module load NiaEnv/2019b intelpython3/2019u5

# Go to SLURM submit directory
cd $SLURM_SUBMIT_DIR

python3 $SCRATCH/peak-patch/python/plotting/run_checks/bin_catalogue_Rsmooth.py $1 # $SCRATCH/runs/s4k_n236_nb20_nt10 
