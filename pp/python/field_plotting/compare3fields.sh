#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=0:15:01 
#SBATCH --job-name plotting
#SBATCH --output=plotting_%j.txt
#SBATCH --mail-type=END

module load NiaEnv/2019b python/3.8.5 # Note that intelpython3 uses numpy version 10.16 but we need 10.17 to use keyword argument offset in np.fromfile
cd $SLURM_SUBMIT_DIR

stdout=std.out

python3 -u $SCRATCH/peak-patch/python/plotting/run_checks/compare3fields.py $1 $2 $3 $4 $5 $6 > $stdout
