#!/bin/bash
DEBUG
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=TIME_LIMIT
#SBATCH --job-name JOB_NAME
#SBATCH --output BASH_OUTPUT

# Load virtual environemnt so we can use astropy
module load NiaEnv/2019b
module load python/3.8.5
# module load gcc/9.2.0

# Change directory to the Niagara submit directory
cd $SLURM_SUBMIT_DIR

# Command line arguments for python script
pksc_file=PKSC_FILE
npz_file=NPZ_FILE
z_chi_tab=Z_CHI_TAB
mass_cutoff=MASS_CUTOFF
zmin=ZMIN
zmax=ZMAX
Omega_m=OMEGAM
Omega_b=OMEGAB
Omega_Lambda=OMEGAL
h=LITTLEH
ns=SPECTRALINDEX
sigma8=SIGMA8

# Change directores to the Peak Patch run output directory
catalogue_dir=$(dirname "$(realpath "$pksc_file")")
cd $catalogue_dir

# Run the pksc2npz python script
python3 PP_DIR/python/scripts/pksc2npz.py $pksc_file $npz_file $z_chi_tab $mass_cutoff $zmin $zmax \
    $Omega_m $Omega_b $Omega_Lambda $h $ns $sigma8
