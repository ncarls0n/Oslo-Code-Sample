#!/bin/bash
#SBATCH --nodes=NNODESMAP_REPLACE
#SBATCH --ntasks=NTASKSMAP_REPLACE
#SBATCH --ntasks-per-node=TPNODEMAP_REPLACE
#SBATCH --cpus-per-task=NCPUSMAP_REPLACE
#SBATCH --time=TLIMITMAP_REPLACE
#SBATCH --job-name mapmaking
#SBATCH --output=mapmaking_%j.txt
#SBATCH --mail-type=END

###########################################################################
# Submission script for multiple serial map-making jobs on SciNet Niagara #
#                                                                         #
# Main script for making maps from lightcone catalogs with pks2cmb        #
#                                                                         #
# Three files will be created in the ./maps subdirectory:                 #
#     make_maps_*.sh                                                      #
#     make_maprunlist_*.sh                                                #
#     make_single_maps.h                                                  #
#                                                                         #
# After lightcones are done (assumed to be in ./output) do:               #
#     cd ./maps                                                           #
#     ./make_maprunlist_*.sh                                              #
#     sbatch make_maps_*.sh                                               #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Load Niagara modules
module load NiaEnv/2019b gnu-parallel intel/2019u4 intelmpi/2019u4

# Go to SLURM submit directory
cd $SLURM_SUBMIT_DIR

# Set OpenMP threading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Use GNU Parallel to run a single map MPI run for each Websky profile
# specified in the parameter file for this run
cat maprun_list | parallel --jobs NTASKSMAP_REPLACE --sshloginfile \
    $SLURM_JOB_NODELIST --workdir $SLURM_SUBMIT_DIR
