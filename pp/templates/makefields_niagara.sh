#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=6:00:00
#SBATCH --job-name plotting
#SBATCH --output=plotting_%j.txt
#SBATCH --mail-type=END

###########################################################################
# This is the submission script for plotting initial conditions fields on #
# the Niagara supercomputer at SciNet. Make fields for a given run in the #
# directory <run_dir>                                                     # 
#                                                                         #
#     cd <run_dir>                                                        #
#     cp <...>/peak-patch/templates/makefields_niagara.sh .               #
#     sbatch makefields_niagara.sh <field_list>                           #
#                                                                         #
# Where <field_list> is "delta", "deltag", "deltang", "zeta", "zetag",    #
# "zetang", "chi" or "all"                                                #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Load Niagara modules needed to make Webskies
module load NiaEnv/2019b python/3.8.5

# Go to SLURM submit directory
cd $SLURM_SUBMIT_DIR

fieldlist=${1:-delta}
rundir=${2:-.}
pkpdir=${3:-$SCRATCH/peak-patch}
stdout=std.out

python3 -u $pkpdir/python/plotting/run_checks/check_fields2.0.py $rundir \
    $fieldlist > $stdout
