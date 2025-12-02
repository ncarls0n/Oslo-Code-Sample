#!/bin/bash
# MOAB/Torque submission script for multiple serial jobs on
# SciNet GPC
#
#PBS -l nodes=NNODESMAP_REPLACE:ppn=8,walltime=TLIMITMAP_REPLACE
#PBS -N mapmaking

######################################################################
#
# Main script for making maps from lightcone catalogs with pks2cmb
#
# Three files will be created in the ./maps subdirectory:
#     make_maps_*.sh
#     make_maprunlist_*.sh
#     make_single_maps.h
#
# After lightcones are done (assumed to be in ./output) do:
#     cd ./maps
#     ./make_maprunlist_*.sh
#     qsub make_maps_*.sh
#    
######################################################################

cd $PBS_O_WORKDIR

module load gnu-parallel

cat maprun_list | parallel -j PPNMAP_REPLACE --sshloginfile $PBS_NODEFILE --workdir $PWD






