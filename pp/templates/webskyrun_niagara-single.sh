#!/bin/bash
#SBATCH --nodes=NNODESMAP_REPLACE
#SBATCH --ntasks=NTASKSMAP_REPLACE
#SBATCH --ntasks-per-node=TPNODEMAP_REPLACE
#SBATCH --cpus-per-task=NCPUSMAP_REPLACE
#SBATCH --time=TLIMITMAP_REPLACE
#SBATCH --job-name PROFILE_REPLACE
#SBATCH --output=mapmaking_%j.txt
#SBATCH --mail-type=END

###########################################################################
# This is the submission script for serial Websky map-making jobs on the  #
# Niagara supercomputer at SciNet. 
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

# Load Niagara modules needed to make Webskies
module load NiaEnv/2019b intel/2019u4 intelmpi/2019u4 fftw/3.3.8 gsl/2.5 \
    cfitsio/3.450 intelpython2/2019u

# Go to SLURM submit directory
cd $SLURM_SUBMIT_DIR

dir_bin=../bin
ppfile=../output/*merge.pksc*
tabfile=../tables/TABFILEMAP_REPLACE
profile=PROFILE_REPLACE
zmin=ZMINMAP_REPLACE
nside=NSIDEMAP_REPLACE
scramble=SCRAMBLEMAP_REPLACE
center=CENTERMAP_REPLACE
npix=NPIXMAP_REPLACE
fov=FOVMAP_REPLACE
zmax=ZMAXMAP_REPLACE
chihview=CHIHVIEWMAP_REPLACE
model=MODELMAP_REPLACE
PSZcut=PSZCUTMAP_REPLACE
seed=SEEDMAP_REPLACE
base=OUTFILEMAP_REPLACE_${profile}_${seed}
stdout=$base.stdout
stderr=$base.stderr

if [ ! -d $base ] ; then mkdir $base; fi

mpirun -np 1 ${dir_bin}/pks2map $ppfile $base/hp $tabfile \
    $profile $zmin $nside $scramble $center $npix $fov $zmax $chihview \
    $model $PSZcut 2>>$stderr 1>>$stdout

# Clean up map directory
cd $base
mv hp_hp.fits $base\_hp.fits
mv hp_fs.map  $base\_fs.map
mv hp_dp.bin  $base\_dp.bin
mv * ../
cd ../
rm -r $base/
