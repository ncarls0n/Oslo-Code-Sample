#!/bin/bash

dir_bin=$1
ppfile=$2
tabfile=$3
profile=$4
zmin=$5
nside=$6
scramble=$7
center=$8
npix=$9
fov=${10}
zmax=${11}
chihview=${12}
model=${13}
PSZcut=${14}
base=${15}

stdout=$base.stdout
stderr=$base.stderr

if [ ! -d $base ] ; then mkdir $base; fi

mpirun -np NCPUSMAP_REPLACE ${dir_bin}/pks2map $ppfile $base/hp $tabfile $profile $zmin $nside $scramble $center $npix $fov $zmax $chihview $model $PSZcut 2>>$stderr 1>>$stdout

# #GET Cls
# cd $base
# ../anafast ../param.anafast 2>>../$stderr 1>>../$stdout
# mv hp_hp.fits $base\_hp.fits
# mv hp_fs.map  $base\_fs.map
# mv hp_dp.bin  $base\_dp.bin
# mv cl.fits    $base\_cl.fits
# mv * ../
# cd ../
# rm -r $base 
