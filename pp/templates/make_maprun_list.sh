#!/bin/bash

dir_bin=../bin
dir_cat=../output

tabfile=../tables/TABFILEMAP_REPLACE
sfrfile=../tables/TABFILESFR_REPLACE
model=MODELMAP_REPLACE
zmin=ZMINMAP_REPLACE
nside=NSIDEMAP_REPLACE
center=CENTERMAP_REPLACE
npix=NPIXMAP_REPLACE
fov=FOVMAP_REPLACE
zmax=ZMAXMAP_REPLACE
chihview=CHIHVIEWMAP_REPLACE
PSZcut=PSZCUTMAP_REPLACE
seed=SEEDMAP_REPLACE

#DELETE OLD FILES
rm -f *stdout *stderr maprun_list

#CREATE MAPTABLE
if [ ! -f $tabfile  ]; then
    ${dir_bin}/make_cmbtable $tabfile $model 
fi

for ppfile in `/bin/ls $dir_cat/*merge*pksc*`
do
    seed=${ppfile: -5}
    for profile in `echo MAPS_REPLACE`
    do 
	for scramble in `echo SCRAMBLEMAP_REPLACE`
	do
	    mapcode=$profile\_$seed
	    if [ $scramble -eq '1' ] ; then mapcode=$mapcode\_shuf; fi

	    base=OUTFILEMAP_REPLACE_$mapcode
	    stdout=$mapcode.stdout
	    stderr=$mapcode.stderr
	    
	    if [ ! -f $ppfile ] ; then exit 1 ; fi

	    echo make_single_map.sh \
		$dir_bin \
		$ppfile \
		$tabfile \
		$profile \
		$zmin \
		$nside \
		$scramble \
		$center \
		$npix \
		$fov \
		$zmax \
		$chihview \
		$model \
		$PSZcut \
		$base \
		>> maprun_list
       	done
    done
    
    ((seed+=2))

done

