#!/bin/bash

#rm -f *.asci *.bin
for f in *merge.pksc*
do 
  python ~/src/peak-patch/python/pksc2asci.py $f foo 
  cat foo.asci >> allhalos.asci
#  cat foo.bin  >> allhalos.bin
done
