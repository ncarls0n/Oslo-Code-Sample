#!/bin/bash

INPUT_LOC="/scratch/m/murray/vasissua/PeakPatch/music/dm_density.bin"
OUTPUT_LOC="outputs/camb_music_noBaryons.csv"
#INPUT_LOC="/scratch/m/murray/vasissua/PeakPatch/pp_runs/run7/fields/Fvec_256Mpc_n256_nb40_nt1"
#OUTPUT_LOC="outputs/camb_pp.csv"
ntiles=1
nbuff=20
boxsize=4096

mkdir outputs
./power_spectrum $INPUT_LOC  $ntiles $nbuff $boxsize $OUTPUT_LOC
