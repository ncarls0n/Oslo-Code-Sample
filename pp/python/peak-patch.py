#!/usr/bin/env python

# System libraries
import os, subprocess, time

# Math and science libraries
import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Peak Patch libraries
from peakpatchtools import PeakPatch
import peakpatchtools as pkp

####################################################################################################
#                                                                                                  #
#                                 PEAK PATCH INITIALIZATION SCRIPT                                 #
#                                                                                                  #
# This is the main script for seting up Peak Patch and Websky runs. All parameters for a given     #
# Peak Patch-Websky run are set in the a parameter file `<run-directory>/param/param.params`. This #
# script reads in those parameters, copies the source files from `peak-patch/src/` to the run      #
# directory, then compiles and links those source code files.                                      #
#                                                                                                  #
# For quickest implementation on SciNet's Niagara supercomputer, execute this peak-patch.py using  #
# the the shell script `peak-patch/setup.sh`. For more information, see the documentation in       #
# `peak-patch/readme.md`.                                                                          #
#                                                                                                  #
# Note that peak-patch.py is meant to be run in the directory of a Peak Patch run.                 #
#                                                                                                  #
# USAGE (if Peak Patch is in direcotory `<...>/peakpatch/`):                                       #
#     cd <run-directory>                                                                           #
#     python3 <...>/peakpatch/python/peak-patch.py ./param/param.params                            #
#                                                                                                  #
# An example parameter file can be found in `peakpatch/example/param/'.                            #
#                                                                                                  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import peak_patch

peak_patch.main()
