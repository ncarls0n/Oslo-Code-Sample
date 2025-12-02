#!/usr/bin/env python

# System libraries
import sys

# Peak Patch libraries
from peakpatchtools import PeakPatch
import peakpatchtools as pkp

# WebSky pks2map setup scripts
import peak_patch as pp

# Load the Peak Patch run
r = PeakPatch( sys.argv[1] )

# Re-compile the WebSky pks2map Fortran script
pphome   = pkp.peak_patch_dir

# # Legacy input, to use python-created .bin file (1) instead of .ini file (0)
# # .ini is preferred since in that case the file gets passed directly into PeakPatch
# legacy_input = 1

# Define command to change directories to Peak Patch, the Peak Patch source files, and the copy of
# the source files that will be put in the run directory
pp_dirs  = pp.PeakPatchLocs(r,pphome)

# # Set of filenames of input/output data, as a class 'Input/Output Data Filenames'
# io_files = pp.IODataFilenames(r)

# Set names of log files, standard output files, as a class 'Logging Variables'
log_vars = pp.LogVars(r)

# # Make list of map profiles to run
# profiles, makemapsname = pp.make_list_of_map_profiles_to_run(r)

pp.compile_the_Websky_pipeline_pks2map( pp_dirs , log_vars, '1' )
