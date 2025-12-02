import numpy as np, sys, os
from peakpatchtools import PeakPatch

# Peak Patch run directory
run_dir = sys.argv[1]

# Field type to measure power spectrum from, allowed types: 'rhog', 'rho', 'zetag', 'zeta' and 'chi'
field_type = sys.argv[2]

# Load Peak Patch run
r = PeakPatch(run_dir)

# Run fortran script to measure power spectrum
r.get_power_spectrum( field_type=field_type )
