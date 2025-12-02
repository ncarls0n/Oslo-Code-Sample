import numpy as np
import sys

"""
How to make me go:

cd peak-patch
module load intelpython3
python3 python/run_checks/catalogue_quickcheck.py <path2run>/<run_name>/output/<boxsize>Mpc_n<nlattice>_nb<nbuff>_nt<ntile>_<raw/merge>.pksc.<seed>

"""
if len(sys.argv) != 2:
    # Error message displayed if you pass the wrong number of arguments.
    print('You\'ve passed me the wrong number of arguments at the command '
        +'line.\n')
    sys.exit(2)

# Reading from command line prompts
catalogue_file = str(sys.argv[1])
in_catalogue   = open( catalogue_file, 'rb' )
Non    = np.fromfile( in_catalogue, dtype=np.int32,count=1 )[0]
RTHmax = np.fromfile( in_catalogue, dtype=np.float32,count=1 )[0]
zin    = np.fromfile( in_catalogue, dtype=np.float32,count=1 )[0]

print()
print(' '+catalogue_file )
print(' Number of halos     Non =', Non    )
print(' Max filter scale RTHmax =', RTHmax )
print(' Initial redshift    zin =', zin    )
print()
