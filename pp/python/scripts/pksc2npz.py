# pksc2npz.py

####################################################################################################
# Imports                                                                                          #
####################################################################################################

import numpy as np
import sys
import os
from peakpatchtools import z_of_r_comoving, load_halos_pksc_to_npz


####################################################################################################
# Read command line arguments                                                                      #
####################################################################################################

# The raw Peak Patch halo catalogue
pksc_file = sys.argv[1]

# The output .npz file
npz_file = sys.argv[2]

# Load the z of chi comoving table
table_file = np.load( sys.argv[3], allow_pickle=True )
z_chi_tab = table_file['z_chi_tab']

# Cutoff filter scale
mass_cutoff = float(sys.argv[4])

# Redshift cutoff
z_range = ( float(sys.argv[5]) , float(sys.argv[6]) )

# Header for the .npz file
cosmo_header = { 'Omega_M' : float(sys.argv[ 7]) ,
                 'Omega_B' : float(sys.argv[ 8]) ,
                 'Omega_L' : float(sys.argv[ 9]) ,
                 'h'       : float(sys.argv[10]) ,
                 'ns'      : float(sys.argv[11]) ,
                 'sigma8'  : float(sys.argv[12])   }


####################################################################################################
# Read in DM halo catalogue from unformatted Peak Patch output                                     #
####################################################################################################

# Physical constants
H_0   = cosmo_header['h']*100
G     = 4.300917270036279e-09 # km^2 Mpc Msol^-1 s^-2
rho_c = 3*H_0**2 / (8*np.pi*G) # Msol Mpc^-3
rho_m = rho_c * cosmo_header['Omega_M'] # Msol Mpc^-3

# Read in Peak Patch dark matter halo catalogue
x,y,z,dx,dy,dz,M,redshift,zform,zcos = load_halos_pksc_to_npz( pksc_file, mass_cutoff, rho_m,
        z_range, z_chi_tab )


####################################################################################################
# Save the catalogue in the .npz format used by LIMLAM Mocker                                      #
####################################################################################################

# Save the halos to a .npz file that can be read in by LIMLAM mocker
np.savez( npz_file, cosmo_header=cosmo_header, x=x, y=y, z=z, vx=dx, vy=dy, vz=dz, M=M,
        zhalo=redshift, zform=zform, zcos=zcos )
