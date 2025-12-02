# pksc2pksc.py
# 
# This script allows you to trim a Peak Patch halo catalogue file, removing
# all halos outside a given range of positions, or above or below a
# specified mass cutoff.
# 
# USAGE: at command line, run
# 
#     python halo_hist2d.py <catalogue_file> <axis>=0 <xmin>='-inf'
#         <xmax>='inf' <mass_cut>='none' <mass_cutoff>=1e11
#
# where
#  - <catalogue_file> is a Peak-Patch formatted halo catalogue ending in
#        `.pksc.<seed_number>`
#  - <axis> is the axis to truncate (0, 1 or 2, default is 0)
#  - <xmin> is the minimum Eulerian position to keep halos from. If set to
#        '-inf', there is no minimum (default is '-inf')
#  - <xmax> is the maximum Eulerian position to keep halos from. If set to
#        'inf', there is no maximum (default is 'inf')
#  - <mass_cut> determines whether to filter masses ('low' for lowpass
#        filter, 'high for highpass filter, 'none' for no mass cutoff,
#        default is 'none')
#  - <mass_cutoff> is mass in solar masses about which to perform cutoff
#        (default is 1e11)
#
# Alternatively, if <mass_cut> is set to 'Rhigh' or 'Rlow' then we instead
# filter by filter scale, instead of <mass_cutoff> being a mass, you should
# enter the lattice size in Mpc.

import numpy as np
import os,sys

# Read halo catalogue file from command line
cat_file = sys.argv[1]
in_cat   = open( cat_file, 'rb' )

# Read catalogue file header
N_halos = np.fromfile( in_cat, dtype=np.int32  , count=1 )[0]
RTHmax  = np.fromfile( in_cat, dtype=np.float32, count=1 )[0]
zin     = np.fromfile( in_cat, dtype=np.float32, count=1 )[0]

# Determine number of columns in halo catalogue table
size_bytes = os.path.getsize(cat_file)
cols       = int( (size_bytes-12)/4/N_halos )

# Read halo catalogue
catalogue = np.reshape( np.fromfile( in_cat, dtype=np.float32, count=N_halos*cols ), (N_halos,cols) )

# Read instructions from command line
axis,xmin,xmax,mass_cut,mass_cutoff = 0,-np.inf,np.inf,'none',1e11
if len(sys.argv)>2: axis        = int       ( sys.argv[2] )
if len(sys.argv)>3: xmin        = np.float32( sys.argv[3] ) 
if len(sys.argv)>4: xmax        = np.float32( sys.argv[4] )
if len(sys.argv)>5: mass_cut    =             sys.argv[5]
if len(sys.argv)>6: mass_cutoff = float     ( sys.argv[6] )

# Cut halos outside of bounds [xmin,xmax]
chi = catalogue[:,axis].T
if xmin > -np.inf:
    chicut = np.where( chi > xmin )[0]
    chi,catalogue = chi[chicut],catalogue[chicut,:]
if xmax < np.inf:
    chicut = np.where( chi < xmax )[0]
    catalogue = catalogue[chicut,:]
del(chi)

# Halo masses
h       = .6735
Omega_m = .3138
rhocrit = 2.775e11*h**2
rho_m   = Omega_m * rhocrit
Rth     = catalogue[:,6]
M       = 4*np.pi/3*Rth**3*rho_m

# Cut halos
if mass_cut=='low':
    Mcut = np.where( M<mass_cutoff )[0]
    catalogue = catalogue[Mcut,:]
elif mass_cut=='high':
    Mcut = np.where( M>mass_cutoff )[0]
    catalogue = catalogue[Mcut,:]
elif mass_cut=='Rlow':
    Rcut = np.where( Rth<mass_cutoff )[0]
    catalogue = catalogue[Rcut,:]
elif mass_cut=='Rhigh':
    Rcut = np.where( Rth>mass_cutoff )[0]
    catalogue = catalogue[Rcut,:]
del(M,Rth)

# Open output file
fileout = cat_file
if xmin > -np.inf     : fileout += '_xmin'+str(xmin)
if xmax <  np.inf     : fileout += '_xmax'+str(xmax)
if   mass_cut=='low'  : fileout += '_M_lt_'+str(np.format_float_scientific(mass_cutoff))+'Msol'
elif mass_cut=='high' : fileout += '_M_gt_'+str(np.format_float_scientific(mass_cutoff))+'Msol'
elif mass_cut=='Rlow' : fileout += '_Rth_lt_'+str(np.format_float_scientific(mass_cutoff))+'Mpc'
elif mass_cut=='Rhigh': fileout += '_Rth_gt_'+str(np.format_float_scientific(mass_cutoff))+'Mpc'
if os.path.exists(fileout): os.system('rm -f '+fileout)
out_file = open(fileout,'wb')

print('\nwriting to\n'+fileout+'\n')

# Write header to file
N_halos  = len(catalogue[:,0])
catalogue = np.reshape( catalogue, (N_halos*cols) )
out_file.write(
    np.int32  ( N_halos ).tobytes() +
    np.float32( RTHmax  ).tobytes() +
    np.float32( zin     ).tobytes() +
    np.float32(catalogue).tobytes()   )
