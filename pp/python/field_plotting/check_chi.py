import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy
from scipy import signal
import pyfftw
import sys

"""
How to make me go:

python3 python/run_checks/check_fields.py <path2runs>

Where there are directories <path2runs>/fnl0, <path2runs>/fnl1e5,
<path2runs>/fnl1e6, <path2runs>/fnl1e7
"""
if len(sys.argv) != 2:
    # Error message displayed if you pass the wrong number of arguments.
    print('You\'ve passed me the wrong number of arguments at the command '
        +'line.')
    sys.exit(2)

# Reading from command line prompts
runs_dir = str(sys.argv[1])+'/' # The directory of the L run
fnl0_dir   = runs_dir+'/fnl0/'
fnl1e5_dir = runs_dir+'/fnl1e5/'
fnl1e6_dir = runs_dir+'/fnl1e6/'
fnl1e7_dir = runs_dir+'/fnl1e7/'

###########################################################################
### Execute useful lines from parameter file                            ###
###########################################################################

def execute_parmeter_file(param_file):
    # Quite literally just goes to the parameter file (whose location is
    # passed via the string param_file and executes a few select lines of 
    # it. Peak patch parameter files are written in python, so I've simply
    # done a line search and if a line starts with the right variable name,
    # it is simply executed. This makes this script a little harder to
    # read, but I don't feel like deeling with all the fucking spacing to
    # make it into a dictionary or something.
    #
    # The variables defined are:
    # seed               initial random field seed number (an integer)
    # boxsize [h^-1 Mpc] sidelength of simulation volume
    # nmesh              number of cells in one tile (including buffers)
    # nbuff              number of cells in buffer
    # ntile              number
    # sigma8             stdev of Gaussian random field smoothed at
    #                        8 h^-1 Mpc
    # Omx                Omega_CDM, energy density fraction of CDM
    # OmB                Omega_B, energy density fraction of baryonic
    #                        matter
    # Omvac              Omega_Lambda, energy denstiy fraction of DE
    # h                  "little h", H_0/100 km/s/Mpc
    # ns                 n_s, spectral index
    # run_name           string used in output files
    # short_name         shorter string used in output files
    # maximum_redshift   z_max, redshift of primordial fields
    # global_redshift    z_0, redshift of Eulerian peaks
    # NonGauss           Type of NonGaussianity to implement
    # fNL                f_NL, nonlinearity coefficent
    #                        delta=delta + f_NL(delta^2-<delta^2>)

    # Read parameter file
    with open(fnl0_dir+'param/param.params') as f:
        params = [i.strip() for i in f.readlines()]

    for line in params:
        if(line[:4]=='seed'  or line[:7]=='boxsize' or line[:5]=='nmesh' or
           line[:5]=='nbuff' or line[:5]=='ntile'  or line[:6]=='sigma8' or
           line[:3]=='Omx'   or line[:3]=='OmB'     or line[:5]=='Omvac' or
           line[:2]=='h '    or line[:2]=='ns'      or 
           line[:8]=='run_name'          or line[:10]=='short_name'      or   
           line[:16]=='maximum_redshift' or line[:15]=='global_redshift' or
           line[:8]=='NonGauss'          or line[:3]=='fNL'              ):
            exec( line, globals() )
    return( seed,nbuff,Omx,h,boxsize,ntile,OmB,ns,nmesh,sigma8,Omvac,
            run_name,short_name,maximum_redshift,global_redshift,NonGauss,
            fNL )


( seed,nbuff,Omx,h,boxsize,ntile,OmB,ns,nmesh,sigma8,Omvac,run_name,
  short_name,maximum_redshift,global_redshift,NonGauss,fNL 
  ) = execute_parmeter_file(fnl0_dir+'param/param.params')

# Additional parameters
nlattice = int( (nmesh-2*nbuff) * ntile + 2*nbuff )
neff     = int( nlattice - 2*nbuff )
Omm      = Omx+OmB        # Omega_m, fraction of energy that can cluster
rhocrit  = 2.775e11*h**2  # critical energy density 3H^2/8piG [Msol Mpc^-3]
rho      = rhocrit*Omm    # average CDM density [Msol Mpc^-3]
deltavir = 200            # <something defining virial collapse I think>
outnum   = 33             # number of columns in <merged_peak_file>
G        = 4.517e-48      # gravitational constant [Mpc^3 Msol^-1 s^-2]
a_latt   = boxsize/neff
R_keep   = 2*RTHmax#5*a_latt # Radius to take filter peaks from
# Note here that we use mass units of solar masses Msol, spatial units of
# megaparsecs Mpc, and time units of seconds s.


# Open delta files
filein      = 'Fvec_'+str(int(boxsize))+'Mpc_n'+str(nlattice)+'_nb'+str(nbuff)+'_nt'+str(ntile)
delta_fnl0_file   = fnl0_dir  +'fields/'+filein
delta_fnl1e5_file = fnl1e5_dir+'fields/'+filein
delta_fnl1e6_file = fnl1e6_dir+'fields/'+filein
delta_fnl1e7_file = fnl1e7_dir+'fields/'+filein
in_delta_fnl0   = open( delta_fnl0_file,   'rb' )
in_delta_fnl1e5 = open( delta_fnl1e5_file, 'rb' )
in_delta_fnl1e6 = open( delta_fnl1e6_file, 'rb' )
in_delta_fnl1e7 = open( delta_fnl1e7_file, 'rb' )

# Make meshgrid for spatial axes of fields
edges = np.linspace( -boxsize/2 , boxsize/2 , neff+1 )
X,Y,Z = np.meshgrid(edges,edges,edges,indexing='ij')
X_slice = X[:,:,0]
Y_slice = Y[:,:,0]
del(X,Y,Z)

# Open merged peak catalogue files
cataloguein = str(int(boxsize))+'Mpc_n'+str(nlattice)+'_nb'+str(nbuff)+'_nt'+str(ntile)+'_merge.pksc.'+str(seed)
cat_fnl0_file   = fnl0_dir  +'/output/'+cataloguein
cat_fnl1e5_file = fnl1e5_dir+'/output/'+cataloguein
cat_fnl1e6_file = fnl1e6_dir+'/output/'+cataloguein
cat_fnl1e7_file = fnl1e7_dir+'/output/'+cataloguein
in_cat_fnl0   = open( cat_fnl0_file,   'rb' )
in_cat_fnl1e5 = open( cat_fnl1e5_file, 'rb' )
in_cat_fnl1e6 = open( cat_fnl1e6_file, 'rb' )
in_cat_fnl1e7 = open( cat_fnl1e7_file, 'rb' )

# For making circles
theta = np.linspace( 0, 2*np.pi, 100, endpoint=True )

# Create figure
fig1,axs1 = plt.subplots(nrows=2, ncols=2)

# Get slice of f_NL=0 field
delta_fnl0 = np.fromfile(in_delta_fnl0,dtype=np.float32,count=-1)
delta_fnl0 = np.reshape(delta_fnl0, (nlattice,nlattice,nlattice), order='F')
delta_fnl0 = delta_fnl0[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
delta_fnl0_slice = delta_fnl0[:,:,int(neff/2)]
del(delta_fnl0)

# Plot colormesh slice in the x-y plane of the delta field
ax1_00 = axs1[0,0].pcolormesh( X_slice, Y_slice, delta_fnl0_slice, cmap='viridis')
ax1_00_pos = make_axes_locatable(axs1[0,0])
ax1_00_colorbar_pos = ax1_00_pos.append_axes('right', size='5%', pad=0.05)
fig1.colorbar(ax1_00, cax=ax1_00_colorbar_pos)
axs1[0,0].set_aspect(1)
axs1[0,0].set_title(r'$\delta_G(x,y,z=0)$')
axs1[0,0].set_xlabel(r'$x$ [$h^{-1}$ Mpc]')
axs1[0,0].set_ylabel(r'$y$ [$h^{-1}$ Mpc]')
del(delta_fnl0_slice)

# Plot Lagrangian space spherical peak patches
Non       = np.fromfile( in_cat_fnl0, dtype=np.int32,   count=1 )[0]
RTHmax    = np.fromfile( in_cat_fnl0, dtype=np.float32, count=1 )[0]
R_keep    = RTHmax 
zin       = np.fromfile( in_cat_fnl0, dtype=np.float32, count=1 )[0]
catalogue = np.fromfile( in_cat_fnl0, dtype=np.float32, count=33*Non )
catalogue = np.reshape( catalogue, (Non,33) )
catalogue = catalogue[:,6:10] # Only interested in Rth and lagrangian position for this case
catalogue = catalogue[catalogue[:,3].argsort()] # sort based on zL
catalogue = catalogue[
    np.searchsorted(catalogue[:,3],-R_keep,side='left') :
    np.searchsorted(catalogue[:,3],R_keep,side='right') ]
for i in range(len(catalogue[:,0])):
    axs1[0,0].plot( catalogue[i,0]*np.cos(theta)+catalogue[i,1],
                    catalogue[i,0]*np.sin(theta)+catalogue[i,2],
                    color='white', alpha=catalogue[i,0]/RTHmax )
del(catalogue)
axs1[0,0].set_xlim( -boxsize/2 , boxsize/2 )
axs1[0,0].set_ylim( -boxsize/2 , boxsize/2 )

# Get slice of f_NL=1e5 field
delta_fnl1e5 = np.fromfile(in_delta_fnl1e5,dtype=np.float32,count=-1)
delta_fnl1e5 = np.reshape(delta_fnl1e5, (nlattice,nlattice,nlattice), order='F')
delta_fnl1e5 = delta_fnl1e5[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
delta_fnl1e5_slice = delta_fnl1e5[:,:,int(neff/2)]
del(delta_fnl1e5)

# Plot colormesh slice in the x-z plane of the delta field
ax1_01 = axs1[0,1].pcolormesh( X_slice, Y_slice, delta_fnl1e5_slice, cmap='viridis')
ax1_01_pos = make_axes_locatable(axs1[0,1])
ax1_01_colorbar_pos = ax1_01_pos.append_axes('right', size='5%', pad=0.05)
fig1.colorbar(ax1_01, cax=ax1_01_colorbar_pos)
axs1[0,1].set_aspect(1)
axs1[0,1].set_title(r'$\delta_{nG,f_{NL}=10^5}(x,y,z=0)$')
axs1[0,1].set_xlabel(r'$x$ [$h^{-1}$ Mpc]')
axs1[0,1].set_ylabel(r'$y$ [$h^{-1}$ Mpc]')
del(delta_fnl1e5_slice)

# Plot Lagrangian space spherical peak patches
Non       = np.fromfile( in_cat_fnl1e5, dtype=np.int32,   count=1 )[0]
RTHmax    = np.fromfile( in_cat_fnl1e5, dtype=np.float32, count=1 )[0]
zin       = np.fromfile( in_cat_fnl1e5, dtype=np.float32, count=1 )[0]
catalogue = np.fromfile( in_cat_fnl1e5, dtype=np.float32, count=33*Non )
catalogue = np.reshape( catalogue, (Non,33) )
catalogue = catalogue[:,6:10] # Only interested in Rth and lagrangian position for this case
catalogue = catalogue[catalogue[:,3].argsort()] # sort based on zL
catalogue = catalogue[
    np.searchsorted(catalogue[:,3],-R_keep,side='left') :
    np.searchsorted(catalogue[:,3],R_keep,side='right') ]
for i in range(len(catalogue[:,0])):
    axs1[0,1].plot( catalogue[i,0]*np.cos(theta)+catalogue[i,1],
                    catalogue[i,0]*np.sin(theta)+catalogue[i,2],
                    color='white', alpha=catalogue[i,0]/RTHmax )
del(catalogue)
axs1[0,1].set_xlim( -boxsize/2 , boxsize/2 )
axs1[0,1].set_ylim( -boxsize/2 , boxsize/2 )

# Get slice of f_NL=1e6 field
delta_fnl1e6 = np.fromfile(in_delta_fnl1e6,dtype=np.float32,count=-1)
delta_fnl1e6 = np.reshape(delta_fnl1e6, (nlattice,nlattice,nlattice), order='F')
delta_fnl1e6 = delta_fnl1e6[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
delta_fnl1e6_slice = delta_fnl1e6[:,:,int(neff/2)]
del(delta_fnl1e6)

# Plot colormesh slice in the y-z plane of the delta field
ax1_10 = axs1[1,0].pcolormesh( X_slice, Y_slice, delta_fnl1e6_slice, cmap='viridis')
ax1_10_pos = make_axes_locatable(axs1[1,0])
ax1_10_colorbar_pos = ax1_10_pos.append_axes('right', size='5%', pad=0.05)
fig1.colorbar(ax1_10, cax=ax1_10_colorbar_pos)
axs1[1,0].set_aspect(1)
axs1[1,0].set_title(r'$\delta_{nG,f_{NL}=10^6}(x,y,z=0)$')
axs1[1,0].set_xlabel(r'$x$ [$h^{-1}$ Mpc]')
axs1[1,0].set_ylabel(r'$y$ [$h^{-1}$ Mpc]')
del(delta_fnl1e6_slice)

# Plot Lagrangian space spherical peak patches
Non       = np.fromfile( in_cat_fnl1e6, dtype=np.int32,   count=1 )[0]
RTHmax    = np.fromfile( in_cat_fnl1e6, dtype=np.float32, count=1 )[0]
zin       = np.fromfile( in_cat_fnl1e6, dtype=np.float32, count=1 )[0]
catalogue = np.fromfile( in_cat_fnl1e6, dtype=np.float32, count=33*Non )
catalogue = np.reshape( catalogue, (Non,33) )
catalogue = catalogue[:,6:10] # Only interested in Rth and lagrangian position for this case
catalogue = catalogue[catalogue[:,3].argsort()] # sort based on zL
catalogue = catalogue[
    np.searchsorted(catalogue[:,3],-R_keep,side='left') :
    np.searchsorted(catalogue[:,3],R_keep,side='right') ]
for i in range(len(catalogue[:,0])):
    axs1[1,0].plot( catalogue[i,0]*np.cos(theta)+catalogue[i,1],
                    catalogue[i,0]*np.sin(theta)+catalogue[i,2],
                    color='white', alpha=catalogue[i,0]/RTHmax )
del(catalogue)
axs1[1,0].set_xlim( -boxsize/2 , boxsize/2 )
axs1[1,0].set_ylim( -boxsize/2 , boxsize/2 )

# Get slice of f_NL=1e7 field
delta_fnl1e7 = np.fromfile(in_delta_fnl1e7,dtype=np.float32,count=-1)
delta_fnl1e7 = np.reshape(delta_fnl1e7, (nlattice,nlattice,nlattice), order='F')
delta_fnl1e7 = delta_fnl1e7[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
delta_fnl1e7_slice = delta_fnl1e7[:,:,int(neff/2)]
del(delta_fnl1e7)

# Plot colormesh slice in the y-z plane of the delta field
ax1_11 = axs1[1,1].pcolormesh( X_slice, Y_slice, delta_fnl1e7_slice, cmap='viridis')
ax1_11_pos = make_axes_locatable(axs1[1,1])
ax1_11_colorbar_pos = ax1_11_pos.append_axes('right', size='5%', pad=0.05)
fig1.colorbar(ax1_11, cax=ax1_11_colorbar_pos)
axs1[1,1].set_aspect(1)
axs1[1,1].set_title(r'$\delta_{nG,f_{NL}=10^7}(x,y,z=0)$')
axs1[1,1].set_xlabel(r'$x$ [$h^{-1}$ Mpc]')
axs1[1,1].set_ylabel(r'$y$ [$h^{-1}$ Mpc]')
del(delta_fnl1e7_slice)

# Plot Lagrangian space spherical peak patches
Non       = np.fromfile( in_cat_fnl1e7, dtype=np.int32,   count=1 )[0]
RTHmax    = np.fromfile( in_cat_fnl1e7, dtype=np.float32, count=1 )[0]
zin       = np.fromfile( in_cat_fnl1e7, dtype=np.float32, count=1 )[0]
catalogue = np.fromfile( in_cat_fnl1e7, dtype=np.float32, count=33*Non )
catalogue = np.reshape( catalogue, (Non,33) )
catalogue = catalogue[:,6:10] # Only interested in Rth and lagrangian position for this case
catalogue = catalogue[catalogue[:,3].argsort()] # sort based on zL
catalogue = catalogue[
    np.searchsorted(catalogue[:,3],-R_keep,side='left') :
    np.searchsorted(catalogue[:,3],R_keep,side='right') ]
for i in range(len(catalogue[:,0])):
    axs1[1,1].plot( catalogue[i,0]*np.cos(theta)+catalogue[i,1],
                    catalogue[i,0]*np.sin(theta)+catalogue[i,2],
                    color='white', alpha=catalogue[i,0]/RTHmax )
del(catalogue)
axs1[1,1].set_xlim( -boxsize/2 , boxsize/2 )
axs1[1,1].set_ylim( -boxsize/2 , boxsize/2 )





fig1.set_size_inches(14,12)
fig1.savefig(runs_dir+'/delta.png',dpi=100)

# plt.show()
