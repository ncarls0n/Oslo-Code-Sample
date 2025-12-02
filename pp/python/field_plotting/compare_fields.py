import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy
from scipy import signal
import pyfftw
import sys

"""
How to make me go:

python3 python/run_checks/check_fields.py runs/<run_1_name> runs/<run_2_name>

"""
if len(sys.argv) != 3:
    # Error message displayed if you pass the wrong number of arguments.
    print('You\'ve passed me the wrong number of arguments at the command '
        +'line.\nCommand line should take the form\ncd <peak_patch_directo'
        +'ry>\npython3.8 python/run_checks/check_fields.py runs/ns128 runs'
        +'/ns128_fNL=0')
    sys.exit(2)

# Reading from command line prompts
run_dir0 = str(sys.argv[1])+'/' # The directory of the L run
run_dir  = str(sys.argv[2])+'/' # The directory of the NL run

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
    with open(run_dir0+'param/param.params') as f:
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


( seed0,nbuff0,Omx0,h0,boxsize0,ntile0,OmB0,ns0,nmesh0,sigma80,Omvac0,
  run_name0,short_name0,maximum_redshift0,global_redshift0,NonGauss0,
  fNL0 ) = execute_parmeter_file(run_dir0+'param/param.params')

( seed,nbuff,Omx,h,boxsize,ntile,OmB,ns,nmesh,sigma8,Omvac,run_name,
  short_name,maximum_redshift,global_redshift,NonGauss,fNL 
  ) = execute_parmeter_file(run_dir+'param/param.params')

# The two fields you're comparing have to have the same spatial and lattice
#  dimensions
if boxsize!=boxsize0 or nmesh!=nmesh0 or nbuff!=nbuff0 or ntile!=ntile0:
    print("I can't compare two runs if they don't have the same spatial an"
        +"d lattice\ndimensions")
    sys.exit(2)
else:
    del(boxsize0,nmesh0,nbuff0,ntile0)

# Computes the LambdaCDM scale factor a(t)
scale_factor  = (1 + maximum_redshift  - global_redshift )**-1
scale_factor0 = (1 + maximum_redshift0 - global_redshift0)**-1

# Additional parameters
nlattice = int( (nmesh-2*nbuff) * ntile + 2*nbuff )
neff     = int( nlattice - 2*nbuff )
Omm      = Omx+OmB        # Omega_m, fraction of energy that can cluster
rhocrit  = 2.775e11*h**2  # critical energy density 3H^2/8piG [Msol Mpc^-3]
rho      = rhocrit*Omm    # average CDM density [Msol Mpc^-3]
deltavir = 200            # <something defining virial collapse I think>
outnum   = 33             # number of columns in <merged_peak_file>
G        = 4.517e-48      # gravitational constant [Mpc^3 Msol^-1 s^-2]
# Note here that we use mass units of solar masses Msol, spatial units of
# megaparsecs Mpc, and time units of seconds s.
Omm0      = Omx0+OmB0
rhocrit0  = 2.775e11*h0**2
rho0      = rhocrit0*Omm0
deltavir0 = 200
outnum0   = 33
G0        = 4.517e-48

filein  = 'Fvec_'+str(int(boxsize))+'Mpc_n'+str(nlattice)+'_nb'+str(nbuff)+'_nt'+str(ntile)
deltafile  = run_dir +'fields/'+filein
deltafile0 = run_dir0+'fields/'+filein
in_delta  = open( deltafile,  'rb' )
in_delta0 = open( deltafile0, 'rb' )

edges = np.linspace( -boxsize/2 , boxsize/2 , neff+1 )
X,Y,Z = np.meshgrid(edges,edges,edges,indexing='ij')

delta0 = np.fromfile(in_delta0,dtype=np.float32,count=-1)
delta0 = np.reshape(delta0, (nlattice,nlattice,nlattice), order='F')
delta0 = delta0[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]

delta = np.fromfile(in_delta,dtype=np.float32,count=-1)
delta = np.reshape(delta, (nlattice,nlattice,nlattice), order='F')
delta = delta[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]

ddelta=delta-delta0

# Field labels
if NonGauss==0:
    field_label='\\delta_L'
else:
    field_label='\\delta_{NL}'

if NonGauss0==0:
    field_label0='\\delta_L'
else:
    field_label0='\\delta_{NL}'

# Create figure
fig1,axs1 = plt.subplots(nrows=2, ncols=2)

# Plot colormesh slice in the x-y plane of the delta field
ax1_00 = axs1[0,0].pcolormesh( X[:,:,0], Y[:,:,0], ddelta[:,:,int(neff/2)], cmap='viridis')
ax1_00_pos = make_axes_locatable(axs1[0,0])
ax1_00_colorbar_pos = ax1_00_pos.append_axes('right', size='5%', pad=0.05)
fig1.colorbar(ax1_00, cax=ax1_00_colorbar_pos)
axs1[0,0].set_aspect(1)
axs1[0,0].set_title(r'$'+field_label+r'(x,y,z=0) - '+field_label0+r'(x,y,z=0)$')
axs1[0,0].set_xlabel(r'$x$ [$h^{-1}$ Mpc]')
axs1[0,0].set_ylabel(r'$y$ [$h^{-1}$ Mpc]')

# Plot colormesh slice in the x-z plane of the delta field
ax1_10 = axs1[1,0].pcolormesh( Z[:,0,:], X[:,0,:], ddelta[:,int(neff/2),:], cmap='viridis')
ax1_10_pos = make_axes_locatable(axs1[1,0])
ax1_10_colorbar_pos = ax1_10_pos.append_axes('right', size='5%', pad=0.05)
fig1.colorbar(ax1_10, cax=ax1_10_colorbar_pos)
axs1[1,0].set_aspect(1)
axs1[1,0].set_title(r'$'+field_label+r'(x,y=0,z) - '+field_label0+r'(x,y=0,z)$')
axs1[1,0].set_xlabel(r'$z$ [$h^{-1}$ Mpc]')
axs1[1,0].set_ylabel(r'$x$ [$h^{-1}$ Mpc]')

# Plot colormesh slice in the y-z plane of the delta field
ax1_11 = axs1[1,1].pcolormesh( Y[0,:,:], Z[0,:,:], ddelta[int(neff/2),:,:], cmap='viridis')
ax1_11_pos = make_axes_locatable(axs1[1,1])
ax1_11_colorbar_pos = ax1_11_pos.append_axes('right', size='5%', pad=0.05)
fig1.colorbar(ax1_11, cax=ax1_11_colorbar_pos)
axs1[1,1].set_aspect(1)
axs1[1,1].set_title(r'$'+field_label+r'(x=0,y,z) - '+field_label0+r'(x=0,y,z)$')
axs1[1,1].set_xlabel(r'$y$ [$h^{-1}$ Mpc]')
axs1[1,1].set_ylabel(r'$z$ [$h^{-1}$ Mpc]')


# Plot a histogram of all delta
max_ddelta = np.max(ddelta)
min_ddelta = np.min(ddelta)
if max_ddelta != min_ddelta:
    num_bins  = int( (max_ddelta-min_ddelta)/10           # number of bins using
        * len(ddelta) / (2*scipy.stats.iqr(ddelta)) +.5)  # Freedman-Draconis rule

    num_bins = 100

    bin_edges = np.linspace(min_ddelta, max_ddelta, num_bins)
else: # keeps the code from crashing if you pass it the same run twice
    num_bins = 10
    bin_edges = np.linspace(-1e-6,1e-6,num_bins)


# Histogram weights
counts,junk  = np.histogram(ddelta, bins=bin_edges)

axs1[0,1].hist(  bin_edges[:-1], bin_edges, weights=counts )
axs1[0,1].set_xlabel(r'$\delta$')
axs1[0,1].set_ylabel(r'counts')

plt.show()
