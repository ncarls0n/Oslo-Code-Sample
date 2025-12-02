import numpy as np
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
from halotools.mock_observables import tpcf
import sys,os,gc
"""
How to make me go:

cd peak-patch
python3 python/run_checks/halo_autocorrelation.py <...>/<run_name>

"""
if len(sys.argv) != 2:
    # Error message displayed if you pass the wrong number of arguments.
    raise ValueError('You\'ve passed me the wrong number of arguments at the command '
        +'line.\nCommand line should take the form\ncd <peak_patch_directo'
        +'ry>\npython3.8 python/run_checks/check_fields.py runs/ns128 runs'
        +'/ns128_fNL=0')

# Reading from command line prompts
run_dir  = str(sys.argv[1])+'/' # The directory of the NL run

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
    with open(param_file) as f:
        params = [i.strip() for i in f.readlines()]

    for line in params:
        if(line[:4]=='seed'  or line[:7]=='boxsize' or line[:5]=='nmesh' or
           line[:5]=='nbuff' or line[:5]=='ntile'  or line[:6]=='sigma8' or
           line[:3]=='Omx'   or line[:3]=='OmB'     or line[:5]=='Omvac' or
           line[:2]=='h '    or line[:2]=='ns'      or 
           line[:8]=='run_name'          or line[:10]=='short_name'      or   
           line[:16]=='maximum_redshift' or line[:15]=='global_redshift' or
           line[:8]=='NonGauss'          or line[:3]=='fNL'              or
           line[:10]=='filterfile'                                       ):
            exec( line, globals() )
    return( seed,nbuff,Omx,h,boxsize,ntile,OmB,ns,nmesh,sigma8,Omvac,
            run_name,short_name,maximum_redshift,global_redshift,NonGauss,
            fNL )

( seed,nbuff,Omx,h,boxsize,ntile,OmB,ns,nmesh,sigma8,Omvac,run_name,
  short_name,maximum_redshift,global_redshift,NonGauss,fNL 
  ) = execute_parmeter_file(run_dir+'param/param.params')
s_latt = boxsize/(ntile*(nmesh-2*nbuff))

# Additional parameters
Omm      = Omx+OmB        # Omega_m, fraction of energy that can cluster

m_per_Mpc = 3.0856775814913673e22 # m/pc
H_0       = h*(100*1e3/m_per_Mpc) # s^-1
G         = 4.51724e-48           # Mpc^3 Msol^-1 s^-2
rhocrit   = 3*H_0**2/(8*np.pi*G)  # Msol Mpc^-3
rho_m     = rhocrit*Omm           # Msol Mpc^-3
# Note that Hubble's constant defined in terms of the dimensionless form h
# is H_0 = h/(100 km/s/Mpc). For G, following convention for h, we use
# digits up to and including the first two digits of the 1-sigma
# uncertainty, i.e. G = (4.51724 +/- 0.000019)*10^-48 Mpc^3 Msol^-1 s^-2.
# rho_crit is the critical energy density for a universe with Lambda=0 and
# k=0, rho_crit = 3H_0^2 / (8pi G) which defines the density parameter
# Omega_j = rho_j / rho_crit, rho_m is comoving matter energy denisty.

# Read in raw and merged halo catalogues
merge_file = '{0}/output/{1}_merge.pksc.{2}'.format(run_dir,run_name,seed)
in_merge   = open(merge_file,'rb')

##M = 4*np.pi/3 * RTH**3 * Omm * 2.775e11*h**2

def halo_mass(R): # R is comoving radius in Mpc
    return 4/3*np.pi*R**3*rho_m # halo mass in Msol

# Read size of halo catalogue in bytes
size_bytes = os.path.getsize(merge_file)

# Read catalogue header
N_halos    = np.fromfile( in_merge, dtype=np.int32,count=1 )[0]
RTHmax     = np.fromfile( in_merge, dtype=np.float32,count=1 )[0]
zin        = np.fromfile( in_merge, dtype=np.float32,count=1 )[0]
m_halo_max = halo_mass(RTHmax)

# All catalogue parameters are 32-bit (or 4-byte) integers or floats, thus
N_columns = int( (size_bytes-12)/4/N_halos ) # is the number of columns
catalogue = np.fromfile(in_merge,dtype=np.float32,count=N_columns*N_halos)
catalogue = np.reshape(catalogue,(N_halos,N_columns))

# Lagrangian positions (in comoving Mpc/h)
xL,yL,zL = catalogue[:,7]/h, catalogue[:,8]/h, catalogue[:,9]/h
lagrange = np.vstack((xL,yL,zL)).T; del(xL,yL,zL);gc.collect()

# Eulerian positions (in comoving Mpc/h)
x,y,z = catalogue[:,0]/h, catalogue[:,1]/h, catalogue[:,2]/h
euler = np.vstack((x,y,z)).T; del(x,y,z,catalogue);gc.collect()

# Halo cross correlation function
r  = np.logspace( np.log(s_latt/h), np.log(boxsize/h), 100 )
xi = tpcf(euler, r, period=boxsize/h)

# Plot n(m)
fig,ax = plt.subplots(nrows=1,ncols=1)
ax.set_yscale('log')
ax.set_xscale('log')
ax.plot( (r[1:]-r[:-1])/2*h , xi, label=r'Peak-Patch')
ax.set_xlabel(r'$r$ $[$Mpc$]$')
ax.set_ylabel(r'$\xi(r)$')
ax.legend()

fig.set_size_inches(14,6)
fig.savefig(run_dir+'/halo_autocorrelation_function.png',dpi=100)

