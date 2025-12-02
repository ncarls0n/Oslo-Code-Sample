import numpy as np
import matplotlib.pyplot as plt 
import sys
import os
"""
How to make me go:

cd peak-patch
python3 python/run_checks/bin_catalogue_Rsmooth.py <path2runs>/<run_1_name>

"""
if len(sys.argv) != 2:
    # Error message displayed if you pass the wrong number of arguments.
    print('You\'ve passed me the wrong number of arguments at the command '
        +'line.\nCommand line should take the form\ncd <peak_patch_directo'
        +'ry>\npython3.8 python/run_checks/check_fields.py runs/ns128 runs'
        +'/ns128_fNL=0')
    sys.exit(2)

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

# Additional parameters
Omm      = Omx+OmB        # Omega_m, fraction of energy that can cluster
rhocrit  = 2.775e11*h**2  # critical energy density 3H^2/8piG [Msol Mpc^-3]
rho      = rhocrit*Omm    # average CDM density [Msol Mpc^-3]
deltavir = 200            # <something defining virial collapse I think>
G        = 4.517e-48      # gravitational constant [Mpc^3 Msol^-1 s^-2]
# Note here that we use mass units of solar masses Msol, spatial units of
# megaparsecs Mpc, and time units of seconds s.

# Read in raw and merged halo catalogues
merge_file = '{0}/output/{1}_merge.pksc.{2}'.format(run_dir,run_name,seed)
in_merge   = open(merge_file,'rb')

#in_raw   = open( '{0}/output/{1}_raw.pksc.{2}'.format(
#    run_dir,run_name,seed),'rb')
#
## Raw catalogue
#Non_r     = np.fromfile( in_raw, dtype=np.int32,count=1 )[0]
#RTHmax_r  = np.fromfile( in_raw, dtype=np.float32,count=1 )[0]
#zin_r     = np.fromfile( in_raw, dtype=np.float32,count=1 )[0]
#catalogue = np.fromfile( in_raw, dtype=np.float32,count=outnum*Non_r )
#catalogue = np.reshape( catalogue, (Non_r,outnum) )
## all we care about for the purposes of this script is R_TH
#RTH_r = catalogue[:,6]
#del(catalogue)
### Estimate of halo masses
##M = 4*np.pi/3 * RTH**3 * Omm * 2.775e11*h**2

size_bytes = os.path.getsize(merge_file)

# Merged catalogue
Non_m     = np.fromfile( in_merge, dtype=np.int32,count=1 )[0]
RTHmax_m  = np.fromfile( in_merge, dtype=np.float32,count=1 )[0]
zin_m     = np.fromfile( in_merge, dtype=np.float32,count=1 )[0]

outnum = int( (size_bytes-12)/4/Non_m )

catalogue = np.fromfile( in_merge, dtype=np.float32,count=outnum*Non_m )
catalogue = np.reshape( catalogue, (Non_m,outnum) )
# all we care about for the purposes of this script is R_TH
RTH_m = catalogue[:,6]
del(catalogue)

# Read filter file
filters = np.fromfile( filterfile, sep=' ')
filters = np.reshape( filters[1:], (int(filters[0]),4) )
filters = filters[:,2]

# Make histogram
fig,ax = plt.subplots(nrows=1,ncols=1)
#counts_r,bin_edges = np.histogram(RTH_r,bins=bins)
#bin_centres      = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2
#axs[0].hist( bin_edges[:-1], bin_edges, weights=counts_r,
#    label='pre-merge, {0} peaks'.format(Non_r) )
#counts_m,temp      = np.histogram(RTH_m,bins=bin_edges)
#ax.hist( bin_edges[:-1], bin_edges, weights=counts_m,
#    label='post-merge, {0} peaks'.format(Non_m) )

ax.set_yscale('log')
ax.set_xscale('log')
ax.hist( RTH_m, bins=np.logspace(np.log10(filters[0]), np.log10(RTHmax_m*1.01), 100 ))

ax.set_xlabel(r'$R_{{TH}}$ $[$Mpc$]$')
ax.set_ylabel(r'counts')
ax.legend()

fig.set_size_inches(14,6)
fig.savefig(run_dir+'/hists.png',dpi=100)



