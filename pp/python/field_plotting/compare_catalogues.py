import numpy as np
import matplotlib.pyplot as plt 
import sys

"""
How to make me go:

cd peak-patch
python3 python/run_checks/compare_catalogues.py <path2runs>/<run_1_name> <path2runs/<run_2_name>

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

# Additional parameters
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

in_cat0 = open( '{0}/output/{1}_merge.pksc.{2}'.format(run_dir0,run_name0,seed0), 'rb' )
in_cat  = open( '{0}/output/{1}_merge.pksc.{2}'.format(run_dir, run_name, seed ), 'rb' )

Non0       = np.fromfile( in_cat0, dtype=np.int32,count=1 )[0]
RTHmax0    = np.fromfile( in_cat0, dtype=np.float32,count=1 )[0]
zin0       = np.fromfile( in_cat0, dtype=np.float32,count=1 )[0]
catalogue0 = np.fromfile( in_cat0, dtype=np.float32,count=outnum0*Non0 )
catalogue0 = np.reshape( catalogue0, (Non0,outnum0) )

Non       = np.fromfile( in_cat, dtype=np.int32,count=1 )[0]
RTHmax    = np.fromfile( in_cat, dtype=np.float32,count=1 )[0]
zin       = np.fromfile( in_cat, dtype=np.float32,count=1 )[0]
catalogue = np.fromfile( in_cat, dtype=np.float32,count=outnum*Non )
catalogue = np.reshape( catalogue, (Non,outnum) )


print( '\tfirst run\tsecond run' )
print( 'Non=\t{0}\t{1}'.format(Non0,Non) )
print( 'RTHmax=\t{0}\t{1}'.format(RTHmax0,RTHmax) )
print( 'zin=\t{0}\t{1}'.format(zin0,zin) )

if Non==Non0:
    print( 'average RTH' )
    print( np.average(catalogue0[:,6]) , np.average(catalogue[:,6]) )
    print( )
    print( 'average fpk' )
    print( np.average(catalogue0[:,10]) , np.average(catalogue[:,10]) )
    print( )
    print( 'first few fpk' )
    print( catalogue0[0,10], catalogue[0,10] )
    print( catalogue0[1,10], catalogue[1,10] )
    print( catalogue0[2,10], catalogue[2,10] )
    print( catalogue0[3,10], catalogue[3,10] )



