#!/usr/bin/env python
import numpy as np
import sys
import os
import subprocess

# Usage 
"""
This python runner script creates the initial conditions for a Peak Patch
run based on a standard Peak Patch parameter file. This allows us to
quickly make Gaussian or non-Gaussian random fields without the memory
intensive halo finding and mergine algorithms needing to be run. 

The script is meant to be run in the RandomField directory with the run in
question passed at command line

cd <path to peak patch>/peak-patch/src/modules/RandomField
python3 make_zeta.py <path to runs>/runs/<run name>/

Where the run directory must have at least one folder <run name>/param/
containing a parameter file called param.params, formatted in the usual way
(i.e. see peak-patch/example/param/param.params).
"""
if len(sys.argv) != 2:
    # Error message displayed if wrong # arguments passed at command line
    print('You\'ve passed me the wrong number of arguments at the command '
        +'line.\nCommand line should take the form\ncd <path-to>/peak-patc'
        +'h\npython3 src/modules/RandomField/make_zeta.py <path to>/pkp-ru'
        +'ns/<run name>')
    sys.exit(2) # sys.exit with arg=2 implies command line error

# Reading from command line prompts
run_dir = str(sys.argv[1]).strip() # The directory of the run
if run_dir[-1] != '/':
    run_dir = run_dir+'/' # Ensure directory name ends in "/"


###########################################################################
### Execute useful lines from parameter file                            ###
###########################################################################

# Open parameter file
with open(f'{run_dir}param/param.params') as f:
    params = [i.strip() for i in f.readlines()]

# Execute the parameter file line-by-line
for line in params:
    exec(line)

# Additional parameters
nlattice = int( (nmesh-2*nbuff) * ntile + 2*nbuff )
neff     = int( nlattice - 2*nbuff )


###########################################################################
### Clean up run directory and make subdirectories                      ###
###########################################################################

# Directories that HomogeneousEllipsoid.f90 makes calls to
dir_ = {'home':'./'}
dir_['peak-patch']      = dir_['home']+'../../../'
dir_['src']             = run_dir+'src/'
dir_['RandomField']     = dir_['src']+'modules/RandomField/'
dir_['External']        = dir_['src']+'modules/External/'
dir_['hpkvd']           = dir_['src']+'hpkvd/'
dir_['GlobalVariables'] = dir_['src']+'modules/GlobalVariables/'
dir_['TabInterp']       = dir_['src']+'modules/TabInterp/'

# Defining names of various files
rawpkoutfile  = 'not in use'#f'{run_dir}output/{run_name}_raw.pksc'
hpkvd_inputs  = f'{run_dir}hpkvd_params.txt'
sigmafile     = 'not in use'#f'{run_dir}{run_name}_sigma.dat'
pkfile        = f'{dir_["peak-patch"]}tables/{pkfile}'
NonGauss3file = f'{dir_["peak-patch"]}tables/FNL_spike_w3_piv12.dat'
NonGauss4file = f'{dir_["peak-patch"]}tables/deltaN-LUT-1.875'

# Removie existing subdirectories and make new ones
os.system('cd {0};if [ -d {1} ];then rm -r -f {1};fi'.format(
    run_dir,'src'))
os.system('cp -r {0} {1}'.format(dir_['peak-patch']+'src',run_dir+'src'))
os.system('cd {0};if [ -d {1} ];then rm -r -f {1};fi;mkdir {1}'.format(
    run_dir,'fields'))

# Edit parameter n1 in src/hpkvd/arrays.f90
os.system("sed 's/N_REPLACE/{0}/g' {1}arrays_gen.f90 > {1}arrays.f90".format(n1,dir_['hpkvd']))

# Format and write hierarchical peak-void parameters to text file so they
# can be read by Fortran codes with subroutine read_parameters
os.system('if [ -f "{0}" ];then rm -f {0};fi'.format(hpkvd_inputs))
proc = (subprocess.Popen(["cat"], stdin=subprocess.PIPE,
    stdout=open(hpkvd_inputs,'w')))
stdout, stderr = (proc.communicate(
  input=(
    '%d\n%d\n%f\n%f\n%d\n%f\n%f\n%f\n%f\n%d\n%d\n%d\n%f\n%f\n%f\n%f\n%f\n%d\n%d\n%d\n%d\n%f\n%f\n%f\n%f\n%d\n%d\n%d\n%d\n%f\n%f\n%f\n%f\n%f\n%f\n%d\n%f\n%d\n%d\n%f\n%f\n%f\n%f\n%d\n%d\n%d\n%d\n%d\n%d\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n'%(
      ireadfield,ioutshear,global_redshift,maximum_redshift,num_redshifts,
      Omx,OmB,Omvac,h,nlx,nly,nlz,dcore_box,dL_box,cenx,ceny,cenz,nbuff,
      next,ievol,ivir_strat,fcoll_3,fcoll_2,fcoll_1,dcrit,iforce_strat,
      TabInterpNx,TabInterpNy,TabInterpNz,TabInterpX1,TabInterpX2,
      TabInterpY1,TabInterpY2,TabInterpZ1,TabInterpZ2,wsmooth,rmax2rs,
      ioutfield,NonGauss,fNL,A_nG,B_nG,R_nG,ilpt,iwant_field_part,largerun,seed,ntasks,nmesh,
      fielddir,densfilein,densfileout,pkfile,filterfile,rawpkoutfile,
      TabInterpFile)
    ).encode()
  ))


###########################################################################
### Compile and link Fortran scripts                                    ###
###########################################################################

# Compilers
machine = subprocess.check_output('hostname')
machines = { b'homes-MacBook-Air.local\n' :'NateMacBookAir',
             b'Nates-MacBook-Pro.local\n' :'NateMacBookPro',
             b'nia-login01.scinet.local\n':'Niagara'       ,   
             b'sheep\n'                   :'Sheep'         }
fortrancompilers = { 'NateMacBookAir':'mpifort'   ,
                     'NateMacBookPro':'gfortran'  ,
                     'Niagara'       :'ifort'     ,   
                     'Sheep'         :'gfortran'  }
#ccompilers = { 'NateMacBookAir':'mpicc',
#               'NateMacBookPro':'mpicc',
#               'Niagara'       :'mpicc',
#               'Sheep'         :'mpicc'}
#cppcompilers = { 'NateMacBookAir':'mpiCC',
#                 'NateMacBookPro':'mpiCC',
#                 'Niagara'       :'mpiCC',
#                 'Sheep'         :'mpiCC'}
if machine in machines: 
    fortrancompiler = fortrancompilers[machines[machine]]
else: 
    fortrancompiler = fortrancompilers['NateMacBookAir']

# Copiler options
fortranoptions = '-O4 -w'
    # -Oi is optimization with -O1 being minimum and -O4 being maximum
    # -w flag inhibits all warnings from the compiler

# HomogeneousEllipsoid dependencies
f90s = { 'intreal_types'      : dir_['External']        ,
         'arrays'             : dir_['hpkvd']           ,
         'params'             : dir_['GlobalVariables'] ,
         'cosmoparams'        : dir_['GlobalVariables'] ,
         'input_parameters'   : dir_['RandomField']     ,
         'cosmology'          : dir_['RandomField']     ,
         'textlib'            : dir_['External']        ,
         'random'             : dir_['RandomField']     ,
         'globalvars'         : dir_['RandomField']     ,
         'pktable'            : dir_['RandomField']     ,
         'grid'               : dir_['RandomField']     ,
        #'mpivars'            : dir_['RandomField']     , # has `include 'mpif.h'`
         'mpivars'            : dir_['External']        , # has `include 'mpif.h'`
        #'timing_diagnostics' : dir_['RandomField']     ,
         'timing_diagnostics' : dir_['External']        ,
         'TabInterp'          : dir_['TabInterp']       ,
         'openmpvars'         : dir_['External']        , # has `include omp_lib, omp_lib.f90 is in gfortran
         'fftw_interface'     : dir_['RandomField']     ,
         'gaussian_field'     : dir_['RandomField']     ,
         'tiles'              : dir_['RandomField']     ,
         'chi2zeta'           : dir_['RandomField']     ,
         'RandomField'        : dir_['RandomField']     }

# Set compiler flags
mpi_include  = ''
fftw_include = ''
fftw_library = ''
fftw_flags   = ''
if machines[machine] == 'NateMacBookAir':
    # Tell fortran where to find mpif.h and fftw header files
    fftw_include='-I/usr/local/Cellar/fftw/3.3.10/include'
    mpi_include ='-I/usr/local/Cellar/open-mpi/4.1.2/include'
    # Tell fortran compiler where to find FFTW and Open-MPI libraries
    fftw_lib_dir='-L/usr/local/Cellar/fftw/3.3.10/lib'
    mpi_lib_dir ='-L/usr/local/Cellar/open-mpi/4.1.2/lib'
    # Tell fortran compiler which FFTW libraries to use
    fftw_libs   ='-lfftw3f_mpi -lfftw3f_omp -lfftw3f_threads -lm'
        # -lfftw3f_mpi     library of 32-bit float FTs with MPI
        # -lfftw3f_omp     library of 32-bit float FTs with OpenMP
        # -lfftw3f_threads library of 32-bit float FTs with threading
        # -lm              math library 
compiler_flags='{0} {1} {2} {3} {4} {5}'.format(fortranoptions,mpi_include,
    mpi_lib_dir,fftw_lib_dir,fftw_include,fftw_libs)

# All fortran scripts to compile and link
scripts_f90 = ''
scripts_mod = ''
for i in f90s:
    scripts_f90 += f90s[i]+i+'.f90 ' # f90s & thier paths from RandomField
    scripts_mod += i+'.mod '         # .mod files created in RandomField

# Clean up existing .mod files
os.system('cd {0};rm -f {1}'.format( dir_['home'], scripts_mod ))

# Change mode to avoid issues with Peak Patch .f90 scripts being executable
os.system('cd {0};chmod -x {1}'.format( dir_['home'], scripts_f90 ))

# Compile and Link .f90 files
make_zeta  = 'make_zeta'
executable = make_zeta
os.system('cd {0};{1} {2} {3} {4}.f90 -o {5}'.format(dir_['home'],
        fortrancompiler,compiler_flags,scripts_f90,make_zeta,executable
        )
    )


###########################################################################
### Run make_zeta                                                       ###
###########################################################################

# Run the executable exe_make_zeta and write run directory to command line
os.system('cd {0};./{1} << EOF\n"{2}"\nEOF'.format(
    dir_['home'], executable, run_dir ))

#in_sigma_8 = open( 'sigma_8.dat', 'rb' )
#sigma_8 = np.fromfile( in_sigma_8, dtype=np.float32, count=-1 )
#print(sigma_8)
#os.system('cd {0};rm -f {0}{1}'.format( dir_['home'], 'sigma_8.dat' ))

# Once execution is complete, clean up stray .mod and .h files
os.system('cd {0};rm -f {1}'.format( dir_['home'], scripts_mod ))
#os.system('cd {0};rm -f {1}'.format( dir_['home'], 'mpif.h'    ))

###########################################################################
### Read in zeta                                                        ###
###########################################################################

'''
# Open zeta file
zetafile = f'{run_dir}fields/zetag_{run_name}'
in_zeta = open( zetafile, 'rb' )

# Make zeta into numpy array
zeta = np.reshape(
           np.fromfile(in_zeta,dtype=np.float32,count=-1),
           (nlattice,nlattice,nlattice),
           order='F'
           )[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
'''
#"""
