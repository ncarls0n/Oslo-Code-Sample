#!/usr/bin/env python

# System libraries
import os, subprocess, time

# Math and science libraries
import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Peak Patch libraries
from peakpatchtools import PeakPatch
import peakpatchtools as pkp

####################################################################################################
#                                                                                                  #
#                                 PEAK PATCH INITIALIZATION SCRIPT                                 #
#                                                                                                  #
# This is the main script for seting up Peak Patch and Websky runs. All parameters for a given     #
# Peak Patch-Websky run are set in the a parameter file `<run-directory>/param/param.params`. This #
# script reads in those parameters, copies the source files from `peak-patch/src/` to the run      #
# directory, then compiles and links those source code files.                                      #
#                                                                                                  #
# For quickest implementation on SciNet's Niagara supercomputer, execute this peak-patch.py using  #
# the the shell script `peak-patch/setup.sh`. For more information, see the documentation in       #
# `peak-patch/readme.md`.                                                                          #
#                                                                                                  #
# Note that peak-patch.py is meant to be run in the directory of a Peak Patch run.                 #
#                                                                                                  #
# USAGE (if Peak Patch is in direcotory `<...>/peakpatch/`):                                       #
#     cd <run-directory>                                                                           #
#     python3 <...>/peakpatch/python/peak-patch.py ./param/param.params                            #
#                                                                                                  #
# An example parameter file can be found in `peakpatch/example/param/'.                            #
#                                                                                                  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#---------------------------------------------------------------------------------------------------
# A class for all PeakPatch-realted directories
class PeakPatchLocs:
    def __init__(self, run, pphome):
        self.run_dir = run.run_dir            # The directory where this Peak Patch run is stored
        self.pphome = pphome                  # The Peak Patch home directory
        self.gotopp = f'cd {self.pphome};'    # Command to change direcotry to Peak Patch home
        self.gotorun = f'cd {self.run_dir};'  # Command to change direcotry to this run
        self.src_dir = f'{self.run_dir}/src'  # Source code directory `src' for this run
        self.ppsrc_dir = f'{self.pphome}/src' # Peak Patch home source code directory `src'
        self.gotosrc = f'cd {self.src_dir};'  # Command to change directory to this run `src'
        self.gotoppsrc = f'cd {self.pphome};' # Command to change directory to Peak Patch home `src'

#---------------------------------------------------------------------------------------------------
# A class for filenames of IO data
class IODataFilenames:
    def __init__(self, run):

        # Directory of the Peak Patch run and the name used for it by scheduler
        self.run_dir       = run.run_dir
        self.run_name      = run.run_name

        # Run subdirectory where interpolation tables are stored
        self.tablddir      = f'{self.run_dir}/tables'

        # Run subdirecotry where Peak Patch initial conditions fields are stored
        self.fielddir      = f'{self.run_dir}/{run.fielddir}'

        # File specifying the bank of filter scales used for smoothing and peak finding
        self.filterfile    = f'{self.run_dir}/{run.filterfile}'

        # Interpolation table for Ellipsoidal collapse calculations
        self.TabInterpFile = f'{self.run_dir}/{run.TabInterpFile}'

        # Pre-merge (raw) and merged DM halo catalogue files
        self.rawpkoutfile  = f'{self.run_dir}/output/{self.run_name}_raw.pksc'
        self.mergepkoutfile= f'{self.run_dir}/output/{self.run_name}_merge.pksc'

        # 
        self.sigmafile     = f'{self.run_dir}/output/{self.run_name}_sigma.dat'

        # File containing tabulated power spectra for use in generating initial conditions fields
        self.pkfile        = f'{self.run_dir}/tables/{run.pkfile}'

        # Files containing preheating non-Gaussianity power spectra
        self.NonGauss3file = f'{self.run_dir}/tables/FNL_spike_w3_piv12.dat' # (used if NonGauss==3)
        self.NonGauss4file = f'{self.run_dir}/tables/deltaN-LUT-1.875'       # (used if NonGauss==4)

        # Transfer functions used in primordial intermittent non-Gaussianity (PING) models
        self.T_Wchi2phi_file     = f'{self.run_dir}/tables/T_Wchi2phi.dat'
        self.T_Wchi2phi_file_LSS = f'{self.run_dir}/tables/T_Wchi2phi_LSS_units.dat'
        self.T_phi2zeta_file     = f'{self.run_dir}/src/instability/transfer_func/transfer_func_log_spaced_k.out'
        self.T_phi2zeta_file_LSS = f'{self.run_dir}/tables/Dphi2Dzeta_Tk_LSS_units.dat'

        # Files containing parameters used in Peak Patch hierarchical peak-void code `hpkvd.f90' and
        # DM halo merging and exclusion code `merge_pkvd.f90'
        self.hpkvd_inputs_bin    = f'{self.run_dir}/hpkvd_params.bin'     # binary passed to hpkvd
        self.hpkvd_inputs_ini    = f'{self.run_dir}/param/parameters.ini' # new .ini param file
        self.hpkvd_inputs_params = f'{self.run_dir}/param/param.params'   # legacy .params paramfile
        self.merge_inputs        = f'{self.run_dir}/merge_params.txt'     # passed to merge_pkvd

#---------------------------------------------------------------------------------------------------
#  A class that encompasses all the logging-related variables
class LogVars:
    def __init__(self, run):

        # Directory of the Peak Patch run and the name used for it by scheduler
        self.run_dir = run.run_dir
        self.run_name = run.run_name

        # Log files for makefiles when linking and compiling:

        # 1) The script that generates tables for the ellipsoidal collapse calculation
        self.makeellipsetablog = f'{self.run_dir}/logfiles/mkellipse_{self.run_name}.log'

        # 2) The script that does the Peak Patch hierarchical peak void calculaiton hpkvd.f90
        self.makehpkvdlog = f'{self.run_dir}/logfiles/mkhpvd_{self.run_name}.log'

        # 3) The merging and exclusion script to find the final DM halo catalogue merge_pkvd.f90
        self.makemergelog = f'{self.run_dir}/logfiles/mkmerge_{self.run_name}.log'

        # 4) The script that generates interpolation tables for WebSky map making
        self.makemaptablelog = f'{self.run_dir}/logfiles/mkmaptable_{self.run_name}.log'

        # 5) The scripts that run the WebSky mapmaking algorithm
        self.makepks2maplog = f'{self.run_dir}/logfiles/mkpks2map_{self.run_name}.log'
        self.makepks2cmblog = f'{self.run_dir}/logfiles/mkpks2cmb_{self.run_name}.log'

        # 6) The script that runs the early-universe simulations for PING model non-Gaussianity
        self.makeinstlog = f'{self.run_dir}/logfiles/mkinst_{self.run_name}.log'

        # Log files for running:

        # 1) The script that generates interpollation tables for the homogeneous ellipsoidal
        #    collapse calculation
        self.ellipsetablog = f'{self.run_dir}/logfiles/ellipse_{self.run_name}.log'

        # 2) The script that runs FFTW to generate initial conditions (IC) fields from power spectra
        self.icslogstdout = f'{self.run_dir}/logfiles/ics_{self.run_name}.stdout'

        # 3) hpkvd.f90 - The script that runs the Peak Patch hierarchical peak void calculation
        self.hpkvdlogstdout = f'{self.run_dir}/logfiles/hpkvd_{self.run_name}.stdout'

        # 4) merge_pkvd.f90 - The script that runs the Peak Patch merging and exclusion algorithm
        self.mergelogstdout = f'{self.run_dir}/logfiles/merge_{self.run_name}.stdout'

# ??? - Nathan doesn't know what this is. Something for running stuff on CITA machines
def handle_nix(run):
    nix = os.getenv('NIX_BUILD')
    if nix == '1':
        run.machine = 'nix'

# On some machines, you MUST compile the code on the machine where you are going to run it
# This function returns a flag that says whether your machine requires compilation on itself
# (0) or you can compile it anywhere (1)
def handle_compileless_machines(run):
    compileless_machines = ['starq']
    if run.machine in compileless_machines:
        return 0
    else:
        return 1

#---------------------------------------------------------------------------------------------------
# This just makes sure run parameters aren't going to cause FFTW or hpkvd.f90 to crash
def check_run_parameters_are_allowed(run):
    if run.ireadfield < 3 and run.next % run.ntasks != 0:
        raise ValueError('FFTW requires that the sidelength of the lattice `next` must be evenly\n'
                +'divisible by the number of tasks run in parallel `ntasks`.')
    if run.global_redshift >= run.maximum_redshift:
        if run.num_redshfits == 1:
            run.maximum_redshift = run.global_redshift + 1
        else:
            raise ValueError('global_redshift must be less than or equal to maximum redshift.\n'
                    +'Review parameter file.')

#---------------------------------------------------------------------------------------------------
def copy_source_files_to_local_run_directory_and_run_make_clean(run, pp_dirs):
    gotosrc = pp_dirs.gotosrc
    src_dir = pp_dirs.src_dir
    ppsrc_dir = pp_dirs.ppsrc_dir

    print('Copying source')
    os.system( 'rm -rf {0}; cp -r {1} {0}'.format(src_dir,ppsrc_dir) )
    os.system( '{0} echo "SYSTYPE=\\"{1}\\"" > Makefile.systype'.format(gotosrc,run.machine) )

#---------------------------------------------------------------------------------------------------
# This function creates subdirectories in the run directory where logfiles, binary executables, 
# interpollation tables and simulation products (outputs) are saved during the Peak Patch
# calculation. If needed, additional directories for Websky maps and initial conditions fields are
# also created.
def create_subdirectories_in_the_run_directory(run, pp_dirs):
    gotorun = pp_dirs.gotorun
    os.system( '{0}if [ ! -e logfiles ] ; then mkdir logfiles ; fi'.format(gotorun) )
    os.system( '{0}if [ ! -e bin ]      ; then mkdir bin      ; fi'.format(gotorun) )
    os.system( '{0}if [ ! -e tables ]   ; then mkdir tables   ; fi'.format(gotorun) )
    os.system( '{0}if [ ! -e output ]   ; then mkdir output   ; fi'.format(gotorun) )
    if run.map_params==1:
        os.system('{0}if [ ! -e maps ]  ; then mkdir maps  ; fi'.format(gotorun) )
    if ( run.ioutfield>0 or run.ireadfield>0 or run.NonGauss>0 ):
        os.system('{0}mkdir -p {1}'.format(gotorun,run.fielddir))

#---------------------------------------------------------------------------------------------------
# This function makes shell scripts for Niagara jobs to generate WebSky mock maps for any response
# functions specified in the parameter file. This is limited to the response functions encompassed
# by the pks2map.f90 script, and does not cover other work done with LIMLAM Movker or XG Paint
# projects, also under the WebSky umbrella.
def make_list_of_map_profiles_to_run(run):
    profiles     = run.maps.split()
    makemapsname = run.maps.split()
    for i in range(len(makemapsname)):
        makemapsname[i] = 'make_'+makemapsname[i]+'_map_'+run.short_name+'.sh'
    return profiles, makemapsname

#---------------------------------------------------------------------------------------------------
# Copy generic Peak Patch run submit scripts and tables to run directory
def copy_pp_run_submit_scripts_and_tables(run, io_files, pp_dirs, batchfname):
    pphome = pp_dirs.pphome
    pkfile = io_files.pkfile
    NonGauss3file = io_files.NonGauss3file
    NonGauss4file = io_files.NonGauss4file

    # Copy Peak Patch run submit script to run directory
    os.system('cp {0}/templates/ppatchrun_{1}-{2}.sh ./{3}'
        .format(pphome,run.machine,run.runtype,batchfname))

    # Copy power spectrum table to run directory
    os.system('cp {0}/tables/{1} {2}'.format( pphome , run.pkfile , pkfile ))

    # Copy tables for pre-heating model non-Gaussianity
    os.system('cp {0}/tables/{1} {2}'.format(pphome,os.path.basename(NonGauss3file),NonGauss3file))
    os.system('cp {0}/tables/{1} {2}'.format(pphome,os.path.basename(NonGauss4file),NonGauss4file))

#---------------------------------------------------------------------------------------------------
# This function runs the python script `filter_gen.py' which generates a bank of filter scales used
# in the smoothed peak finding phase of the Peak Patch calculation
def generate_filter_bank(run, pp_dirs):
    pphome = pp_dirs.pphome
    print('Generating filter bank')
    os.system( 'python3 {0}/python/filter_gen.py {1} {2} {3} > {4}/{5}'
        .format( pphome, run.dL_box, run.n1, run.Rsmooth_max, run.run_dir, run.filterfile ))

#---------------------------------------------------------------------------------------------------
# This function calculates linear-theory matter power spectra for the cosmology specified by
# cosmological parameters set in the parameter file
def calculate_new_power_spectrum(run):
    run.pkfile = '{0}/tables/{1}.dat'.format( run.run_dir, run.short_name )
    if os.path.isfile(run.pkfile):
        os.system('rm -rf {0}'.format( run.pkfile ))
    pkp.powerspectrum_create( outfile = run.pkfile,
        z=0.0, k_min=2*np.pi/run.boxsize/10, k_max=2*np.pi/run.cellsize*10, nkpoints=1000,
        Omega_b   = run.OmB ,  Omega_CDM = run.OmB    ,  Omega_k   = 0.0     ,
        h         = run.h   ,  sigma_8   = run.sigma8 ,  tau       = run.tau ,
        m_nu      = 0.0     ,  n_s       = run.ns     ,  A_s       = run.As   )

#---------------------------------------------------------------------------------------------------
# This function runs the fortran scripts for Early Universe calculations used in primordial
# intermittent non-Gaussianity (PING) model
def perform_calculations_for_early_universe_non_gaussian_ics(run, pp_dirs, io_files, log_vars):

    # Specify directory paths where early universe interpollation tables are to be saved
    src_dir             = pp_dirs.src_dir
    makeinstlog         = log_vars.makeinstlog
    T_phi2zeta_file     = io_files.T_phi2zeta_file
    T_phi2zeta_file_LSS = io_files.T_phi2zeta_file_LSS
    T_Wchi2phi_file     = io_files.T_Wchi2phi_file
    T_Wchi2phi_file_LSS = io_files.T_Wchi2phi_file_LSS

    # Additional fixed parameters
    def kj(k0,kn,j,n): return k0*(kn/k0)**(j/n)
    ins_dir = src_dir+'/instability/'
    lcode        = 2.6259e-52         # Mpc/(10^5 reduced planck masses)
    acode_approx = 3                  # the approximate early u scale factor
    nk   = 100
    k0   = 2 * np.pi            / run.a_e * acode_approx * lcode / run.boxsize / 10
    kn   = (2*run.next+1)*np.pi / run.a_e * acode_approx * lcode / run.boxsize * 10
    dk   = kj(k0,kn,-1,nk)
    kmax = kj(k0,kn,nk+1,nk)

    # Run the early universe instability code, the outputs are saved as the non-Gaussianity
    # parameters from models 5 and 6, we reuse these since they are floats that are already being
    # passed ahead to the Fortran source code.
    inst1, inst2, inst3 = pkp.run_instability( ins_dir   = ins_dir,
                                               a_e       = run.a_e,
                                               boxsize   = run.boxsize,
                                               m_phi     = run.m_phi,
                                               m_chi     = run.m_chi,
                                               phi_w     = run.phi_w,
                                               phi_p     = run.phi_p,
                                               m_tach    = run.m_tach,
                                               vev       = run.vev,
                                               nstep     = 2**15,
                                               stepadapt = 1,
                                               dt0       = 2.0**-5,
                                               nk        = nk,
                                               k0        = k0,
                                               kn        = kn,
                                               save_out  = True,
                                               overwrite = False,
                                               ng_model  = run.NonGauss,
                                               strength  = -1
                                             )

    # For Non-Gaussianity models 7-10
    if run.NonGauss in [ 7, 8, 9, 10 ]:

        # Interpret output from the instability code as logarithm of scale factor at the end of the
        # instability, Hubble parameter at the end of the instability, and the fit parameter
        # relating \chi^2(H_e) to \phi(H_e)
        alpha_e            = inst1
        H_e                = inst2
        quadratic_coupling = inst3

        # The coupling constant for the instability
        lambda_chi = run.m_tach**2 * run.vev**-2

        # Save the quadratic fit for \chi_e to \Delta\phi_e to the non-Gaussianity parameters that
        # get passed along to the fortran code
        run.A_nG = quadratic_coupling
        run.B_nG = alpha_e
        run.R_nG = lambda_chi
        run.A_nG, run.B_nG, run.R_nG = inst1, inst2, inst3
        print( '\\Delta\\phi = {0}\\chi^2 - <\\chi^2>'.format( run.A_nG ) )
        
    # For Non-Gaussianity model 11
    elif run.NonGauss in [ 11 ]:
        
        # Interpret output from the instability code as logarithms of scale factors at beginning
        # (i), middle (p) and end (e) of the instability, similarly for Hubble parameter, and the
        # amplification factor in the \Delta\phi_e to \Delta\zeta_f transfer function
        alpha_i , alpha_p , alpha_e = inst1
        H_i     , H_p     , H_e     = inst2
        strength                    = inst3

        # Save Non-Gaussianity model 11 outputs
        np.savetxt( ins_dir+'ng_model_11_outputs.dat',
                np.array([ alpha_i, alpha_p, alpha_e, H_i, H_e, H_p, strength ]) ) 

        # This will at some point be changed to a table read, but for now it's just a top-hat filter
        # Read in the (W*\chi_e)^2 to \Delta\phi_e transfer function table, we use the same k range
        # as the \Delta\phi_e to \Delta\zeta_f transfer function
        t_Wchi2Dphi       = np.loadtxt(T_phi2zeta_file)
        t_Wchi2Dphi[:,1] *= 0.0

        # The first column is k/a_eH_e, so we first multiply by a_e H_e to get k
        t_Wchi2Dphi[:,0] *= np.exp(alpha_e) * H_e

        # Next, we convert to LSS units
        t_Wchi2Dphi[:,0] *= run.a_e * np.exp( -alpha_e ) / lcode

        # Set up the (W*\chi_e)^2 to \Delta\phi_e transfer function as a step function
        k_instability     = np.exp(alpha_p) * np.sqrt( run.m_tach**2 - run.m_chi**2 )
        k_instability    *= run.a_e * np.exp( -alpha_e ) / lcode
        t_Wchi2Dphi[:,1]  = np.heaviside( 4 * k_instability - t_Wchi2Dphi[:,0] , 0 )
        
        # Save the (W*\chi_e)^2 to \Delta\phi_e transfer function in LSS units
        np.savetxt( T_Wchi2phi_file_LSS, t_Wchi2Dphi )

    ################################################################################################
    ### CONVERT T_{\Delta\phi_e\to\Delta\zeta_f}(k) TO UNITS USED BY Peak Patch, OUTPUT TO  FILE ###
    if run.NonGauss in [ 9, 10, 11 ]:

        # Read in the \Delta\phi_e to \Delta\zeta_f transfer function
        t_Dphi2Dzeta = np.loadtxt(T_phi2zeta_file)

        # The first column is k/a_eH_e, so we first multiply by a_e H_e to get k
        t_Dphi2Dzeta[:,0] *= np.exp(alpha_e) * H_e

        # Next, we convert to LSS units
        t_Dphi2Dzeta[:,0] *= run.a_e * np.exp( -alpha_e ) / lcode

        # Multiply the amplification factor to the \Delta\phi_e to \Delta\zeta_f transfer function
        if run.NonGauss == 11:
            t_Dphi2Dzeta[:,1] *= run.A_nG

            # Add the (W*\chi_e)^2 to \Delta\phi_e transfer function as a third column
            t_Dphi2Dzeta = np.column_stack(( t_Dphi2Dzeta, t_Wchi2Dphi[:,1] ))

        # Save the \Delta\phi_e to \Delta\zeta_f transfer function
        np.savetxt( T_phi2zeta_file_LSS, t_Dphi2Dzeta )

    return H_e

#---------------------------------------------------------------------------------------------------
# For purely Gaussian initial conditions, the additional parameters used by primordial non-
# Gaussianity models are not needed and set to 0 by this function.
def zero_out_nongaussianity_parameters(run):
    # Set unused parameters to 0
    H_e = 0.0

    # For backwards compatibility, any of the more recently developed non-Gaussianity parameters are
    # not included in the parameter file already are set to zero
    nong_params = [ 'A_nG', 'B_nG', 'R_nG', 'm_phi', 'm_chi', 'phi_w', 'phi_p', 'vev', 'm_tach',
                    'a_e', 'chi_seed']
    for nong_param in nong_params:    
        if nong_param not in dir(run):
            vars(run)[nong_param] = 0.0
    return H_e

#---------------------------------------------------------------------------------------------------
# Compile hpkvd (hierarchical peak-finding script) from source
def compile_hpkvd(run, pp_dirs, log_vars, num_compile_jobs):
    src_dir = pp_dirs.src_dir
    gotosrc = pp_dirs.gotosrc
    makehpkvdlog = log_vars.makehpkvdlog

    print('Compiling hpkvd')
    os.system('cd {0}/hpkvd;sed \'s/N_REPLACE/{1}/g\' arrays_gen.f90 > arrays.f90'
              .format(src_dir,run.n1) )
    os.system(gotosrc+'make clean &> /dev/null')
    os.system(gotosrc+'make -j'+num_compile_jobs+' hpkvd &> '+makehpkvdlog)

#---------------------------------------------------------------------------------------------------
# Run hpkvd using .ini parameter file
def run_hpkvd_or_make_input_file_ini(run, io_files, pp_dirs, log_vars):
    gotorun        = pp_dirs.gotorun
    hpkvd_inputs   = io_files.hpkvd_inputs_ini
    hpkvdlogstdout = log_vars.hpkvdlogstdout

    if run.batch==0:
        print('Running hpkvd')    
        os.system('{0}mpirun -np {1} ./bin/hpkvd 1 {2} {3}>{4}'.format(gotorun,run.ntasks,run.seed,
                                                                       hpkvd_inputs,hpkvdlogstdout))

#---------------------------------------------------------------------------------------------------
# Run hpkvd using legacy .params parameter file
def run_hpkvd_or_make_input_file_bin(run, io_files, pp_dirs, log_vars, H_e):

    # Defining some convenience varables from the above classes
    gotorun             = pp_dirs.gotorun
    fielddir            = io_files.fielddir
    pkfile              = io_files.pkfile
    filterfile          = io_files.filterfile
    rawpkoutfile        = io_files.rawpkoutfile
    TabInterpFile       = io_files.TabInterpFile
    hpkvd_inputs        = io_files.hpkvd_inputs_bin
    T_phi2zeta_file_LSS = io_files.T_phi2zeta_file_LSS
    hpkvdlogstdout      = log_vars.hpkvdlogstdout

    # List of hpkvd input parameters, formatted in the approprate format (32-bit numbers, strings
    # converted to bytes preceded by a 32-bit integer telling hpkvd how many bytes the subsequent
    # string will be)
    hpkvd_list = [ np.int32(     run.ireadfield       ) , np.int32(   run.ioutshear               ),
                   np.float32(   run.global_redshift  ) , np.float32( run.maximum_redshift        ),
                   np.int32(     run.num_redshifts    ) , np.float32( run.Omx                     ),
                   np.float32(   run.OmB              ) , np.float32( run.Omvac                   ),
                   np.float32(   run.h                ) , np.int32(   run.nlx                     ),
                   np.int32(     run.nly              ) , np.int32(   run.nlz                     ),
                   np.float32(   run.dcore_box        ) , np.float32( run.dL_box                  ),
                   np.float32(   run.cenx             ) , np.float32( run.ceny                    ),
                   np.float32(   run.cenz             ) , np.int32(   run.nbuff                   ),
                   np.int32(     run.next             ) , np.int32(   run.ievol                   ),
                   np.int32(     run.ivir_strat       ) , np.float32( run.fcoll_3                 ),
                   np.float32(   run.fcoll_2          ) , np.float32( run.fcoll_1                 ),
                   np.float32(   run.dcrit            ) , np.int32(   run.iforce_strat            ),
                   np.int32(     run.TabInterpNx      ) , np.int32(   run.TabInterpNy             ),
                   np.int32(     run.TabInterpNz      ) , np.float32( run.TabInterpX1             ),
                   np.float32(   run.TabInterpX2      ) , np.float32( run.TabInterpY1             ),
                   np.float32(   run.TabInterpY2      ) , np.float32( run.TabInterpZ1             ),
                   np.float32(   run.TabInterpZ2      ) , np.int32(   run.wsmooth                 ),
                   np.float32(   run.rmax2rs          ) , np.int32(   run.ioutfield               ),
                   np.int32(     run.NonGauss         ) , np.float32( run.fNL                     ),
                   np.float32(   run.A_nG             ) , np.float32( run.B_nG                    ),
                   np.float32(   run.R_nG             ) , np.float32( H_e                         ),
                   np.int32(     run.chi_seed         ) , np.int32(   run.ilpt                    ), 
                   np.int32(     run.iwant_field_part ) , np.int32(   run.largerun                ),
                   np.int32(len( fielddir             )), bytes(      fielddir            ,'ascii'),
                   np.int32(len( run.densfilein       )), bytes(      run.densfilein      ,'ascii'),
                   np.int32(len( run.densfileout      )), bytes(      run.densfileout     ,'ascii'),
                   np.int32(len( pkfile               )), bytes(      pkfile              ,'ascii'),
                   np.int32(len( T_phi2zeta_file_LSS  )), bytes(      T_phi2zeta_file_LSS ,'ascii'),
                   np.int32(len( filterfile           )), bytes(      filterfile          ,'ascii'),
                   np.int32(len( rawpkoutfile         )), bytes(      rawpkoutfile        ,'ascii'),
                   np.int32(len( TabInterpFile        )), bytes(      TabInterpFile       ,'ascii')]

    # Write hpkvd inputs to a binary file
    with open(hpkvd_inputs,'wb') as stdout:
        for i in hpkvd_list:
            stdout.write( i )

    # Run hpkvd if batch is set to 0 in the parameter file (the default is batch=1 in which case a
    # job script is made but hpkvd is not submitted).
    if run.batch==0:
        print('Running hpkvd')    
        os.system('{0}mpirun -np {1} ./bin/hpkvd 1 {2} {3}>{4}'.format(gotorun,run.ntasks,
                                                         run.seed,hpkvd_inputs,hpkvdlogstdout))


#---------------------------------------------------------------------------------------------------
# Compile merge_pkvd which runs the exclusion and merge algorithm generating a final halo catalogue
def compile_merge_pkvd(pp_dirs, log_vars, num_compile_jobs):
    gotosrc = pp_dirs.gotosrc
    makemergelog = log_vars.makemergelog

    print('Compiling merge_pkvd')    
    os.system('{0}make -j{1} merge_pkvd &> {2}'.format(gotosrc,num_compile_jobs,makemergelog))

#---------------------------------------------------------------------------------------------------
# Run merging and exclusion algorithm
def run_merge_pkvd_or_create_input_file_(run, pp_dirs, log_vars, io_files):
    gotorun = pp_dirs.gotorun
    mergelogstdout = log_vars.mergelogstdout
    rawpkoutfile = io_files.rawpkoutfile
    merge_inputs = io_files.merge_inputs
    mergepkoutfile = io_files.mergepkoutfile
    
    proc = ( subprocess.Popen( ["cat"] , stdin=subprocess.PIPE , stdout=open(merge_inputs,'w') ) )
    end=""
    stdout, stderr = (proc.communicate(input=(
    "%f\n%d\n%d\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%f\n%d\n%s\n%s\n"%(
        run.boxsize          , run.ntilemerge       , run.ievol            , run.cenx             , 
        run.ceny             , run.cenz             , run.maximum_redshift , run.Omx              ,
        run.OmB              , run.Omvac            , run.h                , run.iZeld            ,
        run.iLexc            , run.iLmrg            , run.iFexc            , run.iFmrg            ,
        run.ioutshear        , run.iwrap            , run.iwant_field_part , run.mlatt            ,
        run.largerun         , rawpkoutfile         , mergepkoutfile
        )).encode()))

    # merge_pkvd_list = [ np.float32( run.boxsize          ), np.int32(   run.ntilemerge ),
    #                     np.int32(   run.ievol            ), np.float32( run.cenx       ),
    #                     np.float32( run.ceny             ), np.float32( run.cenz       ),
    #                     np.float32( run.maximum_redshift ), np.float32( run.Omx        ),
    #                     np.float32( run.OmB              ), np.float32( run.Omvac      ),
    #                     np.float32( run.h                ), np.int32(   run.iZeld      ),
    #                     np.int32(   run.iLexc            ), np.int32(   run.iLmrg      ),
    #                     np.int32(   run.iFexc            ), np.int32(   run.iFmrg      ),
    #                     np.int32(   run.ioutshear        ), np.int32(   run.iwrap      ),
    #                     np.int32(   run.iwant_field_part ), np.float32( run.mlatt      ),
    #                     np.int32(   run.largerun         ),
    #                     # 32-bit str? I don't know? ( rawpkoutfile   )
    #                     # 32-bit str? I don't know? ( mergepkoutfile )
    #                     # Then two strings which I'll have to change merge_pkvd to do the same thing as hpkvd and I don't have the energy for that right now
    # 
    # with open(merge_inputs,'wb') as stdout:
    #     for i in merge_pkvd_list:
    #         stdout.srite(i)
    
    if run.batch==0:
        print('Running merge_pkvd')
        os.system('{0}./bin/merge_pkvd 1>{1}'.format(gotorun,mergelogstdout))


#---------------------------------------------------------------------------------------------------
# Compiles `pks2map.f90' which is the portion of the WebSky pipeline built-in to the Peak Patch
# repository.
def compile_the_Websky_pipeline_pks2map(pp_dirs, log_vars, num_compile_jobs):
    gotosrc        = pp_dirs.gotosrc
    makepks2maplog = log_vars.makepks2maplog

    # Compile pks2map
    print('Compiling pks2map')
    os.system('{0}make -j{1} pks2map &> {2}'.format(gotosrc,num_compile_jobs,makepks2maplog))
    #os.system(gotosrc+'make pks2cmb > '+makepks2cmblog)
    #os.system(gotosrc+'make make_maptable > '+makemaptablelog)


#---------------------------------------------------------------------------------------------------
# This function runs each of the above compilation functions if the appropriate switch is set to 1
# in the parameter file.
def compile_everything(run, io_files, pp_dirs, log_vars, num_compile_jobs, H_e, legacy_input):
    num_compile_jobs = str(num_compile_jobs)

    #--------------------- HPKVD ---------------------------
    if run.compile_hpkvd==1:
        compile_hpkvd(run, pp_dirs, log_vars, num_compile_jobs)
    
    if run.hpkvd_params==1:
        if legacy_input==1:
            run_hpkvd_or_make_input_file_bin(run, io_files, pp_dirs, log_vars, H_e)
        else:
            run_hpkvd_or_make_input_file_ini(run, io_files, pp_dirs, log_vars)
    
    #------------------- MERGE PKVD ------------------------

    # Compile merge_pkvd which runs the exclusion and merge algorithm generating a final halo catalogue
    if run.compile_merge==1:
        compile_merge_pkvd(pp_dirs, log_vars, num_compile_jobs)
    
    if run.merge_params==1:
        run_merge_pkvd_or_create_input_file_(run, pp_dirs, log_vars, io_files)

    #---------------------- WEBSKY ---------------------------
           
    if run.compile_maps==1:
        compile_the_Websky_pipeline_pks2map(pp_dirs, log_vars, num_compile_jobs)


#---------------------------------------------------------------------------------------------------
# Getting compilation steps if need to compile on the target machine
def pass_compilation_to_submit_script( run, pp_dirs, log_vars, io_files, legacy ):

    print('Bypassing compilation - it will be done with the submission script')
    gotosrc = pp_dirs.gotosrc
    gotorun = pp_dirs.gotorun
    if legacy == 1:
        config_file = io_files.hpkvd_inputs_params
    elif legacy == 0:
        config_file = io_files.hpkvd_inputs_ini

    # List to store all relevant compilation steps
    compilation_steps = ['# These are the compilation steps that are required to be done on the target machine']

    if run.compile_hpkvd == 1:
        src_dir = pp_dirs.src_dir
        n1 = run.n1
        makehpkvdlog = log_vars.makehpkvdlog

        # Replicating compile_hpkvd function
        compilation_steps.append(f'cd {src_dir}/hpkvd;sed \'s/N_REPLACE/{n1}/g\' arrays_gen.f90 > arrays.f90')
        compilation_steps.append(f'{gotosrc} make clean &> /dev/null')
        compilation_steps.append(f'{gotosrc} make hpkvd CONFIG_FILE = {config_file} &> {makehpkvdlog}')

    if run.compile_merge == 1:
        makemergelog = log_vars.makemergelog

        compilation_steps.append(f'{gotosrc} make merge_pkvd &> {makemergelog}')

    if run.compile_maps == 1:
        makepks2maplog = log_vars.makepks2maplog

        # Compile pks2map
        compilation_steps.append(f'{gotosrc} make pks2map &> {makepks2maplog}')

    # Return to the run directory
    compilation_steps.append(f'{gotorun}')

    return compilation_steps


#-------------------------------------------------------------------------------------------------
# This function is used to create input files when hpkvd/merge_pkvd weren't compiled
def create_input_files_only(run, io_files, pp_dirs, log_vars):

    # Force NOT to run hpkvd/merge files, since they wouldn't be compiled
    batch_stored_value = run.batch
    run.batch = 1

    # Only creates the input files, does not run appropriate programs
    if run.hpkvd_params==1:
        if legacy_input==1:
            run_hpkvd_or_make_input_file_bin(run, io_files, pp_dirs, log_vars, H_e)
        else:
            run_hpkvd_or_make_input_file_ini(run, io_files, pp_dirs, log_vars)

    if run.merge_params==1:
        run_merge_pkvd_or_create_input_file_(run, pp_dirs, log_vars, io_files)

    run.batch = batch_stored_value 


#-------------------------------------------------------------------------------------------------
# Assemble a submission script to contain the compilation steps
def insert_compilation_steps(batchfname, compilation_steps):

    # Read the batch file content
    with open(batchfname, 'r') as file:
        batch_content = file.readlines()

    # Insert each step from the list above the COMPILATION_REPLACE placeholder
    new_content = []
    for line in batch_content:
        if 'COMPILATION_REPLACE' in line:
            for step in compilation_steps:
                new_content.append(step + '\n')
        else:
            new_content.append(line)

    # Write the modified content back to the batch file
    with open(batchfname, 'w') as file:
        file.writelines(new_content)


#---------------------------------------------------------------------------------------------------
def make_Peak_Patch_run_submit_scripts( run, log_vars, pp_dirs, batchfname, makemapsname, profiles,
        compilation_steps ):

    gotorun        = pp_dirs.gotorun
    pphome         = pp_dirs.pphome
    hpkvdlogstdout = log_vars.hpkvdlogstdout
    icslogstdout   = log_vars.icslogstdout
    mergelogstdout = log_vars.mergelogstdout

    # Use Linux Stream Editor to find and replace placeholder strings in generic submit scripts with
    # run-specific values from parameter file
    print('Creating batch file for parallel run')
    s = gotorun+"sed 's!{0}!{1}!g' {2} > temp;mv temp {2}"
    os.system( s.format( 'NNODES_REPLACE'      , run.nnodes            , batchfname ))
    os.system( s.format( 'TPNODE_REPLACE'      , run.tpnode            , batchfname ))
    os.system( s.format( 'NTASKS_REPLACE'      , run.nnodes*run.tpnode , batchfname ))
    os.system( s.format( 'TLIMIT_REPLACE'      , run.tlimit            , batchfname ))
    os.system( s.format( 'SNAME_REPLACE'       , run.short_name        , batchfname ))
    os.system( s.format( 'LNAME_REPLACE'       , run.run_name          , batchfname ))
    os.system( s.format( 'MNAME_REPLACE'       , run.run_name          , batchfname ))
    os.system( s.format( 'RAPI_REPLACE'        , run.rapi              , batchfname ))
    os.system( s.format( 'NTASKS_REPLACE'      , run.ntasks            , batchfname ))
    os.system( s.format( 'NTASKSMERGE_REPLACE' , run.ntasksmerge       , batchfname ))
    os.system( s.format( 'SEED_REPLACE'        , run.seed              , batchfname ))
    os.system( s.format( 'NCPUS_REPLACE'       , run.ncpus             , batchfname ))
    os.system( s.format( 'BOXSIZE_REPLACE'     , run.boxsize           , batchfname ))
    os.system( s.format( 'NMESH_REPLACE'       , run.nmesh             , batchfname ))
    os.system( s.format( 'NBUFF_REPLACE'       , run.nbuff             , batchfname ))
    os.system( s.format( 'NTILE_REPLACE'       , run.ntile             , batchfname ))
    os.system( s.format( 'NTILE_REPLACE'       , run.ntilemerge        , batchfname ))
    os.system( s.format( 'HPKVDSTDOUT_REPLACE' , hpkvdlogstdout        , batchfname ))
    os.system( s.format( 'ICSSTDOUT_REPLACE'   , icslogstdout          , batchfname ))
    os.system( s.format( 'MERGESTDOUT_REPLACE' , mergelogstdout        , batchfname ))

    # Compilation on a target machine (if needed)
    insert_compilation_steps(batchfname,compilation_steps)

    # Submit Peak Patch run jobscript to scheduler (usually set to 0)
    if run.submit == 1:
        print('Submitting batch file for parallel run')
        os.system(gotorun+run.submit_command+' '+batchfname)

    # Make Websky run submit scripts
    if run.map_params == 1: 
        print('Creating batch file for mapmaking')

        # Copy BBPS pressure profile tables
        os.system('cp '+pphome+'/tables/'+run.tabfile_map+' ./tables') 
        os.system('cp '+pphome+'/tables/'+run.tabfile_sfr+' ./tables')

        # Go to this run's maps directory and copy of the WebSky submit script there
        gotomaps = 'cd {0}/maps;'.format(run.run_dir)
        mapscript = 'make_MODEL_map_{0}.sh'.format(run.short_name)
        os.system( '{0}cp {1}/templates/webskyrun_niagara-single.sh ./{2}'
            .format(gotomaps,pphome,mapscript))

        # Use Steam Editor to find and replace placeholders in generic map-making scripts with run-
        # specific values from parameter file
        m = gotomaps+"sed 's/{0}/{1}/g' {2} > temp;mv temp {2}"
        os.system( m.format( 'NNODESMAP_REPLACE'  , run.nnodes_map   , mapscript ))
        os.system( m.format( 'NTASKSMAP_REPLACE'  , run.ntasks_map   , mapscript ))
        os.system( m.format( 'TPNODEMAP_REPLACE'  , run.ppn_map      , mapscript ))
        os.system( m.format( 'NCPUSMAP_REPLACE'   , run.np_map       , mapscript ))
        os.system( m.format( 'TLIMITMAP_REPLACE'  , run.tlimit_map   , mapscript ))
        os.system( m.format( 'TABFILEMAP_REPLACE' , run.tabfile_map  , mapscript ))
        os.system( m.format( 'ZMINMAP_REPLACE'    , run.zmin_map     , mapscript ))
        os.system( m.format( 'NSIDEMAP_REPLACE'   , run.nside_map    , mapscript ))
        os.system( m.format( 'SCRAMBLEMAP_REPLACE', run.scramble_map , mapscript ))
        os.system( m.format( 'CENTERMAP_REPLACE'  , run.center_map   , mapscript ))
        os.system( m.format( 'NPIXMAP_REPLACE'    , run.npix_map     , mapscript ))
        os.system( m.format( 'FOVMAP_REPLACE'     , run.fov_map      , mapscript ))
        os.system( m.format( 'ZMAXMAP_REPLACE'    , run.zmax_map     , mapscript ))
        os.system( m.format( 'CHIHVIEWMAP_REPLACE', run.chihview_map , mapscript ))
        os.system( m.format( 'MODELMAP_REPLACE'   , run.model_map    , mapscript ))
        os.system( m.format( 'PSZCUTMAP_REPLACE'  , run.PSZcut_map   , mapscript ))
        os.system( m.format( 'SEEDMAP_REPLACE'    , run.seed         , mapscript ))  
        os.system( m.format( 'OUTFILEMAP_REPLACE' , run.run_name     , mapscript ))

        # Make separate Websky submit script for each model
        for i in range(len(makemapsname)):
            os.system( '{0}cp {1} {2}'.format( gotomaps,mapscript,makemapsname[i] ) )
            os.system( m.format( 'PROFILE_REPLACE' , profiles[i] , makemapsname[i] ) )
            #os.system('{0}chmod +x {1}'.format(gotomaps,makemapsname[i]))

        # Delete copy of generic Websky submit script from <run-dir>/map/
        os.system( '{0}rm {1}'.format(gotomaps,mapscript))


####################################################################################################
######################################### MAIN FUNCTION ############################################
####################################################################################################

def main():

    #----------- VARIABLES DEFINITIONS BEGIN ------------

    # Set path to the Peak Patch directory
    pphome = os.environ["PP_DIR"] 

    # Legacy input, to use python-created .bin file (1) instead of .ini file (0)
    # .ini is preferred since in that case the file gets passed directly into PeakPatch
    legacy_input = 1
    
    # Start timer
    starttime=time.time()
    
    # Load Peak Patch run as PeakPatch class object
    run = PeakPatch('.')
    
    # Define command to change directories to Peak Patch, the Peak Patch source files, and the copy of
    # the source files that will be put in the run directory
    pp_dirs = PeakPatchLocs(run, pphome)

    # A check to use nix config if needed
    handle_nix(run)

    # Set of filenames of input/output data, as a class 'Input/Output Data Filenames'
    io_files = IODataFilenames(run)

    # Set names of log files, standard output files, as a class 'Logging Variables'
    log_vars = LogVars(run)
    
    # Set name of Peak Patch job submit script
    batchfname = '{0}_{1}.sh'.format( run.short_name , run.seed )
    
    # Make list of map profiles to run
    profiles, makemapsname = make_list_of_map_profiles_to_run(run)

    #------------- VARIABLES DEFINITIONS END -------------

    # Error handling
    check_run_parameters_are_allowed(run)
    
    # Update short_name to reflect initial seed
    if run.short_name[ run.short_name.rfind('_') : ] != '_'+str(run.seed):
        run.short_name += '_'+str(run.seed)

    #-------------------- COPY SOURCE ---------------------

    if run.compile_hpkvd == 1 or run.compile_merge == 1 or run.compile_maps == 1:    
        copy_source_files_to_local_run_directory_and_run_make_clean(run, pp_dirs)
    
    create_subdirectories_in_the_run_directory(run, pp_dirs)
    
    # Copy generic Peak Patch run submit scripts and tables to run directory
    if run.batch == 1:
        copy_pp_run_submit_scripts_and_tables(run, io_files, pp_dirs, batchfname)

    #-------------- CREATE TABLES/POWERSPEC ----------------
    
    # Generate Filter Bank (list of smoothing scales used in peak finding)
    if run.create_filterbank==1:
        generate_filter_bank(run, pp_dirs)
    
    # If pkfile = 0 or None, a new power spectrum is calculated for this cosmology
    if run.pkfile == '0' or run.pkfile.lower() == 'none' or run.pkfile == None:
        calculate_new_power_spectrum(run)

    immediate_compilation = handle_compileless_machines(run)
    
    #---------------- NONGAUSSIANITIES ---------------------
    
    # Compile and run lattice simulation intitial conditions
    # THIS REQUIRES A FULL REWRITE TO WORK FINE WITH ANY NEW NATE'S EDITS.
    # Brought up to date with ddee992eee71eaebef3ad4a75ebe2ebb815f89fb
    # (first commit of June 12th)
    # Any function should be identical to that commit
    # ALSO, MAKEFILES REQUIRE SOME WORK ON THEM
    if run.NonGauss>=7 and run.NonGauss<=11:
        H_e = perform_calculations_for_early_universe_non_gaussian_ics(run, pp_dirs, io_files, log_vars)
    else:
        H_e = zero_out_nongaussianity_parameters(run)
    
    #---------------- COMPILATION ---------------------

    num_compile_jobs = 1 # Would have to manually adjust this for now unfortunately

    if immediate_compilation==1:
        compile_everything(run, io_files, pp_dirs, log_vars, num_compile_jobs, H_e, legacy_input)
        compilation_steps = []
    else:
        # TODO: Add compilation parallelization here
        compilation_steps = pass_compilation_to_submit_script(run, pp_dirs, log_vars, io_files,
                legacy_input) 
        create_input_files_only(run, io_files, pp_dirs, log_vars)

    #--------------- SUBMISSION SCRIPTS --------------------

    # Untested after rewrite, but should function identically
    if run.batch==1:
        make_Peak_Patch_run_submit_scripts(run, log_vars, pp_dirs, batchfname, makemapsname,
                profiles, compilation_steps)
    
    # End timer and report
    endtime = time.time()
    elapsed = endtime-starttime
    print("\nTime elapsed to run was",elapsed,"seconds\n")
