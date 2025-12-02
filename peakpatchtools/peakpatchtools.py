# peakpatchtools.py

####################################################################################################
####################################################################################################
#####                                                                                          #####
#####                                                                                          #####
#####   `7MM"""Mq.                  `7MM          `7MM"""Mq.         mm         `7MM           #####
#####     MM   `MM.                   MM            MM   `MM.        MM           MM           #####
#####     MM   ,M9 .gP"Ya   ,6"Yb.    MM  ,MP'      MM   ,M9 ,6"Yb.mmMMmm ,p6"bo  MMpMMMb.     #####
#####     MMmmdM9 ,M'   Yb 8)   MM    MM ;Y         MMmmdM9 8)   MM  MM  6M'  OO  MM    MM     #####
#####     MM      8M""""""  ,pm9MM    MM;Mm         MM       ,pm9MM  MM  8M       MM    MM     #####
#####     MM      YM.    , 8M   MM    MM `Mb.       MM      8M   MM  MM  YM.    , MM    MM     #####
#####   .JMML.     `Mbmmd' `Moo9^Yo..JMML. YA.    .JMML.    `Moo9^Yo.`MbmoYMbmd'.JMML  JMML.   #####
#####                                                                                          #####
#####                                                                                          #####
#####   MMP""MM""YMM              `7MM                                                         #####
#####   P'   MM   `7                MM                                                         #####
#####        MM  ,pW"Wq.   ,pW"Wq.  MM  ,pP"Ybd                                                #####
#####        MM 6W'   `Wb 6W'   `Wb MM  8I   `"                                                #####
#####        MM 8M     M8 8M     M8 MM  `YMMMa.                                                #####
#####        MM YA.   ,A9 YA.   ,A9 MM  L.   I8                                                #####
#####      .JMML.`Ybmd9'   `Ybmd9'.JMML.M9mmmP'                                                #####
#####                                                                                          #####
#####                                                                                          #####
####################################################################################################
####################################################################################################
#                                                                                                  #
# Script for post-processing and plotting Peak Patch and WebSky data.                              #
#                                                                                                  #
# The programme is structured around objects of class PeakPatch. All the variables in Peak Patch   #
# parameter files are attributes of PeakPatch objects. For a full list, see the example in         #
# <...>/peakpatch/example/param/param.params. The following functions are attributes as well:      #
#                                                                                                  #
#     add_halos() : load Peak Patch halo catalogue                                                 #
#                                                                                                  #
#     add_field() : load Peak Patch initial conditions fields                                      #
#                                                                                                  #
# See function definitions for full descriptions of usage and specifications of all keyword        #
# arguments.                                                                                       #
#                                                                                                  #
####################################################################################################
#                                                                                                  #
# USAGE: for a Peak Patch run in directory /path/to/peak-patch-run/, run python                    #
#                                                                                                  #
# >>> import os, sys, numpy as np                                                                  #
# >>> sys.path.insert( 0 , '<...>/peakpatch/python' )                                              #
# >>> from peakpatchtools import PeakPatch                                                         #
# >>> run = PeakPatch('/path/to/peak-patch-run/')                                                  #
#                                                                                                  #
####################################################################################################
#                                                                                                  #
# How We Discretise a Volume of Space in the Peak Patch Calculation                                #
#                                                                                                  #
# Peak Patch tries finds the cosmological structure that will form in a cubic region of space of   #
# volume (boxsize)^3 Mpc. The calculation discretises this volume into a cubic lattice of n_eff^3  #
# voxels.                                                                                          #
#     Results from early work on the evolution of structure has shown that on the largest scales,  #
# dark matter halos around galaxies are mostly only affected gravitationally by their nearest      #
# neighbours, meaning that a very large cubic volume of space can be divided into pieces (and the  #
# Peak Patch calculation run in parallel across them) provided a buffer is used. This buffer has   #
# to be at least about 32 Mpc in thickness, so the Peak Patch calculation has to have a shell of   #
# buffer voxels with a thickness n_buff voxels such that the thickness in Mpc r_buff is greater    #
# than or equal to about 32 Mpc.                                                                   #
#         r_buff = n_buff * r_latt > 32 Mpc                                                        #
# When a run is broken up into cubic blocks (we call these cubic sub-lattices "tiles") each tile   #
# must also have a buffer of thickness n_buff voxels, consisting of voxels from the adjacent tile. #
#     The total number of voxels is n^3, where the number when we subract off the exterior buffers #
# is n_eff^3. The number of voxels in a single tile (including buffers) is called n_mesh^3, and    #
# the number of cubic tiles is n_tile^3. Mathematically, these are all related:                    #
#         n = n_eff + 2*n_buff = n_tile * (n_mesh-2*n_buff) + 2*n_buff                             #
#     Whereas we define the discretisation as the total number of voxels n^3 with buffers included #
# (since the discretisation characterises the computational expense of the PEak Patch calculation, #
# it makes sense to include buffers in this number) the total physical volume is boxsize^3, which  #
# excludes the buffers (since this is the volume of what is physically meaningful in the Peak      #
# Patch calculation, it wouldn't make sense to include the buffers, which are tossed out once we   #
# have found halos). From this we can then define a lattice spacing, the physical volume (in Mpc)  #
# of a voxel r_latt^3 defined as follows:                                                          #
#         boxsize = r_latt * n_eff = r_latt * ( n - 2*n_buff )                                     #
# Figure 1 shows a schematic representation of all of these various length scales.                 #
#                                                                                                  #
#         |  |<---------------- n (voxels) ----------------->|                                     #
#         v          |<------- boxsize (Mpc) ------->|                                             #
#         -   --- --- --- --- --- --- --- --- --- --- --- ---   -                                  #
#   r_latt   | b | b | b | b | b | b | b | b | b | b | b | b |  ^                                  #
#   (Mpc) -   --- --- --- --- --- --- --- --- --- --- --- ---   |                                  #
#         ^  | b | b | b | b | b | b | b | b | b | b | b | b |  |                                  #
#  ------ |   --- --- --- --- --- --- --- --- --- --- --- ---   |  ---------                       #
#  ^         | b | b | 1 | 1 | 1 | 1 |2b1|2b1|   |   | b | b |  |          ^                       #
#  |          --- --- --- --- --- --- --- --- --- --- --- ---   |          |                       #
#  |         | b | b | 1 | 1 | 1 | 1 |2b1|2b1|   |   | b | b |  |          |                       #
#  |          --- --- --- --- --- --- --- --- --- --- --- ---   |          |                       #
#  |         | b | b | 1 | 1 | 1 | 1 |2b1|2b1|   |   | b | b |  |          |                       #
#  |          --- --- --- --- --- --- --- --- --- --- --- ---   |          |                       #
# n_eff      | b | b | 1 | 1 | 1 | 1 |2b1|2b1|   |   | b | b |  n (voxels) |                       #
# (voxels)    --- --- --- --- --- --- --- --- --- --- --- ---   |          |                       #
#  |         | b | b |3b1|3b1|3b1|3b1|4b1|4b1|   |   | b | b |  |       boxsize (Mpc)              #
#  |          --- --- --- --- --- --- --- --- --- --- --- ---   |          |                       #
#  |         | b | b |3b1|3b1|3b1|3b1|4b1|4b1|   |   | b | b |  |          |                       #
#  |          --- --- --- --- --- --- --- --- --- --- --- ---   |          |                       #
#  |         | b | b |   |   |   |   |   |   |   |   | b | b |  |          |                       #
#  |          --- --- --- --- --- --- --- --- --- --- --- ---   |          |                       #
#  v         | b | b |   |   |   |   |   |   |   |   | b | b |  |          v                       #
#  ---------  --- --- --- --- --- --- --- --- --- --- --- ---   |  ---------                       #
#            | b | b | b | b | b | b | b | b | b | b | b | b |  |                                  #
#             --- --- --- --- --- --- --- --- --- --- --- ---   |                                  #
#            | b | b | b | b | b | b | b | b | b | b | b | b |  v                                  #
#             --- --- --- --- --- --- --- --- --- --- --- ---   -                                  #
#              ->|   |<- r_latt (Mpc)              ->|       |<- n_buff (voxels)                   #
#            |<------ n_mesh (voxels) ------>|     ->|       |<- r_buff (Mpc)                      #
#                                                                                                  #
#     Legend:                                                                                      #
#     b   - cells in the buffer                                                                    #
#     1   - cells in tile 1                                                                        #
#     ib1 - cells in tile i that are also in the buffer for tile 1, note that some cells           #
#           labelled "1" are also in the buffer for the adjaceent tiles, but I've just             #
#           labelled them 1 for simplicity here.                                                   #
#                                                                                                  #
# Figure 1: For a run with n_tile=2, n_buff=2, n=12, boxsize=8.0 Mpc and therefore r_buff=1 Mpc,   #
# this figure shows a cross section of the Peak Patch discretisation for the field. Note that n    #
# and next are equal.                                                                              #
#                                                                                                  #
####################################################################################################
#                                                                                                  #
# Selected References                                                                              #
#                                                                                                  #
# The Peak Patch and WebSky Papers:                                                                #
#                                                                                                  #
#     1) Peak Patch was introduced in a series of three papers: "ads:1996ApJS..103....1B",         #
#        "ads:1996ApJS..103...41B" and "ads:1996ApJS..103...63B".                                  #
#                                                                                                  #
#     2) Peak Patch was updated for massively parallel computing in arXiv:1810.07727.              #
#                                                                                                  #
#     3) The WebSky Simulations were introduced in arXiv:2001.08787                                #
#                                                                                                  #
#     4) The groundwork for Peak Patch is laid out in BBKS (named for its authors, Bardeen, Bond,  #
#        Kaiser and Szalay) which lays out the theory of Gaussian random field statistics in       #
#        cosmology: "ads:1986ApJ...304...15B". "Excursion sets" which are the precursor to Peak    #
#        Patch's homogeneous ellipsoidal collapse halo finding method, are introduced in a paper   #
#        BCEK (again, for authors Bond, Cold, Efstathiou and Kaiser) "ads:1991ApJ...379..440B".    #
#                                                                                                  #
#     5) The WebSky simulations approximate halos as spherically summetric with a has pressure     #
#        profile similar to the NFW profile usually called a BBPS profile (again, this is named    #
#        for its authors: Battaglia, Bond, Pfrommer and Sievers) which were introduced in a series #
#        of four papers: arXiv:1109.3709, arXiv:1109.3711, arXiv:1209.4082 and arXiv:1405.3346.    #
#                                                                                                  #
# Further reading:                                                                                 #
#                                                                                                  #
#      - There's some additional information about the anisotropic information inherent that can   #
#        be extracted from Peak Patch halo catalogues in paper arXiv:2101.01455.                   #
#                                                                                                  #
#      - Peak Patch was compared to a number of other simulations for generating realisations of   #
#        large-scale dark matter distributions without running full N-body or hydrodynamical       #
#        simulations in a series of papers: arXiv:1806.09477, arXiv:1806.09497, arXiv:1806.09499.  #
#        Among those compared, this paper considers HALOGEN/Patchy, ICE-COLA, simple log-normal    #
#        distribution models, Pinocchio, and a standard 2LPT model. Peak Patch performs very       #
#        favourably.                                                                               #
#                                                                                                  #
#      - For the most up-to date line intensity mapping work work, see Patrick Horlaville's paper
#        arXiv:2309.15733 and check out the GitHub 
#      - For a general review of cosmology that references the Peak Patch picture for structure    #
#        formation, see the review articles:                                                       #
#         o  https://www.cita.utoronto.ca/~bond/cweb/wb1_fulltext.pdf                              #
#         o  https://www.cita.utoronto.ca/~bond/cweb/wb2_fulltext.pdf                              #
#                                                                                                  #
# Title graphic adapted from http://www.patorjk.com/software/taag/.                                #
#                                                                                                  #
####################################################################################################
# created: Nathan J. Carlson,  9 May 2023
# edited:  Nathan J. Carlson,  3 September 2024

# Import standard python modules
import os, subprocess, datetime, time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams['text.usetex']=True

# Import astropy and astropy
# Non-standard python modules like astropy can be accessed on Niagara in a virtual environment (see
# https://docs.scinet.utoronto.ca/index.php/Installing_your_own_Python_Modules for instructions).
# This is a hassle though, so if you just want to do something quick like make a copy of a run
# directory with the PeakPatch.copy() method, if astropy can't be found, a warning is thrown instead
# of an error
import warnings

# This is needed to run on nix. If you don't have nix, ignore this
nix = os.getenv('NIX_BUILD')
if nix:
    nix = int(nix)

try:
    machine = subprocess.check_output('hostname').decode('utf-8')[:-1]
    if machine[:3]+machine[-13:] == 'nia.scinet.local':
        from astropy.config import paths
        import tempfile
        # Create a temporary directory to store config
        temp_config_dir = tempfile.mkdtemp()
        paths.set_temp_config(temp_config_dir)

    from astropy import units as u
    from astropy import constants as const
    from astropy.cosmology import Planck18, z_at_value
    if nix != 1:
        import healpy as hp
except:
    host = os.getenv('HOSTNAME')
    if 'nia-login' in host and '.scinet.local' in host:
        warnings.warn('Could not import astropy, attempting to run PeakPatch without may resul'+
            't in errors.\nIt looks like you are in a Niagara login node, astropy is a non-sta'+
            'ndard module for Niagara, to use\nit you need to pip install it in a virtual envi'+
            'ronment.\nFor more details, see https://docs.scinet.utoronto.ca/index.php/Install'+
            'ing_your_own_Python_Modules',ImportWarning)
    elif host=='nia-dm1' or host=='nia-dm2':
        warnings.warn('Could not import astropy, attempting to run PeakPatch without may result in'+
                'errors.\nIt looks like you are in a Niagara data mover node, PeakPatch is better '+
                'suited for use on the login\nnodes.',ImportWarning)
    elif 'nia' in host and '.scinet.local' in host:
        warnings.warn('Could not import astropy, attempting to run PeakPatch without may result in'+
                'errors.\nIt looks like you are in a Niagara compute node, PeakPatch is better sui'+
                'ted for use on the login\nnodes.',ImportWarning)
    else: warnings.warn('Could not import astropy, attempting to run PeakPatch without may result '+
                'in errors.',ImportWarning)

# Try to locate the Peak Patch directory
try:
    peak_patch_dir = os.getenv('PP_DIR')
except:
    warnings.warn('Could not find the Peak Patch directory in your environment variables. This cou'+
            'ld cause problem\nwith some of peakpatchtools functions. You can fix this in terminal'+
            'by finding the Peak Patch\ndirectory `<...>/peakpatch/` and running the following com'+
            'mand in terminal\n\n    export PP_DIR=<~/...>/peak-patch\n\nIn the meantime, we\'ll t'+
            'ry running without defining the Peak Patch directory.',Warning)
    peak_patch_dir = None

# Class for post-processing on Peak Patch runs.
class PeakPatch:

    # Initializes PeakPatch object by reading in parameter file
    def __init__(self, run_dir=None):

        # If no run directory passed, use the example one
        if run_dir==None:
            run_dir = peak_patch_dir+'/example'

        # Strip any leading or trailing spaces from run_dir passed as argument of __init__()
        else:
            run_dir = run_dir.strip() 

        params_file_params = run_dir + '/param/param.params'
        params_file_ini = run_dir + '/param/parameters.ini'
        params_file = params_file_params

        # Check that run_dir is formatted as a Peak Patch run
        if not os.path.exists(params_file):
            print(f'{params_file} not found... Checking {params_file_ini}')
            params_file = params_file_ini
            if not os.path.exists(params_file):
                raise OSError( ('files not found: {0} or {1}.\n{2} must be formatted as a P'+
                                'eak Patch run directory').format(params_file_params,
                        params_file_ini, run_dir) )
        self.run_dir = os.path.abspath(run_dir) # Make sure run_dir is absolute not relative path
        self.run_dir = self.run_dir.rstrip('/') # Strip trailing / from run_dir

        # Execute parameter file parameters as variables of the class object
        with open(params_file) as f:
            line_counter = 0
            for line in f.readlines():
                line          = line.strip() # Strip white space
                line_counter += 1

                # Ignore lines that are just formatted as python comments or are blank
                if len(line)==0 or line[0]=='#' or line[0]=='[':
                    continue

                # Split the line at each '#', then just take the first of these to drop trailing
                # comments and strip white space again, then split at the '='
                line = line.split('#')[0].strip().split('=')

                if len(line)==2:

                    # strip all whitespace from lines
                    line[0] , line[1] = line[0].strip() , line[1].strip()
                    
                    # Some variable definitions in the parameter file are functions of earlier
                    # variables, so we must evaluate the line as code. The value `line[1]` of each
                    # parameter `line[0]` is saved as an attribute of the class PeakPatch
                    exec( '{0}={1};self.{0}={0}'.format( line[0] , line[1] ) )
                    # It would be preferable to remove the exec statements, but this would require
                    # parsing out each line[1], some of which are expressions, e.g. with import re

                # Raise a warning and skip if line cannot be parsed
                else:
                    warnings.warn( justify_text( ('Line {0} of the parameter file could not be par'+
                                'sed and is being skipped. Prameter files with codec .params shoul'+
                                'd be formatted as python code, make sure there are no python erro'+
                                'rs in this line. Additionally, be sure that if the line is specif'+
                                'ying the value of a string that that string does not contain any '+
                                'illegal characters (e.g. "#", "=").').format(line_counter),
                            width=100 , first_line_buffer='AttributeError' ) , Warning )
 
        # Clean up after all those exec statements by deleting superfluous global variables
        for key in vars(self): exec( 'del( {0} )'.format(key) )

        # Effective voxel dimension when buffers are excluded
        self.neff = self.next - 2*self.nbuff

        # Fundamental and Nyquist wavenumbers
        self.kfund = 2*np.pi/self.boxsize
        self.knyq  = 2*np.pi/self.boxsize*(self.next/2+1)



    def __str__(self):
        return os.path.basename(self.run_dir)



    def __repr__(self):
        output = ('PeakPatch object representing {0}\nwith volume ({1} Mpc)^3, resolution {2}^3'
                 .format(self.run_dir,self.boxsize,self.neff) )
        if hasattr(self,'N_halos'): output+=', {0} halos'.format(self.N_halos)
        field_count=0
        if hasattr(self, 'rhog' ): field_count+=1
        if hasattr(self, 'rho'  ): field_count+=1
        if hasattr(self, 'zetag'): field_count+=1
        if hasattr(self, 'zeta' ): field_count+=1
        if hasattr(self, 'chi'  ): field_count+=1
        if field_count>0: output+=', {0} fields'.format(field_count)
        return output



    def run(self, cmdln=''):
        # This function allows you to quickly run Peak Patch from a paramter file. You must first
        # create a PeakPatch object, then you can run Peak Patch. This function makes use of the
        # setup.sh shell script in the Peak Patch home directory.
        #    
        # Usage:
        # >>> run = Peak_atch( <run> )
        # >>> run.run()
        #    
        # You can also make use of the debug run feature or hold feature (respectively running a
        # small job in the debug queue so that it finishes more quickly, and setting the run up but
        # not submitting it to the job scheduler) with `run.run('debug')` and `run.run('hold')`.

        # Check that machine is Niagara and command line arguments are properly formatted

        if self.machine == 'niagara':
            pkpdir = os.environ['PP_DIR']
            if cmdln.strip().lower() not in [ '' , 'debug' , 'hold' ]:
                raise ValueError('Command line arguments for peakpatch/setup.sh are none, `debug` o'
                                 +'r `hold`.')
        else:
            raise OSError('I haven\'t set this up to work on machines other than Niagara yet.')

        # Run shell script to pass Peak Patch run to the job scheduler
        os.system( 'cd {0};./setup.sh {1} {2}'.format(pkpdir,self.run_dir,cmdln) )



    def autofill( self , autos=True , **kwargs ):
        # This function calculates values for all the parameters in a Peak Patch parameter file that
        # have dependencies on other Peak Patch prameters (such as `ntasks`, the number of parallel
        # computational tasks that a Peak Patch run is broken up into, which depends on `nnodes` the
        # number of parallel computing nodes requested, and `tpnode` the number of tasks run on each
        # compute node).
        # 
        # Optional parameters:
        # 
        #     autos : bool
        #         If `autos` is False then the function just outputs the list of the parameters with
        #         dependences. If `autos` is True, PeakPatch attributes among the parameters whose
        #         values would have been returned are updated with appropriate vlues based on the
        #         other attributes.
        # 
        #     **kwargs : key word arguments (passed as variables treated like dict)
        #         There are additionally some Peak Patch parameters whose values have standard 
        #         formatting. These can be appropriately reformatted by passing them in as kwargs.
        #         Additionally, if the only kwarg is passed in as
        #             defaults=True
        #         then all parameters with defalut forms are reformatted.

        # Peak Patch parameters with dependencies whose values must be calculated from other params
        autofill_params  = [ 'ntasks' , 'ntasks_map' , 'densfilein' , 'densfileout' , 'iZeld' ,
                             'ntilemerge' , 'ntasksmerge' , 'nsub' , 'next' , 'dcore_box' ,
                             'cellsize' , 'buffersize' , 'dL_box' , 'mlatt' ,
                             'nlx' , 'nly' , 'nlz' , 'n1' , 'n2' , 'n3' ]

        # Peak Patch post-processing convenience parameters
        post_proc_params = [ 'neff' , 'kfund' , 'knyq' ]

        # Peak Patch parameters with default values
        if kwargs == {}:
            default_params = []
        else:
            default_params = [ key for key in kwargs ]
            for param in default_params:
                if param not in [ 'run_name' , 'short_name' , 'Omvac' ,
                                  'TabInterpZ1' , 'TabInterpZ2' ]:
                    raise KeyError(param+' must be a Peak Patch parameter.')

        # If autos is False, just return names of autofill parameters
        if autos == False: return autofill_params + post_proc_params + default_params

        # Calculate Peak Patch parameters with dependencies on other params
        else:
            self.ntasks      = self.nnodes * self.tpnode
            self.ntasks_map  = self.ppn_map * self.nnodes_map
            self.densfilein  = self.run_name
            self.densfileout = self.run_name
            self.iZeld       = self.ilpt
            self.ntilemerge  = self.ntile
            self.ntasksmerge = self.ntasks
            self.nsub        = self.nmesh - 2 * self.nbuff
            self.next        = self.nsub * self.ntile + 2 * self.nbuff
            self.dcore_box   = self.boxsize / self.ntile
            self.cellsize    = self.dcore_box / self.nsub
            self.buffersize  = self.cellsize * self.nbuff
            self.dL_box      = self.dcore_box + 2 * self.buffersize
            self.mlatt       = 2.7754e11 * (self.Omx+self.OmB) * self.h**2 * self.cellsize**3
            self.nlx         = self.ntile
            self.nly         = self.ntile
            self.nlz         = self.ntile
            self.n1          = self.nmesh
            self.n2          = self.nmesh
            self.n3          = self.nmesh

            # Variables only in post-processing
            self.neff  = self.next - 2*self.nbuff
            self.kfund = 2*np.pi / self.boxsize
            self.knyq  = 2*np.pi / self.boxsize*(self.next/2+1)

            # Variables with default formatting
            if 'run_name' in default_params   : self.run_name    = '{0}Mpc_n{1}_nb{2}_nt{3}'.format(
                                                                   self.boxsize, self.nmesh, 
                                                                   self.nbuff, self.ntile )
            if 'short_name' in default_params : self.short_name  = '{0}Mpc_nb{1}'.format( 
                                                                   self.boxsize, self.nbuff )
            if 'Omvac' in default_params      : self.Omvac       = 1 - self.Omx - self.OmB
            if 'TabInterpZ1' in default_params: self.TabInterpZ1 = -1 + 1e-4
            if 'TabInterpZ2' in default_params: self.TabInterpZ2 =  1 - 1e-4



    def auto_set_HEC_table(self):

            # Default values
            self.TabInterpNx = 50
            self.TabInterpNy = 20
            self.TabInterpNz = 20
            self.TabInterpX1 =  1.5
            #self.TabInterpX2 =  8.0
            self.TabInterpY1 =  0.0
            self.TabInterpY2 =  0.5
            self.TabInterpZ1 = -1 + 1e-4
            self.TabInterpZ2 =  1 - 1e-4

            # Set max overdensity parameter as a function of box maximum redshift
            r_comoving_max = np.sqrt(   (self.cenx+self.boxsize/2)**2 
                                      + (self.ceny+self.boxsize/2)**2 
                                      + (self.cenz+self.boxsize/2)**2 )

            tab = z_of_r_comoving_table( zmin=1.0e-3, zmax=3402.0, dlog10z=6.538e-3, Omega_r0=0.0,
                    Omega_m0=self.Omx+self.OmB, Omega_k0=0.0, Omega_Lambda=self.Omvac, 
                    H_0=self.h*100, c=299792458.0 )
            z = z_of_r_comoving( r_comoving_max, z_chi_table=tab )

            self.TabInterpX2 = HEC_Frho_of_z( z, Omega_m=self.OmB+self.Omx, Omega_k=0.0,
                                              Omega_Lambda=self.Omvac )

            # Update number of overdensity values in the HEC table
            while(self.TabInterpX2/self.TabInterpX1)**(1/(1-self.TabInterpNx)) > (8/1.5)**(1/49):
                self.TabInterpNx += 1



    def edit(self, autofill_params=True, **kwargs):
        # This function allows you to quickly make edits to an existing Peak Patch run. To make
        # adjustments, you can pass kwargs, a dictionary of any variables you wish to change.
        # 
        # Usage:
        # >>> run  = PeakPatch( <run> )
        # >>> edit = run.edit( seed=23456 , boxsize=1000.0 )
        # 
        # Parameters
        # 
        #     autofill_params : bool
        #         A number of parameters in Peak Patch parameter files are typically dependent on
        #         other parameters. (For instance, the number of parallel tasks in the Peak Patch
        #         calculation, `ntasks' is set to the product of the number of nodes and the number
        #         of tasks per node.) If autofill_params is set to True, such parameters will be set
        #         as they are in the example parameter file. Default is True.
        # 
        #     kwargs : dictionary
        #         Any changes from the run you're copying should be specified in this dictionary.
        #         For instance, if you want to change the seed to 23456 and the boxsize to 1000 Mpc,
        #         but keep all other parameters the same, kwargs would be set to
        #             kwargs = { 'seed' : 23456 , 'boxsize' : 1000.0 }
        #         If the dictionary is empty, no changes are made. Default is empty dictionary: {}.
        #
        # Returns:
        #
        #     PeakPatch(copy_dir) : PeakPatch object
        #         The copied/edited run is then read in as a PeakPatch class object and returned.

        # Check that kwargs is formatted correctly
        if type(kwargs) != dict: raise TypeError('kwargs must be of type dict.')
        if kwargs=={}:
            warnings.warn('You\'ve made no edits. Why are you calling the edit function?',Warning)
            return
        for key in kwargs:
            if not hasattr(self,key):
                raise KeyError( ('key {0} not recognised. Keys must be standard Peak Patch variabl'+
                                 'es.').format(key) )
            if key in self.autofill(False):
                raise KeyError( ('You can\'t arbitrarily reset the following parameters:\n{0}\nas '+
                                 'these are dependent on other parameters. See\n{0}/example/param/'+
                                 'param.params\nfor details.')
                                 .format( self.autofill(False) , peak_patch_dir ) )

            # Edit self
            else: vars(self)[key] = kwargs[key]

        # Update parameters with dependencies
        self.autofill()

        # Read parameter file and identify all lines that need to be edited
        edits = []
        with open(self.run_dir+'/param/param.params') as f:
            for line in f.readlines():
                l = line.strip() # Strip white space

                # Ignore lines that are just comments or blank
                if len(l)==0 or l[0]=='#':
                    continue

                # Split the line at each '#', then just take the first of these to drop trailing
                # comments and strip white space again, then split at the '='
                l = l.split('#')[0].strip().split('=')

                # Keep only lines of the form `key = value`
                if len(l)==2:

                    # strip all whitespace from lines
                    l[0] , l[1] = l[0].strip() , l[1].strip()

                    # Check if the parameter is one that has been edited
                    if ( ( l[0] in kwargs or l[0] in self.autofill(False) )
                            and l[1] != str(vars(self)[l[0]]) ):

                        # Replace the line in the parameter file
                        if type(vars(self)[l[0]]) == str:
                            edits += [[ line.split('#')[0].strip()                ,
                                        "{0}='{1}'".format(l[0],vars(self)[l[0]]) ]]
                        else:
                            edits += [[ line.split('#')[0].strip()                ,
                                        '{0}={1}'  .format(l[0],vars(self)[l[0]]) ]]

        # Edit the necessary lines that have changed
        for j in range(len(edits)):

            # Its generally bad practice to include math in your parameter files, but for backwards
            # compatibility, we've kept it in as a feature. This leads to problems if you have
            # exponentiation though because "**" has special meaning for the sed command. This block
            # of code resolves that issue
            if '**' in edits[j][0]:
                temp = edits[j][0].split('**')
                edits[j][0] = edits[j][0].split('**')[0]
                for k in range(1,len(temp)):
                    edits[j][0] += '\\*\\*' + temp[k]

            # Note, Peak Patch parameters are not allowed to have '#' in their values. 
            if '#' in edits[j][0] or '#' in edits[j][0]:
                raise ValueError( justify_text( ( 'PeakPatch parameter files do not support parame'+
                                                  'ters or values containing number signs ("#").' ),
                                                 width=100, first_line_buffer='ValueError' ) )

            # Make edits
            os.system( ("cd {0}/param;sed \"s#{1}#{2}#g\" param.params > temp;mv temp param.params"
                       ).format( self.run_dir , edits[j][0] , edits[j][1] ) )

 

    def cleanup(self, flag=None):
        # This function mimics the cleanup.sh shell script in the Peak Patch home directory. Like
        # that shell script, it's purpose is for deleting unneccessary data or runs that didn't work
        # so that they can be re-run.
        # 
        # Optoinal parameter:
        # 
        #     flag : string or None
        #         An optional parameter to clean up the run according to different set of rules. The
        #         allowed values are:
        #         None - the default, delete everything except param/param.params
        #         archive - move important run output to a direcotry archive/<subdirectory> (where
        #             <subdirectory> indicates the time and date that the run finished) all of the
        #             runs subdirectories besides param and archive are deleted.
        #         maps - if there is a maps subdirecotory, this flag will return the run directory
        #             to the state just after Peak Patch is run but before pks2map is run, so every-
        #             thing except shell scripts are deleted from the maps subdirectory.
        #         fail - ideal for debugging, this deletes everything in the run directory except
        #             logfiles and param/param.params, useful if you're editing the source code and
        #             introduce a bug that is causing IC fields or halo catalogues not to be
        #             calculated properly.
        # 
        # Note that it is best practices NOT to use the cleanup() function with argument 
        # flag='archive' in conjunction with the edit() function, as the latter will change the
        # parameter file and the archive function is meant purely for running the same run multiple
        # times from the same parameter file (e.g. for comparing source code edits).

        # If no flag passed, just delete everything except parameter file
        if flag == None:
            os.system( 'cd {0};rm -rf !("param");cd param;rm -rf !("param.params")'.format(
                    self.run_dir) )

        # If flag is set to archive, archive old run data in an archive subdirectory and clear
        # everything else
        elif flag == 'archive':
            
            # Read metadate of output directory to figure out when it was created
            output_time     = os.path.getctime( '{0}/output/'.format( self.run_dir) )

            # Reformat time that output was created into a string in the form `formatted_time' =
            # <hours>:<minutes>:<seconds>_<day>-<month>-<year>
            output_time_str = datetime.fromtimestamp( output_time )
            formatted_time  = output_time_str.strftime( '%H:%M:%S_%d-%m-%Y' )

            # Move all run relevant run information to a new subdirecotry `formatted_time` and
            # delete everything else
            os.system( ( 'cd {0};rm -rf !("param"|"output"|"maps"|"logfiles"|"tables"|"archive");m'+
                         'kdir -p archive/{1};mv !("param") arcive/{1}/.;cd param;rm -rf !("param.'+
                         'params")' ).format( self.run_dir, formatted_time ) )

        # Just clear the maps directory, reset so that the run directory looks as it did just after
        # running Peak Patch but before running pks2map
        elif flag == 'maps':
            
            # Check that there is a maps directory, if there is, clear everything but shell scripts
            if os.path.isdir( '{0}/maps/'.format(self.run_dir) ):
                os.system( 'cd {0}/maps;rm -rf !(*.sh)'.format(self.run_dir) )

            # If there is no maps direcotry, throw an error
            else:
                raise SystemError('The directory {0}/maps/ does not exist.'.format(self.run_dir))
            
        # If flag is set to fail, remove everything except for logfiles, useful for a run that
        # didn't work properly e.g. if you get bad halo catalogues or IC fields with NaNs.
        elif flag == 'fail':
            os.system( 'cd{0};rm -rf !("param"|"logfiles")'.format(self.run_dir) )

        # If the flag is not recognized, throw an error
        else:
            raise ValueError('Flag {0} not recognized, aborting.'.format(flag))


 
    def copy(self, copy_dir=None, autofill_params=None, **kwargs):
        # This function allows you to quickly make a copy of a Peak Patch run. To make adjustments 
        # to the copy, you can pass the optional argument kwargs, a dictionary of any variables you
        # wish to change.
        # 
        # Usage:
        # >>> run  = PeakPatch( <run> )
        # >>> copy = run.copy()
        # >>> edit = run.copy( seed=23456 , boxsize=1000.0 )
        # 
        # Optional parameters:
        # 
        #     copy_dir : string or None
        #         The directory in which to save the copy of your run. If None, it will create a
        #         copy in a directory called self.run_dir+'_<#>' where <#> is the lowest nonzero
        #         positive integer for which a directory does not already exist.
        # 
        #     autofill_params : bool
        #         A number of parameters in Peak Patch parameter files are typically dependent on
        #         other parameters. (For instance, the number of parallel tasks in the Peak Patch
        #         calculation, `ntasks' is set to the product of the number of nodes and the number
        #         of tasks per node.) If autofill_params is set to True, such parameters will be set
        #         as they are in the example parameter file. Default is False if no kwargs are
        #         passed or True otherwise.
        # 
        #     kwargs : dictionary
        #         Any changes from the run you're copying should be specified in this dictionary.
        #         For instance, if you want to change the seed to 23456 and the boxsize to 1000 Mpc,
        #         but keep all other parameters the same, kwargs would be set to
        #             kwargs = { 'seed' : 23456 , 'boxsize' : 1000.0 }
        #         If the dictionary is empty, no changes are made. Default is empty dictionary: {}.
        # 

        # Check that kwargs is a dictionary
        if type(kwargs) != dict: raise TypeError('kwargs must be of type dict.')

        # If no copy_dir is given, make one up
        if copy_dir == None:
            i=1
            while os.path.exists( self.run_dir+'_'+str(i) ): i+=1
            copy_dir = self.run_dir+'_'+str(i)

        # Otherwise, copy_dir must be a string
        elif type(copy_dir) != str:
            raise TypeError('copy_dir must be string or None.')

        # If target directory for the copied run exists, it is overwritten
        os.system( 'rm -rf {0}'.format( copy_dir ) )

        # Make the copy direcotry and it's parameter file and load it as a PeakPatch object
        os.system( 'mkdir -p '+copy_dir+'/param/' )
        os.system( 'cd {0}/param/;cp {1}{2} .'.format(copy_dir,self.run_dir,'/param/param.params') )
        copys_self = PeakPatch(copy_dir)

        # Make any edits specified in the key word arguments
        if len(kwargs) > 0:
            if autofill_params == None or autofill_params == True:
                copys_self.edit(**kwargs)
            else:
                copys_self.edit(autofill_params=False,**kwargs)

        # Autofill if autofill_params set to True
        else:
            if autofill_params == True:
                copys_self.autofill()

        # Return the copied run as PeakPatch object
        return copys_self



    def diff(self, other, dirs=False):
        # This function allows you to quickly see what the difference is between two runs's
        # parameter files.
        # 
        # Usage:
        # >>> run1 = PeakPatch( <run1> )
        # >>> run2 = PeakPatch( <run2> )
        # >>> run1.diff(run2)
        # 
        # Parameters:
        #     
        #     other : PeakPatch object
        #     The second Peak Patch run to compare with.
        #     
        #     dirs : bool
        #     If True it compares attribute `run_dir', the directory where the run is stored. If
        #     False then `run_dir' is ignored. Usually you just want to know if the cosmologies
        #     differ so dirs is set to False by default.

        def greaterlen(maxlens,lens):
            # input must be two lists of 3 integers
            if lens[0] > maxlens[0]: maxlens[0] = lens[0]
            if lens[1] > maxlens[1]: maxlens[1] = lens[1]
            if lens[2] > maxlens[2]: maxlens[2] = lens[2]
            return maxlens

        # Find parameters that differ between self and other and save them in the dict diffs
        diffs = {}
        maxlens = [5,len(str(self)),len(str(other))]
        for key in vars(self):
            
            # Don't count the run directories as differeing if dirs is set to False
            if key!='run_dir' or (key=='run_dir' and dirs==True):

                # Find any parameters whose values are not the same in self and other
                if key in vars(other):
                    if ( np.array(vars(self)[key]) == np.array(vars(other)[key]) ).all() == False:
                        if ( ( np.array(vars(self)[key]).shape == () and
                               type(vars(self)[key]) not in (list,tuple,dict)
                             ) or len(str(vars(self)[key]))+len(str(vars(other)[key])) < 50 ):
                            diffs[key] = [ vars(self)[key] , vars(other)[key] ]
                            maxlens = greaterlen( maxlens, [ len(key)                   ,
                                                             len(str(vars( self)[key])) ,
                                                             len(str(vars(other)[key])) ] )
                        else:
                            diffs[key] = [ 'x' , 'y ≠ x' ]
                            maxlens = greaterlen( maxlens, [len(key),5] )

                # Find any parameters that are in self but not in other
                else:
                    diffs[key] = [ vars(self)[key] , None ]
                    maxlens = greaterlen( maxlens, [ len(key),len(str(vars(self)[key])),0 ] )

        # Find any parameters that are in other but not in self
        for key in vars(other):
            if key not in vars(self) and ( key!='run_dir' or (key=='run_dir' and dirs==True) ):
                diffs[key] = [ None , vars(other)[key] ]
                maxlens = greaterlen( maxlens, [ len(key),0,len(str(vars(other)[key])) ] )

        # Nicely display diffs
        if diffs == {}:
            return 'Identical parameter files'
        else:
            headers=['     ',str(self),str(other)]
            print( f'{headers[0]:<{maxlens[0]}}  |  {headers[1]:<{maxlens[1]}}  |  {headers[2]}' )
            print( '-'*(2+maxlens[0]) + '+' + '-'*(4+maxlens[1]) + '+' + '-'*(2+maxlens[2]) )
            for key in diffs:
                if key != 'run_dir' or ( key == 'run_dir' and dirs == True ): print(
                    f'{key:<{maxlens[0]}}  |  {diffs[key][0]:<{maxlens[1]}}  |  {diffs[key][1]}' )

        return diffs



    ################################################################################################
    ###                                                                                          ###
    ###                            LOADING Peak Patch AND WebSky DATA                            ###
    ###                                                                                          ###
    ################################################################################################
    # The scripts in this scetion load Peak-Patch/WebSky data products (including initial conditions
    # cosmic fields, dark matter halo catalogues generated by Peak Patch and sky maps generated by
    # WebSky) into objects of class PeakPatch. Options exist to only read subsets of the data for
    # the purposes of economizing RAM when manipulating the large datasets in question. Functions
    # included are:
    # 
    #     add_halos() : load Peak Patch halo catalogue
    # 
    #     halos_to_npz() : convert raw binary halo file to python-formatted .npz used by LIM/LAM
    #         Mocker to make LIM mocks.
    #     
    #     add_field() : load Peak Patch initial conditions fields
    #     
    # See bolow for more information on their usage.

    def add_halos(self,catalogue_file=None,mass_cut='highpass',mass_cutoff=None,
            lim=[None, np.array( [[-np.inf,np.inf]]*3,dtype=np.float32 )]):
        # Assign a Peak Patch halo catalogue to the PeakPatch object, load the halos in from file.
        # 
        # Usage:
        # 
        # >>> run = PeakPatch( <run> )
        # >>> 
        # 
        # Optional Parameters:
        # 
        #     catalogue_file : string or None
        #         File name of Peak-Patch formatted dark matter halo catalogue.
        # 
        #     mass_cut : string or None
        #         Option to filter halos by mass. 'High' or 'h' for a highpass filter, 'Low' or 'l'
        #         for low pass filter. If None, no mass cut is done. Default None.
        # 
        #     mass_cutoff : float or None
        #         If mass_cut != None, then mass_cutoff is the cutoff mass. If None, it will try to
        #         determine the cutoff mass M corresponding to 1.5 times the filter radius of the
        #         lattice scale R_latt, where M = 4/3 \pi \rho_{crit} \Omega_m * R_{latt}^3
        # 
        #     lim : [ logical or str , 3x2 NumPy array ]
        #         Option to limit the bounds of the halo catalogue. If lim[0]==False, then all halos
        #         are kept. If lim[0]='Lagrangian' or 'l' then Lagrangian halo positions ioutside of
        #         lim[1][0,0] < xL < lim[1][0,1], lim[1][1,0] < yL <lim[1][1], lim[1][2,0] < zL <
        #         lim[1][2,1] are dumped. Or similarly for Eulerian positions if lim[0]='Eulerian' 
        #         or 'e'.
 
        # Try to load halo catalogue file automatically
        if catalogue_file==None:
            catalogue_file='{0}/output/{1}_merge.pksc.{2}'.format( self.run_dir, self.run_name,
                                                                   self.seed )
            if not os.path.exists(catalogue_file):
                raise OSError( ('file not found: {0}\nUnable to auto-locate halo catalogue file'
                               ).format(catalogue_file) )

        # Try to locate halo catalogue file passed to function
        else:
            if not os.path.exists(catalogue_file):
                raise OSError( ('file not found: {0}').format(catalogue_file) )

        # Open Peak Patch halo catalogue file
        self.catalogue_file = catalogue_file
        c_in = open( catalogue_file , 'rb' )

        # Read halo catalogue header
        self.N_halos  = np.fromfile( c_in, dtype=np.int32  , count=1 )[0]
        self.R_th_max = np.fromfile( c_in, dtype=np.float32, count=1 )[0]
        self.z_obs    = np.fromfile( c_in, dtype=np.float32, count=1 )[0]

        # Raise Error if catalogue is just a header without any halos
        if self.N_halos == 0 or os.path.getsize(catalogue_file) <= 12:
            raise OSError(('The catalogue you have passed: {0}\ncontains no halos. This error is ty'
                          +'pical of a run with a bug during peak\nfinding, meaning no peaks corres'
                          +'ponding to halos could be identified. Check\nlog files for more informa'
                          +'tion:\n{1}/logfiles/.').format(catalogue_file,self.run_dir) )

        # Read number of columns in catalogue
        self.cols  = int( (os.path.getsize(catalogue_file)-12)/4/self.N_halos )

        # Read in Peak Patch halo catalogue (dimension cols x N_halos)
        c = np.reshape( np.fromfile( c_in, dtype=np.float32,
                        count=self.N_halos*self.cols ), (self.N_halos,self.cols) )

        # Throw out halos outside Lagrangian domain
        if lim[0]:
            if   lim[0][0].lower()[0]=='l': q = 0
            elif lim[0][0].lower()[0]=='e': q = 7
            if lim[1].shape==(3,2):
                for j in range(3):
                    if lim[1][j,0] > -np.inf: c=c[ np.where( c[:,j+q] > lim[1][j,0] )[0] , : ]
                    if lim[1][j,1] < +np.inf: c=c[ np.where( c[:,j+q] < lim[1][j,1] )[0] , : ]

        # Save halo catalogue limits for use in plotting functions
        if lim[0]: self.halo_lim = lim[1]
        else     : self.halo_lim = np.array( [[-np.inf,np.inf]]*3,dtype=np.float32 )

        # Constants & cosmological parameters to convert Peak Patch filter scales to masses
        self.H_0      = self.h * (100 * u.km / u.s / u.Mpc)                       # Hubble const
        self.G        = const.G.to( u.Mpc**3 * u.M_sun**-1 * u.s**-2 )            # Gravity const
        self.rho_crit = (3*self.H_0**2/(8*np.pi*self.G)).to( u.M_sun * u.Mpc**-3 )# Critical density
        self.rho_m    = self.rho_crit.value * ( self.Omx + self.OmB )             # Matter density
        # All in units of solar masses and Mpc

        # Perform mass cut
        if mass_cut:

            # Determine whether to do high-pass or low-pass mass filter
            if 'l' in str(mass_cut).lower(): mass_cut='l' # low-pass
            else:                            mass_cut='h' # high-pass

            # Interpolate filter scale cut based on 1.5 x lattice scale
            if not mass_cutoff: mass_cutoff = 1.5 * self.cellsize
            # Convert mass cutoff in solar masses to Peak Patch filter scale cutoff in Mpc
            else: mass_cutoff = ( mass_cutoff/(4/3*np.pi*self.rho_m) )**(1./3.)

            # Perform cut
            if   mass_cut=='l': c = c[ np.where( c[:,6] < mass_cutoff )[0] , : ]
            elif mass_cut=='h': c = c[ np.where( c[:,6] > mass_cutoff )[0] , : ]

        # Allocate halo arrays
        if self.cols >=  3: self.x ,self.y ,self.z  = np.ndarray.view( c[:,  : 3].T )
        if self.cols >=  6: self.dx,self.dy,self.dz = np.ndarray.view( c[:, 3: 6].T )
        if self.cols >=  7: self.R_th               = np.ndarray.view( c[:,  6  ]   )
        if self.cols >=  7: self.M                  = 4./.3*np.pi*self.R_th**3*self.rho_m
        if self.cols >= 10: self.xL,self.yL,self.zL = np.ndarray.view( c[:, 7:10].T )
        if self.cols >= 11: self.f_pk               = np.ndarray.view( c[:,  10 ]   )
        if self.cols >= 13: self.e_v,self.p_v       = np.ndarray.view( c[:,11:13].T )
        if self.cols >= 19:(self.epsxx,self.epsyy,self.epszz,self.epsyz,self.epsxz,self.epsxy
                           ) = np.ndarray.view( c[:,13:19].T )
        if self.cols >= 32:(self.F_d2  , self.zform  , self.gradx  , self.grady  , self.gradz   ,
                            self.gradfx, self.gradfy , self.gradfz , self.Rfclv  , self.FcollvRf,
                            self.F_d2Rf, self.gradrfx, self.gradrfy, self.gradrfz
                           ) = np.ndarray.vew( c[:,19:33].T )

        # Refresh halo count
        self.N_halos = len(c[:,0])



    def add_halo_r_comoving( self, overwrite=False ):
        # This function calculates the comoving distance to each halo. Note that you must run
        # add_halos() before running this to add halo Eulerian positions.
        #    
        # Optional Parameters
        #    
        #     overwrite : bool, default False
        #         Whether to overwrite r_comoving if it has already been calculated.

        # Make sure halos correctly loaded
        if not hasattr(self,'x') or not hasattr(self,'y') or not hasattr(self,'z'):
            raise AttributeError('make sure to call add_halos() method before this method.')

        # Comoving distance to halo
        if not hasattr(self,'r_comoving') or ( hasattr(self,'r_comoving') and overwrite==True ):
            self.r_comoving = np.sqrt( self.x**2 + self.y**2 + self.z**2 )



    def add_N_eff( self, N_eff=None, overwrite=False ):
        # Set the value of the effective number of neutrino species (which is not currently a
        # default Peak Patch parameter). The default is to set it to the Planck 2018 value 
        # N_eff = 3.046.
        #    
        # Optional Parameters
        #    
        #     N_eff : float or None, default None
        #    
        #     overwrite : bool, default False
        #         Whether to overwrite an existing value.
     
        # if no value passed, use the default value from Planck 2018 results, N_eff=3.046
        if N_eff == None and ( not hasattr(self,'N_eff') 
                               or (hasattr(self,'N_eff') and overwrite==True) ):
            self.N_eff = 3.046

        # otherwise, assign the value passed to the function
        elif not hasattr(self,'N_eff') or ( hasattr(self,'N_eff') and overwrite==True ):
            self.N_eff = N_eff



    def add_mnu( self, mnu=None, overwrite=False ):
        # Set the value of the total mass of massive neutrino species (which is not always a
        # parameter included in the Peak Patch parameter file). If no value is passed, the default
        # value of 0.0 is used.
        # 
        # Optional Parameters
        # 
        #     mnu : float or None, default None
        #         The value to set self.mnu to. The default is None, in which case it will be set to
        #         0.0.
        # 
        #     overwrite : bool, default False
        #         Whether to overwrite an existing value.

        # If no values passed, use default value of 0.0
        if mnu == None and ( not hasattr(self,'mnu') or (hasattr(self,'mnu') and overwrite==True) ):
            self.mnu = 0.0

        # otherwise, assign the value passed to the function
        elif not hasattr(self,'mnu') or ( hasattr(self,'mnu') and overwrite==True ):
            self.mnu = mnu



    def add_halo_redshifts( self, z_max=3402.0, mnu=None, N_eff=None, overwrite=False ):
        # This function calculates the redshift of each halo in a in Peak Patch halo catalogue.
        # 
        # Optional Parameters
        # 
        #     z_max : float, default 3402.0
        #         Astropy requires a maximum redshift for its calculation of the redshift as a
        #         function of comoving distance. The default value is the redshift of matter
        #         radiation equality. Typically if you go to higher redshift than this, the
        #         equations of motion become more compicated as you need to consider regimes with
        #         multiple dominant types of energy.
        # 
        #     mnu : float or None, default None
        #         The total mass of massive neutrino species. Default is None in which case 0.0 will
        #         be used.
        # 
        #     N_eff : float or None, default None
        #         The effective number of neutrino species (which is not currently a default Peak
        #         Patch parameter. If set to None, which is the default, the value from the Planck
        #         2018 results N_eff = 3.046 is used.
        # 
        #     overwrite : bool, default False
        #         Whether to overwrite z_halo if it has already been calculated.

        # Make sure halos correctly loaded
        if not hasattr(self,'x') or not hasattr(self,'y') or not hasattr(self,'z'):
            raise AttributeError('make sure to call add_halos() method before this method.')

        # Add r_comoving if it hasn't already been added
        if not hasattr(self,'r_comoving'):
            self.add_halo_r_comoving(overwrite=overwrite)

        # Add the total mass of massive neutrino species
        self.add_mnu(mnu, overwrite=False)

        # Add effective number of massive neutrino species N_eff
        self.add_N_eff(N_eff, overwrite=False)

        # Halo redshift (note that we just use the effective number of Neutrino species from Planck
        # 2018 as this is not a stock Peak Patch parameter)
        tab = z_of_r_comoving_table( zmin=1.0e-3, zmax=z_max, dlog10z=6.538e-3, Omega_r0=0.0,
                Omega_m0=self.Omx+self.OmB, Omega_k0=0.0, Omega_Lambda=self.Omvac, H_0=self.h*100,
                c=299792458.0 )
        self.z_halo = z_of_r_comoving( self.r_comoving, z_chi_table=tab )



    def add_merged_catalogue_file(self):
        # This method finds the merged catalogue file for this Peak Patch run. This is useful if you
        # want to make formatted .npz files for instance to run LIM/LAM Mocker and get LIM Mocks.

        # Check that the output directory exists
        if not os.path.exists('{0}/output'.format(self.run_dir)):
            raise FileNotFoundError('This run has no output folder.')

        # Guess catalogue file name based on typical Peak Patch formatting
        self.catalogue_file='{0}/output/{1}_merge.pksc.{2}'.format( self.run_dir, self.run_name,
                self.seed )

        # If that format can't be found, see if you can find any file with merge in its name
        if not os.path.exists(self.catalogue_file):
            for f in os.listdir(self.run_dir+'/output'):
                if 'merge' in f:
                    self.catalogue_file = '{0}/output/{1}'.format( self.run_dir, f )
                    break
        if not os.path.exists(self.catalogue_file):
            raise OSError( ('file not found: {0}\nUnable to auto-locate halo catalogue file'
                           ).format(self.catalogue_file) )



    def add_raw_catalogue_files(self):
        # This method finds all raw (pre-merge) Peak Patch catalogue files in a run directory.

        # Check that the output directory exists
        if not os.path.exists('{0}/output'.format(self.run_dir)):
            raise FileNotFoundError('This run has no output folder.')

        # Guess raw catalogue file name based on typical Peak Patch formatting and find all files
        raw='{0}_raw.pksc.{1}'.format( self.run_name, self.seed )
        if os.path.exists( '{0}/output/{1}'.format( self.run_dir, raw ) ):
            self.raw_catalogue_files = []
            for f in os.listdir(self.run_dir+'/output'):
                if raw in f:
                    self.raw_catalogue_files += [ '{0}/output/{1}'.format( self.run_dir, f ) ]

        # If no raw catalogue files found with standard name, just add anything with raw in the
        # filename
        else:
            if any( 'raw' in f for f in os.listdir(self.run_dir+'/output') ):
                self.raw_catalogue_files = []
                for f in os.listdir(self.run_dir+'/output'):
                    if 'raw' in f:
                        self.raw_catalogue_files += [ '{0}/output/{1}'.format( self.run_dir, f ) ]
            else:
                raise FileNotFoundError('could not find raw catalogues.')



    def halos_to_field( self, save_out=True, overwrite=False, mass_bins=None, binning='log',
            filter_kernel='scipy', N=None, lims=None, **kwargs ):
        # mass_bins can be bin edges or an integer or None.

        # filter_kernel = 'top-hat', 'gaussian', 'einasto', 'nfw'

        field_file = '{0}/tables/halo_field.npz'.format(self.run_dir)
        if overwrite == False:
            if os.path.exists( field_file ):
                npzfile = np.load( field_file, allow_pickle=True )
                f = [ npzfile['field'], npzfile['axes'] ]
                self.halo_field = f
                return f[0], f[1]

        # Set limits on the halos to be converted to field
        if lims is None:
            lims = np.array([[ self.cenx-self.boxsize/2, self.cenx+self.boxsize/2 ],
                             [ self.ceny-self.boxsize/2, self.ceny+self.boxsize/2 ],
                             [ self.cenz-self.boxsize/2, self.cenz+self.boxsize/2 ] ])

        # Apply halo limits
        where = np.argwhere( (self.x>=lims[0,0]) & (self.x<=lims[0,1]) )[:,0]
        M,R_th = self.M[where], self.R_th[where]
        x,y,z  = self.x[where], self.y[where], self.z[where]
        where = np.argwhere( (y>=lims[1,0]) & (y<=lims[1,1]) )[:,0]
        M,R_th,x,y,z = M[where], R_th[where], x[where], y[where], z[where]
        where = np.argwhere( (z>=lims[2,0]) & (z<=lims[2,1]) )[:,0]
        M,R_th,x,y,z = M[where], R_th[where], x[where], y[where], z[where]
        
        # Set the number of mass smoothing bins
        if mass_bins is None:
            mass_bins = 20

        # Make bins yourself
        if isinstance(mass_bins,int):
            if binning=='log':
                mass_bins = np.logspace( np.log10(np.min(self.M)) ,
                                         np.log10(np.max(self.M)) , mass_bins+1 )
            else:
                mass_bins = np.linspace( np.min(self.M) , np.max(self.M) , mass_bins+1 )

        # 
        if N is None:
            if lims is None:
                N = (self.neff,self.neff,self.neff)
            else:
                N = [ int((lims[j,1]-lims[j,0])/self.cellsize) for j in range(3) ]
        elif isinstance(N,int):
            N = (N,N,N)

        # If this flag is set to True, we use scipy's Gaussian field routine rather than the python
        # FFTW wrapper
        if filter_kernel == 'scipy':
            from scipy.ndimage import gaussian_filter
        else:
            import pyfftw

        # 
        rho_m = np.zeros( (N[0],N[1],N[2]) )
        for j in range(len(mass_bins)-1):
            if j==0:
                where = np.argwhere( (M>=mass_bins[0]) & (M<=mass_bins[  1]) )[:,0]
            else:
                where = np.argwhere( (M> mass_bins[j]) & (M<=mass_bins[j+1]) )[:,0]

            f = np.histogramdd(
                    sample=np.array([ x[where] , y[where] , z[where] ]).T,
                    bins=( np.linspace( lims[0,0], lims[0,1], N[0]+1 ),
                           np.linspace( lims[1,0], lims[1,1], N[1]+1 ),
                           np.linspace( lims[2,0], lims[2,1], N[2]+1 ) ),
                    weights=M[where]
                    )[0]

            if filter_kernel == 'scipy':
                
                # First grab mean Lagrangian space top-hat filter radius of all halos in the mass bin
                R_th_j = np.mean( R_th[where] )

                # The Gaussian and top-hat filters follow the relation W_G(k;R_G) = W_TH(k;2 R_th)
                R_g_j = 2 * R_th_j

                # Get the Eulerian space radius by assuming the density of the resultant peak is 200 times the initial
                R_g = 200**(-1/3) * R_g_j

                # Filter the matter field
                f = gaussian_filter( f, sigma = np.array([ R_g * N[j] / (lims[j,1]-lims[j,0])
                                                           for j in range(3) ]),
                                     order=0, truncate=5.0 )

                rho_m += f

            else:

                # Initialise FFTW
                nyquist         = int(N[3]/2)+1
                f_rfftwf        = pyfftw.empty_aligned( (N[0],N[1],N[2]     ), dtype='float32'   )
                fk_rfftwf       = pyfftw.empty_aligned( (N[0],N[1],N[3]//2+1), dtype='complex64' )
                f_irfftwf       = pyfftw.empty_aligned( (N[0],N[1],N[2]     ), dtype='float32'   )
                rfftwf_forward  = pyfftw.FFTW( f_rfftwf, fk_rfftwf, axes=(0,1,2),
                                           direction='FFTW_FORWARD', threads=4 )
                rfftwf_backward = pyfftw.FFTW( fk_rfftwf, f_irfftwf, axes=(0,1,2),
                                           direction='FFTW_BACKWARD', threads=4 )

                # Fill initial array with values
                f_rfftwf[:] = f
                #for ii in range(len(f_rfftwf)):
                #    for jj in range(len(f_rfftwf)):
                #        for kk in range(len(f_rfftwf)):
                #            f_rfftwf[ii,jj,kk] = f[ii,jj,kk]
    
                # Perform FFT to get Fourier transform of f
                f_k = rfftw_forward()

                kbins = np.arange(dk,dk*(klen+2),dk) #np.logspace( np.log10(dk), np.log10(kmax), neff )
                P_k   = np.zeros(len(kbins))
                P_0   = np.zeros(len(P_k  ))

                for ii in range(N[0]):
                    if ii<nyquist: iii = ii
                    else:          iii = nyquist-ii
                    for jj in range(N[1]):
                        if jj<nyquist: jjj = jj
                        else:          jjj = nyquist-jj
                        for kk in range(N[2]):
                            if kk<nyquist: kkk = kk
                            else:          kkk = nyquist-kk

                            p = np.sqrt( iii**2 + jjj**2 + kkk**2 )
                            l = int(np.floor(p))
                            c = (1-np.array([ l-p , l-p+1 ])**2)**2

                            P_k[l:l+2] = P_k[l:l+2] + c * np.abs(f_k[i,j,kk])**2
                            P_0[l:l+2] = P_0[l:l+2] + c
                            # That measures P(k), which we don't actually need
   
                            # What we actually need is to multiply f_k by a window function
                            # and then ifft back, then we add together the different things
                            # need to find a way around these for loops

                # Transform back to get smoothed field
                f_W = rfftw_backward()

        # Bin edges for discretisation of the halo field
        axes = [ np.linspace( lims[0,0], lims[0,1], N[0]+1 ),
                 np.linspace( lims[1,0], lims[1,1], N[1]+1 ),
                 np.linspace( lims[2,0], lims[2,1], N[2]+1 )  ]

        # Convert from mass to density
        voxel_volume = ( ( axes[0][1] - axes[0][0] ) *
                         ( axes[1][1] - axes[1][0] ) *
                         ( axes[2][1] - axes[2][0] )   )
        rho_m /= voxel_volume

        # Save the halo field as an attribute of the PeakPatch class object self
        self.halo_field = [ rho_m , axes ]

        # Save the halo field as a file
        if save_out == True:
            np.savez( field_file, field=rho_m, axes=axes, allow_pickle=True )

        # Return the halo field
        return rho_m, axes



    def plot_halo_field_slab_diff( minuend, subtrahend, fig, ax, plot_type='pcolormesh', plane='xz',
            intercept=None, lims=None, **kwargs ):

        # 
        if not ( minuend.halo_field[1] == subtrahend.halo_field[1] ).all():
            raise ValueError('The two fields must have halo fields with the same axes.')

        # 
        field  = []
        field += [ minuend.halo_field[0] - subtrahend.halo_field[0] ]
        field += [ minuend.halo_field[1] ]
        return plot_halo_field_slab( minuend, fig, ax, field=field, plot_type=plot_type,
                plane=plane, intercept=intercept, lims=lims, **kwargs )



    def plot_halo_field_slab( self, fig, ax, field=None, plot_type='pcolormesh', plane='xz',
            intercept=None, lims=None, **kwargs ):
        # Function for plotting Peak Patch initial conditions fields from attributes rho, rhog,
        # zeta, zetag, and chi or from a field passed at command line. Note that these fields should
        # be formatted as PeakPatch field slice attributes are (see add_field above).
        # 
        # Parameters
        # 
        #     fig : matplotlib figure
        #         The figure in which to plot the field slice.
        # 
        #     ax : matplotlib.axes._subplots.Axessubplot
        #         An axes subplot object (e.g. one made using `fig,ax=plt.subplots`) in which to
        #         plot the 2D histogram, one of the axes contained within the figure, fig.
        #
        #     field : None or [ np.ndarray, np.ndarray ], default None
        #         If None, slices are plotted from the field attributes of the PeakPatch object (see
        #         parameter field_type above for more on this). Otherwise it must be a list format-
        #         ted like a PeakPatch field attributes.
        # 
        #     plot_type : str, default 'pcolormesh'
        #         The type of plot to make: 'pcolormesh' for plt.pcolormesh, 'contour' for
        #         plt.contour.
        # 
        #     plane : str, default 'xz'
        #         If self has a full 3D field loaded instead of a single slice, then a plane of that
        #         full field is chosen. This parameter is only used if field[1]==None (meaning that
        #         a 3D field has been loaded). The default value is 'xz' to plot the x-z plane. The
        #         other allowed values are 'xy' and 'yz'.
        # 
        #     **kwargs magic variables passed as variables and treated as a dict
        #         Allows for setting of key word arguemnets in the matplotlib pcolormesh method,
        #         such as vmin and vmax so taht you can keep colour bars consistent between plots.
        #         Supported kwargs for this function:
        #         
        #         - xlabel : str
        #               The label for the horizontal axis. If not set, use the default:
        #               "r'comoving distance ${0}~[\mathrm{{Mpc}}]$'"
        #               formatted according to the slice passed. This option is useful if you want
        #               to plot several panes with the same axis type to suppress labels and tidy up
        #               your plots.
        #         
        #         - ylabel : string
        #               The label for the vertical axis. If not set, use the default:
        #               "r'comoving distance ${0}~[\mathrm{{Mpc}}]$'"
        #               formatted according to the slice passed. This option is useful if you want
        #               to plot several panes with the same axis type to suppress labels and tidy up
        #               your plots.
        #         
        #         - xlim : tuple of floats, length 2
        #               Axis limits for the horizontal axis.
        #         
        #         - ylim : tuple of floats, length 2
        #               Axis limits for the vertical axis
        #         
        #         - vmin : float
        #               Lower limit for the colour bar.
        #         
        #         - vmax : float
        #               Upper limit for the colour bar.
        #         
        #         - cbarlabel : str
        #               The label for the colour bar. If none passed, a label is chosen based on the
        #               field passed. Like with the axis label kwargs, this option is usefull if you
        #               want to suppress automatic labels.
        # 
        #         - cbarformat : str
        #               Set to None to use automatic colour bar formatting, 'scinot' for scientific
        #               notation, or 'log' for log scale.

        # Check that the halo field has been generated
        if not hasattr(self,'halo_field'):
            self.halos_to_field( save_out=True, overwrite=False, mass_bins=None, binning='log',
                    filter_kernel='scipy', N=None, lims=None, **kwargs )

        # Select field if none given as argument
        if field is None:
            field  = []
            field += [ self.halo_field[0] ]
            field += [ self.halo_field[1] ]

        # 
        if lims is None:
            lims = np.array([ [ field[1][0][0], field[1][0][-1] ],
                              [ field[1][1][0], field[1][1][-1] ],
                              [ field[1][2][0], field[1][2][-1] ] ])

        # Check that limits and intercepts are in bounds
        if ( lims[0,0] < field[1][0][0] or lims[0,1] > field[1][0][-1] or 
                lims[1,0] < field[1][1][0] or lims[1,1] > field[1][1][-1] or
                lims[2,0] < field[1][2][0] or lims[2,1] > field[1][2][-1] ):
            raise ValueError('Limits out of bounds.')

        # Set the indices of the axis limits
        ilo, ihi = 0, len(field[1][0])
        jlo, jhi = 0, len(field[1][1])
        klo, khi = 0, len(field[1][2])
        for i in range(0, ihi):
            ilo = i
            if field[1][0][ilo] >= lims[0,0]:
                break
        for i in range(ihi, ilo, -1):
            ihi = i
            if field[1][0][ihi-1] <= lims[0,1]:
                break
        for j in range(0, jhi):
            jlo = j
            if field[1][1][jlo] >= lims[1,0]:
                break
        for j in range(jhi, jlo, -1):
            jhi = j
            if field[1][1][jhi-1] <= lims[1,1]:
                break
        for k in range(0, khi):
            klo = k
            if field[1][2][klo] >= lims[2,0]:
                break
        for k in range(khi, klo, -1):
            khi = k
            if field[1][2][khi-1] <= lims[2,1]:
                break

        # Consider only the portion of the field within the limits
        f  = []
        f += [ field[0][ ilo:ihi, jlo:jhi, klo:khi ] ]
        f += [ [ field[1][0][ilo:ihi],
                 field[1][1][jlo:jhi],
                 field[1][2][klo:khi]  ] ]
 
        # Used later to set default axis labels
        l = r'comoving distance ${0}~[\mathrm{{Mpc}}]$'

        # Setup for a yz plane
        if plane == 'yz':

            # Horizontal and vertical axis labels and bin edges for the plot
            xlabel , ylabel = l.format('y') , l.format('z')
            axes = ( f[1][1] , f[1][2] )

            # the projectino of the map
            voxel_width = r[1][0][ 1] - f[1][0][0]
            slab_depth  = f[1][0][-1] - f[1][0][0]
            proj = np.sum( f[0], axis=0 ).T * voxel_width / slab_depth
            # By convention, all field slices are viewed as if from "above", in other words,
            # the normal vector to the plane being plotted (in this case the yz plane)
            # oriented along the Cartesian basis vector (in this case x) points to the
            # observer. To achieve this, the transpose is taken above.

        # Setup for an xz plane
        elif plane == 'xz':

            # Horizontal and vertical axis labels and bin edges for the plot
            xlabel , ylabel = l.format('z') , l.format('x')
            axes = ( f[1][2] , f[1][0] )

            # set field[0] to just the relevant slice of the field
            voxel_width = f[1][1][ 1] - f[1][1][0]
            slab_depth  = f[1][1][-1] - f[1][1][0]
            proj = np.sum( f[0], axis=1 ) * voxel_width / slab_depth
            # By convention, all field slices are viewed as if from "above", in other words,
            # the normal vector to the plane being plotted (in this case the xz plane)
            # oriented along the Cartesian basis vector (in this case y) points to the
            # observer. To achieve this, the transpose is NOT taken above.  

        # Setop for an xy plane
        else:

            # Horizontal and vertical axis labels and bin edges for the plot
            xlabel , ylabel = l.format('x') , l.format('y')
            axes = ( f[1][0] , f[1][1] )

            # set field[0] to just the relevant slice of the field
            voxel_width = f[1][2][ 1] - f[1][2][0]
            slab_depth  = f[1][2][-1] - f[1][2][0]
            proj = np.sum( f[0], axis=2 ).T * voxel_width / slab_depth
            # By convention, all field slices are viewed as if from "above", in other words,
            # the normal vector to the plane being plotted (in this case the xy plane)
            # oriented along the Cartesian basis vector (in this case z) points to the
            # observer. To achieve this, the transpose is taken above.  
                    

        # Overwrite default axis labels if something else is passed in **kwargs
        if 'xlabel' in kwargs: xlabel = kwargs['xlabel']
        if 'ylabel' in kwargs: ylabel = kwargs['ylabel']

        # Set the colour map
        if 'cmap' in kwargs:

            # Use Planck colour map
            if kwargs['cmap'] == 'planck': cmap = planck_cmap()

            # Use Earth Tones colour map
            elif kwargs['cmap'] == 'earth_tones': cmap = earth_tones_cmap()

            # Use matplotlib default colour maps
            else: cmap = kwargs['cmap']

        # Use the default colour map (viridis)
        else:
            cmap = 'inferno'

        # Make colour mesh object
        if plot_type == 'pcolormesh':

            ax.set_aspect(1)
            cax  = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)

            # Put colour bar axis in scientific notation with power mulitplying the scale above
            if 'cbarformat' in kwargs: cbarformat = kwargs['cbarformat']
            else                     : cbarformat = 'log'

            # Make the colour map
            if 'vmin' in kwargs and 'vmax' in kwargs:
                vmin, vmax = kwargs['vmin'], kwargs['vmax']
            else:
                vmin, vmax = np.min(proj), np.max(proj)

            if vmin > np.min(proj) and vmax < np.max(proj):
                extend = 'both'
            elif vmin > np.min(proj) and vmax >= np.max(proj):
                extend = 'min'
            elif vmin <= np.min(proj) and vmax < np.max(proj):
                extend = 'max'
            else:
                extend = 'neither'

            # 
            if cbarformat == 'log':
                if vmin <= 0:
                    vmin = np.min( proj[ proj > 0 ] )/10
                from matplotlib.colors import LogNorm
                proj[ proj <= 0 ] = np.min( proj[ proj > 0 ] )/10
                pcm = ax.pcolormesh( axes[0] , axes[1] , proj , cmap=cmap,
                                     norm=LogNorm(vmin=vmin, vmax=vmax) )
            else:
                pcm = ax.pcolormesh( axes[0] , axes[1] , proj , cmap=cmap,
                                     vmin=vmin , vmax=vmax )
            cbar = fig.colorbar( pcm, cax=cax, extend=extend )

            if cbarformat == 'scinot':
                from matplotlib.ticker import ScalarFormatter
                formatter = ScalarFormatter(useMathText=True)
                formatter.set_powerlimits((0,0))
                cbar.ax.yaxis.set_major_formatter(formatter)

            # Colour bar labels
            if 'cbarlabel' in kwargs:
                cbarlabel = kwargs['cbarlabel']
            else:
                cbarlabel = r'$\rho_\mathrm{halo} ~ [M_\odot \, \mathrm{Mpc}^{-3}]$'
            if 'rotation' in kwargs: rotation = kwargs['rotation']
            else:                    rotation = 90
            cbar.set_label( cbarlabel, rotation=rotation, labelpad=5 )

            # Set axis labels
            if xlabel != '': ax.set_xlabel(xlabel)
            if ylabel != '': ax.set_ylabel(ylabel)

            # Define the type of axes object used by make_fig() to make figure files
            axtype = 'Peak Patch halo field slab'


            # Returns pcolormesh and labels
            return { 'fig'           : fig       ,
                     'ax'            : ax        ,
                     'axtype'        : axtype    ,
                     'pcolormesh'    : pcm       ,
                     'cax'           : cax       ,
                     'cbar'          : cbar      ,
                     'xlabel'        : xlabel    ,
                     'ylabel'        : ylabel    ,
                     'colorbarlabel' : cbarlabel }
        
        # Make contours
        elif plot_type == 'contour':
            axtype = 'Peak Patch halo field slab contour'
            topo = ax.contour( axes[0] , axes[1] , proj ) # , levels = ...

            return { 'ax'     : ax     ,
                     'axtype' : axtype ,
                     'contour': topo   }



    def halos_to_npz(self, fileout=None, overwrite=False, zform_max=3402, mnu=None):
        # Convienence function for converting Peak Patch's native unformatted binary halo catalogue
        # files to .npz python readable files that are used in LIM/LAM Mocker to produce WebSky line
        # intensity mock maps. This function must be called after a halo catalogue has been read in
        # using the method add_halos().
        # 
        # Usage:
        # 
        # >>> from peakpatchtools import PeakPatch
        # >>> run = PeakPatch( <run> )
        # >>> run.add_halos()
        # >>> run.halos_to_npz()
        # 
        # Optional Parameters:
        # 
        #     fileout : string or None (default None)
        #         File name of .npz formatted dark matter halo catalogue to be output. The default
        #         is None in which case the file name is set to the same as the unformatted binary
        #         file followed by the file extension .npz.
        # 
        #     overwrite : bool (default False)
        #         If the file specified by fileout exists and overwrite is set to True, an existing
        #         .npz file is overwritten. If set to False, no additional .npz is produced and a
        #         message is printed to screen indicating as such.
        # 
        #     zform_max : float (default 3402)
        #         This floating point number serves two purposes. The redshift interpollator
        #         requires a maximum redshift, which is set here. Additionally, the .npz files read
        #         in by the LIM/LAM Mocker code expect a parameter zform indicating the redshift
        #         that each DM halo forms at, however, the formation redshift is not always included
        #         in Peak Patch halo catalogues, so if no formation redshift is included in the
        #         existing halo catalogue, this value is used for all halos. The default value is
        #         the redshift at which matter domination begins (effectively the earliest a halo
        #         could form under most conditions).
        # 
        #     mnu : float (default None)
        #         If no neutrino mass is specified in the parameter file, it can be 

        # Automatically set the output file name
        if fileout == None:
            fileout = self.catalogue_file + '.npz'
        
        # Check that halo catalogue has been properly read in
        attributes = [ 'OmB', 'Omx', 'Omvac', 'h', 'ns', 'sigma8', 'x', 'y', 'z', 'N_halos', 'M',
                'xL', 'yL', 'zL', 'dx', 'dy', 'dz' ]
        for attribute in attributes:
            if not hasattr(self,attribute):
                raise AttributeError(justify_text(('catalogue missing attribute '+attribute+'. Mak'+
                                        'e sure to call method add_halos() before halos_to_npz().'),
                                width=100, first_line_buffer='AttributeError' ))

        # Neutrino mass is not actually necessary for the Peak Patch calculation and is sometimes
        # omitted from the parameter file, for this reason we set it to zero here
        if not hasattr(self,'mnu'): self.mnu=0.0

        # Check that fileout is formatted correctly
        if not isinstance(fileout,str):
            raise TypeError('fileout has to be a string.')
        elif fileout[-4:] != '.npz':
            raise ValueError('fileout must end in ".npz".')

        # Check if fileout exists and overwrite if overwrite==True
        if os.path.exists(fileout) and overwrite==False:
            print('.npz file detected and overwrite set to False, exitting.')
            return
        elif os.path.exists(fileout) and overwrite==True:
            os.system('rm '+fileout)

        # Header for .npz halo file
        cosmo_header = { 'Omega_M' : self.OmB + self.Omx ,
                         'Omega_B' : self.OmB            ,
                         'Omega_L' : self.Omvac          ,
                         'h'       : self.h              ,
                         'ns'      : self.ns             ,
                         'sigma8'  : self.sigma8          }

        # Get redshift for each halo
        self.add_halo_redshifts( z_max=3402.0 )

        # If formation redshifts are included in the catalogue, include these, otherwise, set
        # formation redshift to arbitrarily large value (default is redshift of beginning of matter
        # domination)
        if hasattr(self,'zform'):
            zform = self.zform
        else:
            zform = [ zform_max for i in range(self.N_halos) ]

        # Save the halos as .npz
        np.savez( fileout, Nhalo=self.N_halos, M=self.M, x=self.x, y=self.y, z=self.z,
                xL=self.xL, yL=self.yL, zL=self.zL, vx=self.dx, vy=self.dy, vz=self.dz,
                cosmo_header=cosmo_header, zform=zform, zhalo=self.z_halo )

        # Store the output file as an object in the Peak Patch class
        self.catalogue_file_npz = fileout



    def pksc_to_npz( self, pksc_file=None, npz_file=None, Mmin=None, z_range=None, parallel_cores=40
            ):
        # Script for converting the unformatted binary dark matter halo catalogues used by Peak
        # Patch to NumPy formatted .npz files that are used in LIMLAM Mocker to make mock
        # observations for line intensity mapping experiments.
        # 
        # Arguments:
        # 
        #     pksc_file : string or None, default None
        #         The unformatted Peak Patch dark matter halo catalogue file. If None, the catalogue
        #         file will be located automatically for the PeakPatch object.
        # 
        #     npz_file : string or None, default None
        #         The file at which to save the .npz halo catalogue. If None, the file will just be
        #         the pksc_file followed by ".pnz".
        # 
        #     Mmin : float or None, default None
        #         A minimum mass cut to perform on the unformatted catalogue such that halos below a
        #         given mass are not saved to the new catalogue. If a float is passed, this is
        #         interpretted as the cutoff mass in solar masses. If None is passed, the default
        #         mass cut is applied (any halo whose Lagrangian radius is less than 1.5 times the
        #         fundamental resolution of the box).
        # 
        #     z_range : NumPy array or None, default None
        #         The .npz catalogues also include halo redshifts. Only halos within this redshift
        #         range are included. If None, then all halos are added to the .npz file regardless
        #         of their redshift.
        # 
        #     parallel_cores : integer, default 40
        #         The number of parallel cores to use in the the pksc to npz calculation. The
        #         default is 40, the number of cores on each Niagara node.
        # 
        # Usage:
        #     
        #     This function is meant mostly to be run in parallel on Niagara using the script
        #     pksc2npz.sh in the scripts subdirectory of the Peak Patch repository:
        #
        #     $ cd .../peakpatch
        #     $ ./output/pksc2npz /path/to/run/dir/

        #from multiprocessing import Pool
        #from astropy import config
        #config.set_temp_config( sys.argv[1] )

        # Generate interpolation table for redshift to comoving distance calculation
        z_chi_tab = z_of_r_comoving_table( zmin=1.0e-3, zmax=3402.0, dlog10z=6.538e-3, Omega_r0=0.0,
                Omega_m0=self.Omx+self.OmB, Omega_k0=0.0, Omega_Lambda=self.Omvac, H_0=self.h*100,
                c=299792458.0 )

        # Set the input halo file
        if pksc_file == None:
            if hasattr(self,'mergepkoutfile'):
                pksc_file = self.mergepkoutfile
            else:
                import glob
                pksc_file = glob.glob( self.run_dir + '/output/*merge*' )
                if len(pksc_file) > 0:
                    pksc_file = pksc_file[0]
                else:
                    raise AttributeError('could not auto-locate halo catalogue file.')

        # If no .npz file is passed, automatically set it
        if npz_file == None:
            npz_file = pksc_file+'.npz'

        # Header for .npz file
        cosmo_header = { 'Omega_M' : self.Omx+self.OmB , 'Omega_B' : self.OmB ,
                'Omega_L' : self.Omvac , 'h' : self.h , 'ns' : self.ns , 'sigma8' : self.sigma8 }

        # Cosmological parameters
        H_0   = self.h * 100 # km/s/Mpc
        G     = 4.300917270036279e-09 # km^2 Mpc Msol^-1 s^-2
        rho_c = 3*H_0**2 / (8*np.pi*G) # Msol Mpc^-3
        rho_m = rho_c * ( self.Omx + self.OmB ) # Msol Mpc^-3

        # Interpolate filter scale cut based on 1.5 x lattice scale
        if not Mmin:
            mass_cutoff = 1.5 * self.cellsize
            
        # Convert mass cutoff in solar masses to Peak Patch filter scale cutoff in Mpc
        else:
            mass_cutoff = ( Mmin / (4/3*np.pi*rho_m) )**(1./3.)

        # Parallelize the funciton poolfunc_pksc_to_npz() defined below across `parallel_cores` CPUs
        # on the cluster
        # Pool( parallel_cores ).starmap(poolfunc_pksc_to_npz, [( pksc_file , npz_file , mass_cutoff ,
        #                         cosmo_header, rho_m , z_range, z_chi_tab )] )
        poolfunc_pksc_to_npz( pksc_file , npz_file , mass_cutoff,  cosmo_header, rho_m , z_range, 
                z_chi_tab )



    def pksc_to_npz_niagara( self, mass_cutoff=None, halo_zmin=-1, halo_zmax=-1,
            table_zmin=1.0e-3, table_zmax=3402.0, table_dlog10z=6.538e-3, tlimit='1:00:00',
            submit=True ):
        # Script for converting the unformatted binary dark matter halo catalogues used by Peak
        # Patch to NumPy formatted .npz files that are used in LIMLAM Mocker to make mock
        # observations for line intensity mapping experiments.
        # 
        # Arguments:
        # 
        #     mass_cutoff : float or None, default None
        #         A minimum mass cut to perform on the unformatted catalogue such that halos below a
        #         given mass are not saved to the new catalogue. If a float is passed, this is
        #         interpretted as the cutoff mass in solar masses. If None is passed, the default
        #         mass cut is applied (any halo whose Lagrangian radius is less than 1.5 times the
        #         fundamental resolution of the box).
        # 
        #     halo_zmin : float, default -1
        #     halo_zmax : float, default -1
        #         The .npz catalogues also include halo redshifts. Only halos within this redshift
        #         range are included. If either bound is set to -1, then no bound is used.
        # 
        #     table_zmin : float, default 1.0e-3
        #     table_zmax : float, default 3402.0
        #     table_dlog10z : float, default 6.538e-3
        #         Bounds of the interpollation table used to determine redshift as a function of
        #         comoving distance.
        # 
        #     tlimit : str, default '1:00:00'
        #         The time limit for the Niagara job. The default corresponds to a 1 hour job. If
        #         the time limit is set to less than or equal to 1 hour then the job will be
        #         sumbitted to the debug queue.
        # 
        #     submit : bool, default True
        #         Whether to submit the jobscript to the Niagara scheduler.
        # 
        # Usage:
        #     
        #     This function is meant mostly to be run on Niagara using the script. Change
        #     directories to the Peak Patch repository home directory and then run the following:
        #     
        #     $ ml python
        #     $ python3
        #     >>> import peakpatchtools as pkp
        #     >>> r = PeakPatch( '/path/to/peak/patch/run/' )
        #     >>> r.pksc_to_npz_niagara()
        # 
        #     This makes a local copy of the template submit script and runs it on Niagara. The
        #     template script runs the function load_halos_pksc_to_npz() and executes the python
        #     script peakpatch/python/scripts/pksc2npz.py.

        # Set the Peak Patch DM halo catalogue file
        if not hasattr(self,'catalogue_file'):
            self.add_merged_catalogue_file()

        # Set the name of the .npz output file
        npz_file = self.catalogue_file+'.npz'

        # Interpolate filter scale cut based on 1.5 x lattice scale
        if mass_cutoff is None:
            mass_cutoff = 1.5*self.cellsize

        else:
            # Cosmological parameters
            H_0   = self.h * 100 # km/s/Mpc
            G     = 4.300917270036279e-09 # km^2 Mpc Msol^-1 s^-2
            rho_c = 3*H_0**2 / (8*np.pi*G) # Msol Mpc^-3
            rho_m = rho_c * ( self.Omx + self.OmB ) # Msol Mpc^-3

            # Convert mass cutoff in solar masses to Peak Patch filter scale cutoff in Mpc
            mass_cutoff = ( Mmin / (4/3*np.pi*rho_m) )**(1./3.)

        # Get the run time limit in hours, minutes and seconds
        seconds = [ float(j) for j in tlimit.split(':')[::-1] ]
        seconds = int( np.sum( np.array([ seconds[j]*60**j for j in range(len(seconds)) ]) ) )
        hours, minutes, seconds = int((seconds//60.)//60.) , int((seconds//60.)%60) , int(seconds%60.)
        tlimit = '{0}:{1:02}:{2:02}'.format(hours,minutes,seconds)

        # Debug flag
        if hours <= 1:
            debug_string = '#SBATCH -p debug'
        else:
            debug_string = ''

        # Generate interpolation table for redshift to comoving distance calculation
        z_chi_tab = z_of_r_comoving_table( zmin=table_zmin, zmax=table_zmax, dlog10z=table_dlog10z,
                Omega_r0=0.0, Omega_m0=self.Omx+self.OmB, Omega_k0=0.0, Omega_Lambda=self.Omvac,
                H_0=self.h*100, c=299792458.0 )

        # Save the interpollation table of redshift and comoving distance values to a file
        z_chi_tab_file = self.run_dir+'/output/z_chi_tab.npz'
        np.savez( z_chi_tab_file, z_chi_tab=z_chi_tab )    

        # Niagara job I/O
        jobname                = 'pkp2npz'
        bashoutput             = './pksc2npz.out'
        submit_dir             = self.run_dir+'/output/'
        submit_script_template = peak_patch_dir+'/templates/pksc2npz_niagara_submit_script.sh'
        submit_script          = submit_dir+'/pksc2npz_niagara_submit_script.sh'
 
        # Make Niagara jobscript
        os.system( 'cd {0};cp {1} {2}'.format(submit_dir,submit_script_template,submit_script) )

        # Set Niagara I/O values in jobscript
        os.system( 'sed -i "s!DEBUG!{0}!g" {1}'      .format(debug_string  ,submit_script) )
        os.system( 'sed -i "s!TIME_LIMIT!{0}!g" {1}' .format(tlimit        ,submit_script) )
        os.system( 'sed -i "s!JOB_NAME!{0}!g" {1}'   .format(jobname       ,submit_script) )
        os.system( 'sed -i "s!BASH_OUTPUT!{0}!g" {1}'.format(bashoutput    ,submit_script) )
        os.system( 'sed -i "s!PP_DIR!{0}!g" {1}'     .format(peak_patch_dir,submit_script) )
 
        # Set cosmology parameters in Niagara jobscript
        os.system( 'sed -i "s!PKSC_FILE!{0}!g" {1}'    .format(self.catalogue_file,submit_script) )
        os.system( 'sed -i "s!NPZ_FILE!{0}!g" {1}'     .format(npz_file           ,submit_script) )
        os.system( 'sed -i "s!Z_CHI_TAB!{0}!g" {1}'    .format(z_chi_tab_file     ,submit_script) )
        os.system( 'sed -i "s!MASS_CUTOFF!{0}!g" {1}'  .format(mass_cutoff        ,submit_script) )
        os.system( 'sed -i "s!ZMIN!{0}!g" {1}'         .format(halo_zmin          ,submit_script) )
        os.system( 'sed -i "s!ZMAX!{0}!g" {1}'         .format(halo_zmax          ,submit_script) )
        os.system( 'sed -i "s!OMEGAM!{0}!g" {1}'       .format(self.Omx+self.OmB  ,submit_script) )
        os.system( 'sed -i "s!OMEGAB!{0}!g" {1}'       .format(self.OmB           ,submit_script) )
        os.system( 'sed -i "s!OMEGAL!{0}!g" {1}'       .format(self.Omvac         ,submit_script) )
        os.system( 'sed -i "s!LITTLEH!{0}!g" {1}'      .format(self.h             ,submit_script) )
        os.system( 'sed -i "s!SPECTRALINDEX!{0}!g" {1}'.format(self.ns            ,submit_script) )
        os.system( 'sed -i "s!SIGMA8!{0}!g" {1}'       .format(self.sigma8        ,submit_script) )

        # submit the job to a scheduler
        if submit == True:
            os.system( 'cd {0};sbatch ./{1}'.format( submit_dir, submit_script.split('/')[-1] ) )
        else:
            print( 'Job script ready to submit to Niagara scheduler with')
            print( 'cd {0}'.format( submit_dir ) )
            print( 'sbatch ./{0}'.format( submit_script.split('/')[-1] ) )



    def add_field(self,field_file=None,field_type=None,plane=None,intercept=0,buffers=False):
        # Load a Peak Patch initial conditions field.
        # 
        # Optional Parameters:
        # 
        #     field_file : string or None
        #         File name of Peak-Patch formatted initial conditions field. If None it will try to
        #         find a field of type field_type.
        # 
        #     field_type : string or None
        #         Option for the type of field to read. Allowed types are:
        #         'rhog'  - Gaussian overdensity field
        #         'rho'   - overdensity field
        #         'zetag' - Gaussian comoving curvature perturbation
        #         'zeta'  - comoving curvature perturbation
        #         'chi'   - transverse inflationary field
        # 
        #     plane : string or None
        #         We often plot slices of the Peak Patch initial conditions fields, so to save RAM 
        #         this option allows you to just read in a slice of the file, making use of byte
        #         offsets when reading the unformatted binary file. Select 'xy', 'xz', or 'yz' for
        #         respective planes. If None is passed, the whole field is read. Default is None.
        # 
        #     intercept : float
        #         If plane!=none, this is the value in Mpc of the intercept with plane <plane>.
        # 
        #     buffers : logical or int
        #         If False, buffers are excised, if True, the whole field is output. If buffers is
        #         an integer, a buffer of size buffers is excised (this allows you to trim the size
        #         of the field conveniently. Default is False.
        # 
        # USAGE - assuming you have loaded the PeakPatch repository:
        # >>> run = PeakPatch('/path/to/peak-patch-run/')
        # >>> run.add_field()
        # 
        # Note on formatting: Peak Patch initial conditions field files are formatted in a slightly peculiar way. 
        # 

        field_types = ['rhog','rho' ,'zetag','zeta'  ,'chi']
        field_names = ['rhog','Fvec','zetag','zetang','chi']
        def fname(ftype): return '{0}/fields/{1}_{2}'.format( self.run_dir, ftype, self.run_name )
        def field_t2n(ftype): return field_names[ field_types.index(ftype) ]

        # Assume Gaussian overdeensity field if no field type or field file specified
        if not field_type and not field_file:
            field_type,field_name='rhog','rhog'

        # Add support for passing "delta" or "matter" instead of "rho" as field type
        if   field_type.lower() in ('delta' ,'matter' ): field_type,field_name='rho' ,'rho'
        elif field_type.lower() in ('deltag','matterg'): field_type,field_name='rhog','rhog'

        # Try to read field type from file name if only field file given
        elif not field_type:
            for i in range(len(field_types)):
                if os.path.basename(field_file) == os.path.basename(fname(field_names[i])):
                    field_type = field_types[i]
                    field_name = field_names[i]
                    break

        # If field type given, determine which allowed type it refers to
        else:
            if type(field_type)==str:
                field_type = field_type.lower()
                if field_type not in field_types:
                    raise OSError(('Field type {0} not recognized. Allowed types are:\n\'{1}\', \''+
                        '{2}\', \'{3}\', \'{4}\', \'{5}\'').format( field_type,
                        field_types[0],field_types[1],field_types[2],field_types[3],field_types[4]))
                else:
                    field_name = field_t2n(field_type)
            elif type(field_type)==int and field_type in range(5):
                field_type = field_types[field_type]
                field_name = field_t2n(field_type)
            else:
                raise OSError(('Field type {0} not recognized. Allowed types are:\n\'{1}\', \'{2}'+
                    '\', \'{3}\', \'{4}\', \'{5}\'').format( field_type,
                    field_types[0],field_types[1],field_types[2],field_types[3],field_types[4]))

        # Auto-locate field file if none passed
        if not field_file:
            field_file = fname(field_name)
            if not os.path.exists(field_file):
                if field_type=='zeta' and os.path.exists(fname('zetag')): field_file=fname('zetag')
                elif os.path.exists(fname('Fvec')):                       field_file=fname('Fvec')
                else:
                    raise OSError(('file not found: {0}\nUnable to auto-locate field file.'
                                  ).format(field_file) )

        # Load entire field if no plane specified
        if not plane:
            field_in = open( field_file, 'rb' )
            field = np.reshape( np.fromfile( field_in, dtype=np.float32, count=-1 ),
                                (self.next,self.next,self.next), order='F' )
            #field = np.transpose( field, (2,1,0) )

            # Trim off the buffers if buffers is False
            if not buffers: field = field[ self.nbuff : -self.nbuff ,
                                           self.nbuff : -self.nbuff ,
                                           self.nbuff : -self.nbuff  ]

            # Trim off a buffer region of width buffer (in cells) if buffer is an integer
            elif isinstance(buffers,int): field = field[ buffers : -buffers ,
                                                         buffers : -buffers ,
                                                         buffers : -buffers  ]

        # Load slice of field
        else:

            # Check that plane is the correct type
            if type(plane) == int:
                if   plane == 0: plane = 'yz'
                elif plane == 1: plane = 'xz'
                elif plane == 2: plane = 'xy'
                else:
                    raise ValueError(('plane={0} not supported. If type(plane)==int, plane must be'+
                        ' 0, 1 or 2.').format(plane))
            elif type(plane) != str:
                raise TypeError('plane={0} not supported. plane must be of type int or str.'
                        .format(plane))

            field = np.array([])

            # To save RAM, read only the desired data from unformatted binary files using byte
            # offsets. Voxels that are adjacent in x have an offset of 4 (because 4B/32-b float).
            # Bytes offsets for voxels that are adjacent in y and z are:
            y_offset, z_offset = 4*self.next , 4*self.next**2

            # Read in slice of the x-y plane
            if 'z' not in plane.lower():

                # Determine lattice index of the intercept with the x-y plane
                if intercept   < self.cenz-self.boxsize/2 : intercept = 0
                elif intercept > self.cenz+self.boxsize/2 : intercept = self.neff
                else: intercept = int( (intercept-self.cenz+self.boxsize/2) / self.cellsize )

                # byte offset to first element of x-y plane with z intercept <intercept>
                if not buffers:
                    offset = 4*(self.next**2*(self.nbuff+intercept)+self.next*self.nbuff+self.nbuff)
                    for j in range(self.neff):
                        f = np.fromfile( field_file, dtype=np.float32, count=self.neff,
                                    offset=offset+j*y_offset )
                        field = np.concatenate( (field,f) )
                elif isinstance(buffers,int):
                    offset = 4*(self.next**2*(buffers   +intercept)+self.next*buffers   +buffers   )
                    for j in range(self.next-2*buffers):
                        f = np.fromfile( field_file, dtype=np.float32, count=self.next-2*buffers,
                                         offset=offset+j*y_offset )
                        field = np.concatenate( (field,f) )
                else:
                    offset = 4*(self.next**2*intercept)
                    for j in range(self.next):
                        f = np.fromfile( field_file, dtype=np.float32, count=self.next,
                                         offset=offset+j*y_offset )
                        field = np.concatenate( (field,f) )

            # Read in slice of the x-z plane
            elif 'y' not in plane.lower():

                # Determine lattice index of the intercept with the x-y plane
                if intercept   < self.ceny-self.boxsize/2 : intercept = 0
                elif intercept > self.ceny+self.boxsize/2 : intercept = self.neff
                else: intercept = int( (intercept-self.ceny+self.boxsize/2) / self.cellsize )

                # byte offset to first element of x-z plane with y intercept <intercept>
                if not buffers:
                    offset = 4*(self.next**2*self.nbuff+self.next*(self.nbuff+intercept)+self.nbuff)
                    #print(f"{offset}")
                    for k in range(self.neff):
                        f = np.fromfile( field_file, dtype=np.float32, count=self.neff,
                                         offset=offset+k*z_offset )
                        field = np.concatenate( (field,f) )
                elif isinstance(buffers,int):
                    offset = 4*(self.next**2*buffers   +self.next*(buffers   +intercept)+buffers   )
                    for k in range(self.next-2*buffers):
                        f = np.fromfile( field_file, dtype=np.float32, count=self.next-2*buffers,
                                         offset=offset+k*z_offset )
                        field = np.concatenate( (field,f) )
                else:
                    offset = 4*(self.next*intercept)
                    for k in range(self.next):
                        f = np.fromfile( field_file, dtype=np.float32, count=self.next,
                                         offset=offset+k*z_offset )
                        field = np.concatenate( (field,f) )

            # Read in slice of the y-z plane
            elif 'x' not in plane.lower():
 
                # Determine lattice index of the intercept with the x-y plane
                if intercept   < self.cenx-self.boxsize/2 : intercept = 0
                elif intercept > self.cenx+self.boxsize/2 : intercept = self.neff
                else: intercept = int( (intercept-self.cenx+self.boxsize/2) / self.cellsize )

                # byte offset to first element of x-z plane with y intercept <intercept>
                if not buffers:
                    offset = 4*(self.next**2*self.nbuff+self.next*self.nbuff+self.nbuff+intercept)
                    for k in range(self.neff):
                        for j in range(self.neff):
                            f = np.fromfile( field_file, dtype=np.float32, count=1,
                                             offset=offset+j*y_offset+k*z_offset )
                            field = np.concatenate( (field,f) )
                elif isinstance(buffers,int):
                    offset = 4*(self.next**2*buffers   +self.next*buffers   +buffers   +intercept)
                    for k in range(self.next-2*buffers):
                        for j in range(self.next-2*buffers):
                            f = np.fromfile( field_file, dtype=np.float32, count=1,
                                             offset=offset+j*y_offset+k*z_offset )
                            field = np.concatenate( (field,f) )
                else:
                    offset = 4*intercept
                    for k in range(self.next):
                        for j in range(self.next):
                            f = np.fromfile( field_file, dtype=np.float32, count=1,
                                             offset=offset+j*y_offset+k*z_offset )
                            field = np.concatenate( (field,f) )

            # If slice could not be recognized
            else:
                raise ValueError(('plane={0} not recognized. plane must be one of the following st'+
                    'rings\n[\'yz\',\'xz\',\'xy\'] (not case sensitive).').format(plane))

            if not buffers:
                field = np.reshape(field,(self.neff,self.neff),order='C')
            elif isinstance(buffers,int):
                field = np.reshape(field,(self.next-2*buffers,self.next-2*buffers),order='C')
            else:
                field = np.reshape(field,(self.next,self.next),order='C')

            # Transpose field to express as xz field with z on horizontal axis
            if 'y' not in plane.lower():
                field = field.T
            # By convention, all field slices are viewed as if from "above", in other words, the
            # normal vector to the plane being plotted (in this case the xz plane) oriented along
            # the Cartesian basis vector (in this case y) points to the observer. To achieve this,
            # the transpose is taken above must be taken here for the xz plane, but not for the yz
            # and xy planes. 

        # Set the bin edges for the discretised fields

        # If the argument `buffers' is False, then the default Peak Patch buffers are trimmed off,
        # and the field has physical size (self.boxsize)^3
        if not buffers:
            field_edges = np.array( [
                np.linspace( self.cenx-self.boxsize/2 , self.cenx+self.boxsize/2 , self.neff+1 ) ,
                np.linspace( self.ceny-self.boxsize/2 , self.ceny+self.boxsize/2 , self.neff+1 ) ,
                np.linspace( self.cenz-self.boxsize/2 , self.cenz+self.boxsize/2 , self.neff+1 ) ] )

        # If the argument `buffers' is an integer, then that many cells are trimed off of each edge
        # of the field, and the field has physical size (self.boxsize + 2*self.cellsize*buffers)^3
        elif isinstance(buffers,int):
            field_edges = np.array( [
                np.linspace( self.cenx-self.boxsize/2-self.buffersize+self.cellsize*buffers,
                             self.cenx+self.boxsize/2+self.buffersize-self.cellsize*buffers, 
                             self.next+1-2*buffers ),
                np.linspace( self.ceny-self.boxsize/2-self.buffersize+self.cellsize*buffers,
                             self.ceny+self.boxsize/2+self.buffersize-self.cellsize*buffers,
                             self.next+1-2*buffers ),
                np.linspace( self.cenz-self.boxsize/2-self.buffersize+self.cellsize*buffers,
                             self.cenz+self.boxsize/2+self.buffersize-self.cellsize*buffers,
                             self.next+1-2*buffers ) ]) 

        # If the argument `buffers' is True, then the no buffers are trimmed off and the field has
        # physical size (self.boxsize + 2*self.cellsize+self.nbuff)^3
        else:
            field_edges = np.array( [
                np.linspace( self.cenx-self.boxsize/2-self.buffersize , 
                             self.cenx+self.boxsize/2+self.buffersize , self.next+1 ) ,
                np.linspace( self.ceny-self.boxsize/2-self.buffersize , 
                             self.ceny+self.boxsize/2+self.buffersize , self.next+1 ) ,
                np.linspace( self.cenz-self.boxsize/2-self.buffersize , 
                             self.cenz+self.boxsize/2+self.buffersize , self.next+1 ) ] )

        # Assign field to correpsonding PeakPatch attribute
        if   field_type == field_types[0] : self.rhog  = [ field , plane , field_edges ]
        elif field_type == field_types[1] : self.rho   = [ field , plane , field_edges ]
        elif field_type == field_types[2] : self.zetag = [ field , plane , field_edges ]
        elif field_type == field_types[3] : self.zeta  = [ field , plane , field_edges ]
        elif field_type == field_types[4] : self.chi   = [ field , plane , field_edges ]



    def add_field_slice(self,field_file=None,field_type=None,plane='xz',intercept=0,buffers=False):
        # Convenience function that points to add_field, but assumes a slicing. See add_field()
        # above for full description of parameters.

        # Check that plane is correctly formatted
        if type(plane) != str:
            raise TypeError('plane must be a string. Note, this is a convenience function\nthat ca'+
                            'lls an essential function add_field(), and allows for less options\nt'+
                            'han its the essential function.')
        plane = plane.lower().strip()
        if plane not in ('xy','xz','yz'):
            raise ValueError('plane must be \'xy\' or \'xz\' or \'yz\'.')

        # Call essential function
        self.add_field( field_file = field_file , field_type = field_type , plane = plane ,
                        intercept  = intercept  , buffers    = buffers )



    def add_field_mean(self,field_type='rhog',overwrite=False):
        # Calculate the mean from initial conditions fields.
        # 
        # Optional Parameters:
        # 
        #     field_type : string or None
        #         Option for the type of field to read. Allowed types are:
        #         'rhog'  - Gaussian overdensity field
        #         'rho'   - overdensity field
        #         'rhong' - just the non-Gaussian part of the overdensity field
        #         'zetag' - Gaussian comoving curvature perturbation
        #         'zeta'  - comoving curvature perturbation
        #         'zetang'- just the non-Gaussian part of the zeta field
        #         'chi'   - transverse inflationary field
        #
        #     overwrite : bool
        #         if a mean has already been calculated, overwrite determines whether to recalculate
        #         it. if true, the existing value will be overwritten, if false it will not.
        # 
        # USAGE - assuming you have loaded the PeakPatch repository:
        # >>> run = PeakPatch('/path/to/peak-patch-run/')
        # >>> run.add_field()
        # >>> run.add_field_mean()

        # Supported values of field_type
        field_types = ['rhog','rho','rhong' , 'zetag','zeta','zetang'  ,'chi']

        # Assume Gaussian overdeensity field if no field type or field file specified
        if field_type not in field_types:
            warnings.warn('field_type "{0}" not recognised. Proceding with field_type="rhog".'
                .format(field_type) , Warning)
            field_type='rhog'

        # Check that the necessary fields have been read already
        if field_type in ['rhog','rhong'] and not hasattr(self,'rhog'):
            raise KeyError('Couldn\'t find field self.rhog, load it using add_field_mean().')
        elif field_type in ['rho','rhong'] and not hasattr(self,'rho'):
            raise KeyError('Couldn\'t find field self.rho, load it using add_field_mean().')
        elif field_type in ['zetag','zetang'] and not hasattr(self,'zetag'):
            raise KeyError('Couldn\'t find field self.zetag, load it using add_field_mean().')
        elif field_type in ['zeta','zetang'] and not hasattr(self,'zeta'):
            raise KeyError('Couldn\'t find field self.zeta, load it using add_field_mean().')
        elif field_type == 'chi' and not hasattr(self,'chi'):
            raise KeyError('Couldn\'t find field self.chi, load it using add_field_mean().')

        # Check that the mean hasn't already been calculated
        if ( not hasattr(self,field_type+'_mean') or 
             ( hasattr(self,field_type+'_mean') and overwrite==True ) ):

            # Calculate the mean
            if   field_type == 'rhog'   : self.rhog_mean   = np.mean( self.rhog[0]                 )
            elif field_type == 'rho'    : self.rho_mean    = np.mean( self.rho[0]                  )
            elif field_type == 'rhong'  : self.rhong_mean  = np.mean( self.rho[0]  - self.rhog[0]  )
            elif field_type == 'zetag'  : self.zetag_mean  = np.mean( self.zetag[0]                )
            elif field_type == 'zeta'   : self.zeta_mean   = np.mean( self.zeta[0]                 )
            elif field_type == 'zetang' : self.zetang_mean = np.mean( self.zeta[0] - self.zetag[0] )
            elif field_type == 'chi'    : self.chi_mean    = np.mean( self.chi[0]                  )

        else:
            print( 'Existing mean found. Not recalculating.' )



    def add_field_stdev(self,field_type='rhog',overwrite=False):
        # Calculate the standard deviation from initial conditions fields. Uses the other PeakPatch
        # method add_field_mean.
        # 
        # Optional Parameters:
        # 
        #     field_type : string or None
        #         Option for the type of field to read. Allowed types are:
        #         'rhog'  - Gaussian overdensity field
        #         'rho'   - overdensity field
        #         'rhong' - just the non-Gaussian part of the overdensity field
        #         'zetag' - Gaussian comoving curvature perturbation
        #         'zeta'  - comoving curvature perturbation
        #         'zetang'- just the non-Gaussian part of the zeta field
        #         'chi'   - transverse inflationary field
        # 
        #     overwrite : bool
        #         if a mean has already been calculated, overwrite determines whether to recalculate
        #         it. if true, the existing value will be overwritten, if false it will not.
        # 
        # USAGE - assuming you have loaded the PeakPatch repository:
        # >>> run = PeakPatch('/path/to/peak-patch-run/')
        # >>> run.add_field()
        # >>> run.add_field_stdev()

        # Supported values of field_type
        field_types = ['rhog','rho','rhong' , 'zetag','zeta','zetang'  ,'chi']

        # Calculate the mean if it has not already been calculated
        self.add_field_mean( field_type=field_type, overwrite=overwrite )

        # Check that the standard deviation hasn't already been calculated
        if ( not hasattr(self,field_type+'_stdev') or 
             ( hasattr(self,field_type+'_stdev') and overwrite==True ) ):

            # Calculate the mean
            if   field_type == 'rhog'   :
                self.rhog_stdev   = np.sqrt( np.mean(self.rhog[0]**2) - self.rhog_mean**2 )
            elif field_type == 'rho'    :
                self.rho_stdev    = np.sqrt( np.mean(self.rho [0]**2) - self.rho_mean **2 )
            elif field_type == 'rhong'  :
                self.rhong_stdev  = np.sqrt( np.mean((self.rho[0]-self.rhog[0])**2) 
                                             - self.rhong_mean**2 )
            elif field_type == 'zetag'  :
                self.zetag_stdev  = np.sqrt( np.mean(self.zetag[0]**2) - self.zetag_mean**2 )
            elif field_type == 'zeta'   :
                self.zeta_stdev   = np.sqrt( np.mean(self.zeta[0]**2) - self.zeta_mean**2 )
            elif field_type == 'zetang' :
                self.zetang_stdev = np.sqrt( np.mean((self.zeta[0]-self.zetag[0])**2)
                                             - self.zetang_mean**2 )
            elif field_type == 'chi'    :
                self.chi_stdev    = np.sqrt( np.mean(self.chi[0]**2) - self.chi_mean**2 )

        else:
            print( 'Existing standard deviation found. Not recalculating.' )



    def add_field_covering_factor(self,field_type='rhog',nsigma=3.,sigma=None,overwrite=False):
        # Calculate a covering factor to find the amount of a given IC field that is greater than a
        # specified number of standard deviations from the mean value. Uses the other PeakPatch
        # method add_field_stdev.
        # 
        # Optional Parameters:
        # 
        #     field_type : string or None
        #         Option for the type of field to read. Allowed types are:
        #         'rhog'  - Gaussian overdensity field
        #         'rho'   - overdensity field
        #         'rhong' - just the non-Gaussian part of the overdensity field
        #         'zetag' - Gaussian comoving curvature perturbation
        #         'zeta'  - comoving curvature perturbation
        #         'zetang'- just the non-Gaussian part of the zeta field
        #         'chi'   - transverse inflationary field
        # 
        #     nsigma : float
        #         The covering factor determines the fraction of the field that is greater than
        #         nsigma times the standard deviation of the field set by sigma. Default value is
        #         3.0.
        #
        #     sigma : float or str or None
        #         The standard deviation to use for the field, the covering determines what fraction
        #         of the field has a value greater than nsigma*sigma. The default value for sigma is
        #         None, in which case the standard deviation of the field is used by emplying the
        #         method self.add_field_stdev(). Other values are 'g' to use the standard deviation
        #         of the corresponding Gaussian field, 'nong' to use the standard deviation of the
        #         corresponding non-Gaussian field, or 'total' to use the standard deviation of the
        #         total field, or you can just pass it a floating point value.
        #
        #     overwrite : bool
        #         if a mean has already been calculated, overwrite determines whether to recalculate
        #         it. if true, the existing value will be overwritten, if false it will not.
        # 
        # USAGE - assuming you have loaded the PeakPatch repository:
        # >>> run = PeakPatch('/path/to/peak-patch-run/')
        # >>> run.add_field()
        # >>> run.add_field_covering_factor()

        # Figure out what sigma should be
        if sigma == None:
        
            # Determine the standard deviation if if hasn't been determined yet
            self.add_field_stdev(field_type=field_type,overwrite=overwrite)
            if   field_type == 'rhog'   : sigma = self.rhog_stdev
            elif field_type == 'rho'    : sigma = self.rho_stdev
            elif field_type == 'rhong'  : sigma = self.rhong_stdev
            elif field_type == 'zetag'  : sigma = self.zetag_stdev
            elif field_type == 'zeta'   : sigma = self.zeta_stdev
            elif field_type == 'zetang' : sigma = self.zetang_stdev
            elif field_type == 'chi'    : sigma = self.chi_stdev

        # If sigma is passed as a string
        elif type(sigma)==str:

            # Use standard deviation of the gaussian field
            if sigma.strip().lower()[0] == 'g':
                if   field_type in ['rhog','rho','rhong']:
                    self.add_field_stdev(field_type='rhog',overwrite=overwrite)
                    sigma = self.rhog_stdev
                elif field_type in ['zetag','zeta','zetang']:
                    self.add_field_stdev(field_type='zetag',overwrite=overwrite)
                    sigma = self.zetag_stdev
                elif field_type == 'chi':
                    self.add_field_stdev(field_type='chi',overwrite=overwrite)
                    sigma = self.chi_stdev

            # Use standard deviation of the non-Gaussian field
            elif sigma.strip().lower()[0] == 'n':
                if   field_type in ['rhog','rho','rhong']:
                    self.add_field_stdev(field_type='rhong',overwrite=overwrite)
                    sigma = self.rhong_stdev
                elif field_type in ['zetag','zeta','zetang']:
                    self.add_field_stdev(field_type='zetang',overwrite=overwrite)
                    sigma = self.zetang_stdev
                elif field_type == 'chi':
                    self.add_field_stdev(field_type='chi',overwrite=overwrite)
                    sigma = self.chi_stdev

            # Use standard deviation of the non-Gaussian field
            elif sigma.strip().lower()[0] == 't':
                if   field_type in ['rhog','rho','rhong']:
                    self.add_field_stdev(field_type='rho',overwrite=overwrite)
                    sigma = self.rho_stdev
                elif field_type in ['zetag','zeta','zetang']:
                    self.add_field_stdev(field_type='zeta',overwrite=overwrite)
                    sigma = self.zeta_stdev
                elif field_type == 'chi':
                    self.add_field_stdev(field_type='chi',overwrite=overwrite)
                    sigma = self.chi_stdev

            # Otherwise throw an error
            else: raise ValueError('sigma = "{0}" not recognised.'.format(sigma))

        # Otherwise, sigma must be a float
        elif type(sigma) != float:
            raise TypeError('sigma must be of type float, None, or a string of specified format')

        # Get the covering (an array of dimension equal to the field array's dimension 
        if   field_type == 'rhog'   :
            self.rhog_covering        = np.where( self.rhog[0] >= nsigma*sigma , 1., 0. )
            self.rhog_covering_factor = np.mean( self.rhog_covering )
        elif field_type == 'rho'    :
            self.rho_covering        = np.where( self.rho[0]  >= nsigma*sigma , 1., 0. )
            self.rho_covering_factor = np.mean( self.rho_covering )
        elif field_type == 'rhong'  :
            self.rhong_covering        = np.where( self.rho[0]-self.rhog[0] >= nsigma*sigma ,1.,0. )
            self.rhong_covering_factor = np.mean( self.rhong_covering )
        elif field_type == 'zetag'  :
            self.zetag_covering        = np.where( self.zetag[0] >= nsigma*sigma , 1.,0.)
            self.zetag_covering_factor = np.mean( self.zetag_covering )
        elif field_type == 'zeta'   :
            self.zeta_covering        = np.where( self.zeta[0]  >= nsigma*sigma , 1.,0.)
            self.zeta_covering_factor = np.mean( self.zeta_covering )
        elif field_type == 'zetang' :
            self.zetang_covering        = np.where( self.zeta[0]-self.zetag[0]>=nsigma*sigma,1.,0. )
            self.zetang_covering_factor = np.mean( self.zetang_covering )
        elif field_type == 'chi'    :
            self.chi_covering        = np.where( self.chi[0] >= nsigma*sigma , 1., 0. )
            self.chi_covering_factor = np.mean( self.chi_covering )



    def add_halo_covering_factor(self,mass_range=[None,None],overwrite=False):
        # Calculate a covering factor to find the amount of the simulation volume that is in a halo.
        # Rather can also constrain that to halos within a certain mass range.
        # Uses the other PeakPatch method add_halos.
        #    
        # Optional Parameters:
        #    
        #     mass_range : list of float or None
        #         Range of masses to consider for covering factor (in solar masses). If set to
        #         [None,None] then all halos are taken into account. Otherwise the first item is the
        #         lower bound and the second the upper bound.
        #    
        #     overwrite : bool
        #         if a mean has already been calculated, overwrite determines whether to recalculate
        #         it. if true, the existing value will be overwritten, if false it will not.
        #    
        # USAGE - assuming you have loaded the PeakPatch repository:
        # >>> run = PeakPatch('/path/to/peak-patch-run/')
        # >>> run.add_halos()
        # >>> run.add_halo_covering_factor()

        # Make sure stuff is properly formatted
        if type(mass_range) != list:
            raise TypeError('`mass_range\' must be type list.')
        else:
            if len(mass_range) != 2: raise IndexError('`mass_range\' must have length 2.')

        # Mass cut
        if mass_range == [None,None]: R_th = self.R_th
        elif mass_range[0] == None  : R_th = self.R_th[ np.argwhere( self.M < mass_range[1] ).T[0] ]
        elif mass_range[1] == None  : R_th = self.R_th[ np.argwhere( self.M > mass_range[0] ).T[0] ]
        else: R_th = self.R_th[ np.argwhere((self.M < mass_range[1])&(self.M>mass_range[0]) ).T[0] ]
        
        # Calculate covering factor
        self.lagrangian_halo_covering_factor = np.sum( 4/3*np.pi*R_th**3 ) / self.boxsize**3
        self.eulerian_halo_covering_factor   = self.lagrangian_halo_covering_factor/200   



    def add_healpix_fullsky_map(self,hp_file=None,response_function='tsz'):

            # Phyysical constants
            h     = const.h.to(   u.J * u.s ).value # Planck constant in joules * second
            k_B   = const.k_B.to( u.J / u.K ).value # Boltzmann constant in joules / kelvin

            # CMB temperature in kelvin
            T_CMB = 2.7255 # K
            def x(nu):    return h*nu/(k_B*T_CMB)
            def dT_y(nu): return T_CMB*(x(nu)*(np.exp(x(nu))+1)/(np.exp(x(nu))-1)-4)

            # Read HEALPix map from .fits file
            self.hp_fits_file = ('{0}/maps/{1}_{2}_{3}_hp.fits'.format(
                    self.run_dir, self.run_name, response_function, self.seed ) )
            if nix == 1:
                print('Healpy is currently not supported on nix')
            else:
                self.hp_map = hp.read_map( self.hp_fits_file )

            # Then mapview does a bunch of scaling i.e. by np.log10(map+1e-9) for tsz and cib as well as smotthing
            #     map = hp.smoothing( map , fwhm=np.radians(1./60.*15.) )
            # there's a note in mapview that says that this scaling of 15 arcminutes looks good for lensing convergence adn CIB. Presumably we could swap that out for a kwarg, maybe if you give like 'chord' or 'so' then it will calcualte the beam angle accordingly, or else it can take a floating point number or auto smooth to match what's done in mapview.

            # Note that you should probably have separate hp_map objects for each response function... can read response func from file name of .fits file...



    def add_flatsky_map(self,response_function=None,filein=None):
        if response_function==None and filein==None:
            raise ValueError('response_function and filein can\'t both be None.')
        elif filein==None:
            filein = '{0}/maps/{1}_kap_{2}_fs.map'.format( self.run_dir, self.run_name, self.seed )
        elif response_function==None:
            if   'tsz' in filein.lower(): response_function = 'tsz'
            elif 'ksz' in filein.lower(): response_function = 'ksz'
            elif 'kap' in filein.lower(): response_function = 'kap'
            elif 'cib' in filein.lower(): response_function = 'cib'

        read_file = open(filein,'rb')

        npix_RA = np.fromfile( read_file, dtype=np.int32, count=1 )[0]
        npix_dec = np.fromfile( read_file, dtype=np.int32, count=1 )[0]

        fov_RA  = float( np.fromfile( read_file, dtype=np.float32, count=1 )[0] /2/np.pi*3600 )
        fov_dec = float( np.fromfile( read_file, dtype=np.float32, count=1 )[0] /2/np.pi*3600 )

        fs_map = np.reshape( np.fromfile( read_file, dtype=np.float32, count=npix_RA*npix_dec ),
                     (npix_RA,npix_dec) )

        vars(self)['fs_map_'+response_function] = [ npix_RA,npix_dec,fov_RA,fov_dec,np.log10(fs_map+1e-9) ]



    def plot_healpix_fullsky_map(self,response_function='tsz'):
        # Actually, response function should be built into the map...

        # So mapview makes you a fig then populates with healpix map
        figsize_inch = 24,16
        fig = plt.figure(figsize=figsize_inch,dpi=dpi)
        if nix == 1:
            print('Healpy is currently not supported on nix')
        else:
            hp.mollview(self.hp_map, fig=fig.number, xsize=figsize_inch[0]*dpi, cmap='viridis', cbar=None, title="", min=vmin, max=vmax) 
        # Where vmin and vmax need to be set based on response function which should be readable from fits file



    def plot_flatsky_map(self,fig,ax,response_function='tsz', vmin=None, vmax=None, cmap=None ):

        ax.set_xlabel(r'RA  (deg)')
        ax.set_ylabel(r'dec (deg)')

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)

        ax.tick_params('both', length=10, width=4, which='major')
        ax.tick_params('both', length=6, width=3, which='minor')

# 41 #Set font sizes etc.
# 42 plt.rcParams.update({'axes.labelsize': 36, 'font.size': 36,
# 43                     'xtick.labelsize': 36,
# 44                     'ytick.labelsize': 36, 'axes.linewidth': 4,
# 45                     'xtick.labelsize': 36})
# 46 
# 47 font = {'family' : 'normal',
# 48         'weight' : 'normal',
# 49         'size'   : 24}
# 50 matplotlib.rc('font', **font)
# 51 
# 52 
# 53 #change font
# 54 plt.rc('text', usetex=True)
       # 55 plt.rc('font',**{'family':'serif','serif':['Palatino']})

        if   response_function=='tsz': cblab = 'log Compton-y'
        elif response_function=='ksz': cblab = r'$\Delta T_{\textrm kSZ}\ [\mu K]$'
        elif response_function=='kap': cblab = r'$\kappa$'
        elif response_function=='cib': cblab = 'log MJy/sr'
        elif response_function=='cmb': cblab = r'$\mu K$'

        if   vmin == None and response_function == 'tsz': vmin = -9.05
        elif vmin == None and response_function == 'ksz': vmin = -1.02e-5
        elif vmin == None and response_function == 'kap': vmin = 0.
        elif vmin == None and response_function == 'cib': vmin = -5.75
        elif vmin == None and response_function == 'cmb': vmin = -560
        elif vmin == False: vmin = np.min( vars(self)['fs_map_'+response_function][4] )

        if   vmax == None and response_function == 'tsz': vmax = -3.92
        elif vmax == None and response_function == 'ksz': vmax = 1.02e-5
        elif vmax == None and response_function == 'kap': vmax = 0.067
        elif vmax == None and response_function == 'cib': vmax = -1.25
        elif vmax == None and response_function == 'cmb': vmax = 560
        elif vmax == False: vmin = np.max( vars(self)['fs_map_'+response_function][4] )

        if   cmap == None and response_function == 'tsz': cmap = 'magma'
        elif cmap == None and response_function == 'ksz': cmap = 'coolwarm'
        elif cmap == None and response_function == 'kap': cmap = 'jet'
        elif cmap == None and response_function == 'cib': cmap = 'hot'
        elif cmap == None and response_function == 'cmb': cmap = planck_cmap()

        im = ax.imshow( vars(self)['fs_map_'+response_function][4],
                        vmin=vmin, vmax=vmax,cmap=cmap,
                        extent=[ -vars(self)['fs_map_'+response_function][2]/2,
                                  vars(self)['fs_map_'+response_function][2]/2,
                                 -vars(self)['fs_map_'+response_function][3]/2,
                                  vars(self)['fs_map_'+response_function][3]/2 ],
                        origin='upper' )

        # Add colorbar
        cbaxes = fig.add_axes([0.925,0.1,0.05,0.8])
        cb     = plt.colorbar(im, cax=cbaxes)
        cb.set_label(cblab)

        # 
        return fig,ax



    def add_IC_power_spectra(self,spectrum_file=None,spectrum_type=None):
        # Loads specified initial conditions power spectrum or transfer function for a run
        # associated with a Peak Patch object.
        # 
        # Optional Parameters:
        # 
        #     spectrum_file : string or None
        #         File name of Peak-Patch formatted initial conditions power spectrum. If None it
        #         will try to find a spectrum file of type spectrum_type.
        # 
        #     spectrum_type : string or None
        #         Option for the type of power spectrum or transfer function to read. Allowed types
        #         are:
        #              'P_deltadelta'                 - z=0 matter density power spectrum
        #              'P_zetazeta'                   - primordial zeta power spectrum
        #              'T_zeta2delta'                 - zeta to delta transfer function
        #              'P_chaoticbilliardspreheating' - preheating power spectrum with chaotic
        #                                               billiards effects
        #              'P_ping'                       - inflationary power spectrum with primordial
        #                                               intermittent non-Gaussianity caused by an
        #                                               instability in the potential surface
        # 
        #     All power spectra are in dimensional form
        #         P_ff(k) = 2\pi^2/k^3 \frac{ d\langle \tilde{f}^2(k) \rangle }{ d ln k }
        #     with units of Mpc^3 assuming h from Planck 2018 (this should be updated to read in the
        #     value of h from the parameter files but hasn't been as of September 2023). All
        #     transfer functions are in the form
        #         T_f2g(k) = k^{-2} \sqrt{ P_gg(k) / P_ff(k) }
        #
        # USAGE - assuming you have loaded the PeakPatch repository:
        # >>> run = PeakPatch('/path/to/peak-patch-run/')
        # >>> run.add_add_IC_power_spectra()
        # 
        # Notes on formatting
        # Peak Patch initial conditions fields are generated from power spectrum and transfer
        # function files that are formatted in a slightly peculiar way. Here are a few important
        # notes on the formats used:
        # 
        # - Units
        #   Peak Patch always uses distance units of Mpc, not Mpc/h. This is a source for confusion
        #   because the LSS community very often uses Mpc/h. You can always divide out the h in the
        #   parameter file, but keep in mind that all calculations of dynamics are assuming this
        #   value of h so by just dividing h out, you're still not truly agnostic of its value.
        # 
        # - Power spectra files
        #   Files tabulating power spectra and transfer functions have some formatting querks as
        #   well, these files are generated using CAMB and the function powerspectrum_create()
        #   below, which saves them to the `peakpatch/tables/' directory. They are formatted with 4
        #   columns as follows
        #       column 1: wave number, k, [Mpc^{-1}]
        #       column 2: linear matter power spectrum at z=0 in units of Mpc^3
        #                 P_peakpatch = (2\pi)^{-3} P_{ff}(k)]
        #       column 3: transfer function to primordial zeta power in units of Mpc^3
        #                 T_peakpatch = (2\pi)^{-3/2} T_{\zeta \to f}(k)
        #       column 4: linear theory power spectrum of primordial field sourcing non-
        #                 Gaussianities, P_{\chi\chi}(k), [Mpc^3] in the same format as the density
        #                 power.
        # where
        #     f = \delta \bar{\rho} = \rho - \bar{\rho}
        # is the density field minus it's average (or you can look at this as the dimensionless
        # overdensity paramter $\delta$ times the average density $\bar{\rho}$), 
        #     P_{ff}(k) = \langle | \tilde{f}(\mathbf{k}) |^2 \rangle
        # is the linear matter power spectrum in the typical form, taken as the Fourier transform of
        # $f$ averaged over spherical shells of each wavenumber $k$ (note that the Peak Patch linear
        # matter power spectrum has an extra factor of $(2\pi)^{-3}$ multiplying it), and
        #     T_{\zeta to f} = \sqrt{ P_{ff}(k) / P_{\zeta\zeta}(k) }
        # is the transfer function that relates the linear matter and primordial $\zeta$ power
        # spectra. Finally, the field $P_{\chi\chi}$ is the power spectrum for a second early 
        # universe field, either a transverse field during multi-field inflation or an effective
        # field during reheating.

        # Check that spectrum_type is nonetype or an allowed string
        if spectrum_type == None:
            if spectrum_file == None:
                raise Warning('assuming matter power spectrum.')
                spectrum_type = 'P_deltadelta'
                spectrum_file = '{0}/tables/{1}'.format( self.run_dir, self.pkfile )
            else:
                if self.pkfile in spectrum_file:
                    raise Warning('assuming matter power spectrum.')
                    spectrum_type = 'P_deltadelta'
                    spectrum_file = '{0}/tables/{1}'.format( self.run_dir, self.pkfile )
                else:
                    spectrum_file = '{0}/tables/{1}'.format( self.run_dir, 'p_chichi.dat' )

        elif isinstance(spectrum_type,str):
            spectrum_type = ( spectrum_type.replace(' ','').replace('\n','').replace('\t','')
                                           .replace('_','').lower() )
            if spectrum_type in [ 'p','pm','pk','pofk','pmofk','plin','pdelta','pdeltadelta',
                                  'plinear','linearmatterpower' ]:
                spectrum_type='P_deltadelta'
            elif spectrum_type in [ 'pzetazeta' ]:
                spectrum_type='P_zetazeta'
            elif spectrum_type in [ 'tzeta2delta' ]:
                spectrum_type='T_zeta2delta'
            elif spectrum_type in [ 'pchaoticbilliardspreheating' ]:
                spectrum_type='P_chaoticbilliardspreheating'
            elif spectrum_type in [ 'ping','pping','pchichi','spike','pspike' ]:
                spectrum_type = 'P_ping'


        # Check that spectrum_file is nonetype or a string pointing to a file that exists
        if spectrum_file == None:
            if ( spectrum_type=='P_deltadelta' or spectrum_type=='P_zetazeta' or 
                 spectrum_type=='T_zeta2delta' or spectrum_type=='P_chaoticbilliardspreheating' ):
                spectrum_file = '{0}/tables/{1}'.format( self.run_dir, self.pkfile )
            elif ( spectrum_type=='P_ping' ):
                if self.NonGauss < 11:
                    spectrum_file = '{0}/tables/{1}'.format( self.run_dir, 'p_chichi.dat' )
                else:
                    H = [ 'i', 'p', 'e' ]
                    spectrum_file = [ '{0}/tables/p_chichi_{1}.dat'
                                      .format(self.run_dir,H[j]) for j in range(3) ]
            else:
                raise ValueError('spectrum_type '+sepectrum_type+' not recognised.\nUnable to autom'
                    +'atically asign spectrum_file.')
            if self.NonGauss < 11 and not os.path.exists(spectrum_file):
                raise ValueError('Could not locate spectrum_file\n'+spectrum_file+'\nUnable to auto'
                    +'matically asign spectrum_file.')
            elif self.NonGauss == 11:
                for f in spectrum_file:
                    if not os.path.exists(f):
                        raise ValueError('Could not locate spectrum_file\n'+spectrum_file+'\nUnabl'+
                                'e to automatically asign spectrum_file.')
        elif isinstance(spectrum_file,str):
            if not os.path.exists(spectrum_file):
                raise ValueError('spectrum_file must be the path to a file that exists.')
        else:
            raise TypeError('spectrum_file must be string or nonetype.')


        # Read in power spectra
        if ( spectrum_type == 'P_deltadelta' or spectrum_type=='P_zetazeta' or
             spectrum_type=='T_zeta2delta' or spectrum_type=='P_chaoticbilliardspreheating' ):

            # Load spectrum file
            spectra_data = np.loadtxt( spectrum_file, delimiter=' ' )
            
            # Save matter power spectrum
            if spectrum_type=='P_deltadelta':
                self.k_deltadelta = spectra_data[:,0]
                self.P_deltadelta = spectra_data[:,1]

            # Save primordial Gaussian zeta power spectrum
            elif spectrum_type=='P_zetazeta':
                self.k_zetazeta = spectra_data[:,0]
                self.P_zetazeta = spectra_data[:,1] * spectra_data[:,2]**-2.

            # Save matter power transfer function
            elif spectrum_type=='T_zeta2delta':
                self.k_zeta2delta = spectra_data[:,0]
                self.T_zeta2delta = spectra_data[:,2] * spectra_data[:,0]**-2.

            # Save primordial power spectrum from chaotic billiards during preheating
            elif spectrum_type=='P_chaoticbilliardspreheating':
                self.k_chaoticbilliardspreheating = spectra_data[:,0]
                self.P_chaoticbilliardspreheating = spectra_data[:,3]
            
        elif spectrum_type == 'P_ping' and self.NonGauss < 11:

            # Length of P(k) table = file size in B / (4 B/32-b float) / 2 columns
            count = int( os.path.getsize(spectrum_file)/4/2 )

            # Read in the PING power spectrum P_chichi(k)
            p = np.reshape( np.fromfile(spectrum_file,dtype=np.float32,count=2*count) , (count,2) )
            self.k_ping        = p[:,0]
            self.p_chichi_ping = p[:,1]

        elif spectrum_type == 'P_ping':

            # Length of P(k) table = file size in B / (4 B/32-b float) / 2 columns
            counts = [ int( os.path.getsize(spectrum_file[j])/4/2 ) for j in range(3) ]

            # Read in the PING power spectra
            p_i, p_p, p_e = [ np.reshape( np.fromfile( spectrum_file[j],
                                                       dtype=np.float32,
                                                       count=2*counts[j]
                                                     ), (counts[j],2) ) for j in range(3) ]

            # Add powers to PeakPatch object
            self.k_i_ping , self.p_chichi_i_ping = p_i[:,0] , p_i[:,1]
            self.k_p_ping , self.p_chichi_p_ping = p_p[:,0] , p_p[:,1]
            self.k_e_ping , self.p_chichi_e_ping = p_e[:,0] , p_e[:,1]



    ################################################################################################
    ###                                                                                          ###
    ###                      POST-PPROCESSING Peak Patch AND WebSky RESULTS                      ###
    ###                                                                                          ###
    ################################################################################################
    # The functions in this section are used for post processing the Peak Patch and WebSky data
    # prodcuts loaded into a PeakPatch object using the load functions above. Included functions
    # are:
    #
    # Halo mass functions and histograms:
    # 
    #     halo_hist() : make histograms of Dark Matter halos catalogued by Peak Patch
    #     
    #     halo_hist2d() : make 2d histograms of DM halos, analagous to slices of ICs fields
    #     
    #     halo_hist3d() : make 3d histograms of DM halos, analagous to ICs fields
    # 
    #     hmf() : make a halo mass function of various forms
    # 
    #     hmzf() : make a 2D halo function of mass and redshift
    # 
    # Halo correlation functions
    # 
    # ICs fields:
    # 
    # Power spectra:
    # 
    #     get_power_spectrum() : calculates the power spectrum from a given field
    # 
    # WebSky maps:
    # 
    # WebSky power spectra:
    # 
    # See bolow for more information on their usage.

    ################################################################################################
    ###   Halo Mass Functions and Histograms   #####################################################

    # Make a histogram of halos binned by mass
    # def halo_hist(self, ... )



    def halo_hist2d(self, axis=0, bins=None, coord='e', bin_lim=None, mass_weighted=False,
            binaxes=['lin','lin'], smoothing=None, smoothing_scale=None):
        # Make a 2D histogram of halos added to PeakPatch object with self.add_halos(...)
        #
        # Optional parameters:
        # 
        #     axis : int or str
        #         The axis to project along to make a 2D histogram from halos distributed in 3D.
        #         Accpets 'x', 'y', 'z' or 0, 1, 2. Default is 0.
        # 
        #     bins : int or [int,int] or NumPy array or [NumPy array,NumPy array] or None
        #         Number of bins if type(bins)==int, the number of bins along each axis if bins is a
        #         list of int, symmetric bin edges if bins is a NumPy array or assymetric bins if
        #         bins is a list of two NumPy arrays. If bins==None, it will try to automatically
        #         bin at the lattice scale of Peak Patch initial conditions fields. Default is None.
        #
        #     coord : str
        #         Whether to bin by Eulerian postions (set to 'e', 'euler', 'eulerian') or
        #         Lagrangian positions (set to 'l', 'lagrange', 'lagrangian'). Not case sensitive.
        #         Default is 'e'.
        # 
        #     bin_lim : str or ndarray or None
        #         The obious choice for limits of the bins is the to run over the comoving volume of
        #         the box. In Peak Patch parlance, this would be x ∈ cenx ± boxsize/2 (and similarly
        #         for y and z). However, if the we're binning the Eulerian halo coordinates, some
        #         halos may move beyond this region, so this option allows you to choose how things
        #         are binned. If bin_lim is a string (not case-sensitive) set to 'l', 'lagrange' or
        #         'lagrangian' then the bounds are lagrangian bounds of the lagrangian space volume.
        #         If bin_lim is set to 'e', 'euler' or 'eulerian', then bins are set to the smallest
        #         and largest Eulerian halo positions along each axis. If bin_lim is set to None, it
        #         will make a best judgement, setting it to 'l' if coord=='l' or the min/max of the
        #         lagrangian minima/maxima and Euleran halo position minima/maxima. If bin_lim is an
        #         array of shape (2,) or (2,2) then its elements will be taken as the bounds.
        # 
        #     mass_weighted : logical
        #         Binning weighted by halo mass if True, binning is just counts if False. Default is
        #         False.
        # 
        #     binaxes : [ str , str ]
        #         if bins is passed as an int or [int,int] then binaxes determines whether bins
        #         should be linearly or logarithmically spaced. binaxes must be None or a list of
        #         two strings, with the allowed strings being 'lin' for linear spacing and 'log' for
        #         logarithmic spacing.
        # 
        #     smoothing : None or str
        #         Because halo hists just plot halo centres, so for larger halos, histograms end
        #         up being too sharp. In principle, this should be addressed by applying to each
        #         halo an effective width, but this is computationally expensive. Smoothing provides
        #         a quick approximate way of taking care of this, by smoothing the field at the
        #         scale of the largest halos. Smoothing can either be set to None for no smoothing,
        #         'gauss' for smoothing with a Gaussian filter kernel, or 'top-hat' for smoothing
        #         with a top-hat filter kernel.
        # 
        #     smoothing_scale : None or float
        #         The smoothing scale to use if smoothing is not None. This can either be a float,
        #         representing a length scale in Mpc, or None to calculate automatically. The auto-
        #         matically calculated method smooths at the scale of the largest halos. In Lagran-
        #         gian space, this is just R_th,max. In eulerian space, this will be divided by
        #         200^1/3 to represent the virialised halos. For Gaussian filter, a factor of 2 is
        #         included (e.g. see the last paragraph in section 3.2 of Clusters and the Theory
        #         of the Cosmic Web: https://www.cita.utoronto.ca/~bond/cweb/wb1_fulltext.pdf).
        # 
        # Returns:
        #     H : ndarray, shape(nx,ny)
        #         Histogram object.
        # 
        #     xedges : ndarray, shape(nx)
        #         Array contianing the bin edges along the first axis of H.
        # 
        #     yedges : ndarray, shape(ny)
        #         Array containing the bin edges along the second axis of H.

        # Check the axis along which to project halos into the plane of the 2D histogram
        if type(axis)==str:
            if 'x' in axis.lower(): axis=0
            if 'y' in axis.lower(): axis=1
            if 'z' in axis.lower(): axis=2
        if axis not in range(3): raise ValueError(('axis {0} not recognized, axis must be one of in'
                                                  +' \{0,1,2\}').format(axis))

        # Indices of the two axes not being projected with index[0] < index[1] these are ordered in
        # the same sense as Peak Patch IC fields, such that (axis,index[0],index[1]) is an even
        # permutation of (0,1,2)
        index = [ (axis-1)%3, (axis+1)%3 ]
        index = index[ :: levi_civita( 3, axis, index[0], index[1] ) ]

        # Checks coord set to Eulerian or Lagrangian histogram
        if type(coord) != str: raise TypeError('coord must be of type str.')

        # Define pointer(-ish) to Eulerian or Lagrangian positions of halos
        coord = coord.lower()
        if   coord=='e' or coord=='euler'    or coord=='eulerian'   :
            xyz = [np.ndarray.view(self.x),np.ndarray.view(self.y),np.ndarray.view(self.z)]
            x,y = np.ndarray.view(xyz[index[0]]) , np.ndarray.view(xyz[index[1]])
        elif coord=='l' or coord=='lagrange' or coord=='lagrangian' :
            xyz = [np.ndarray.view(self.xL),np.ndarray.view(self.yL),np.ndarray.view(self.zL)]
            x,y = np.ndarray.view(xyz[index[0]]) , np.ndarray.view(xyz[index[1]])
        else: raise ValueError('coord must be set to \'e\' or \'l\'.')

        # Sets weighting by masss
        if type(mass_weighted) != bool: raise TypeError('mass_weighted must be of type bool.')
        else:
            self.mass_weighted2d = mass_weighted
            if mass_weighted == True : weights = self.M
            else                     : weights = None

        # Check that binning is done correctly, binning can be linear or log
        if type(binaxes) != list: raise TypeError('binaxes must be a list.')
        elif len(binaxes) != 2  : raise ValueError('binaxes must be a list of length 2.')
        elif type(binaxes[0]) != str or type(binaxes[1]) != str:
            raise TypeError('binaxes must be list of stirngs.')
        else:
            binaxes = [ binaxes[0].lower() , binaxes[1].lower() ]
            if not ( binaxes[0] in ['lin','log'] and binaxes[1] in ['lin','log'] ):
                raise ValueError('Elements of binaxes must be either \'lin\' or \'log\'.')

        # If NumPy array(s) passed in bins, this supercedes binaxes
        if type(bins) != np.ndarray and type(bins) != int and type(bins) != list and bins != None:
            raise TypeError('bins must be of type int, np.ndarray, list or NoneType.')
        elif type(bins) == np.ndarray: binaxes=None
        elif type(bins) == list:
            if len(bins) != 2: raise ValueError('if bins is type List, it must have length 2.')
            elif type(bins[0]) != type(bins[1]):
                raise TypeError('bins[0] and bins[1] must have same type.')
            elif type(bins[0]) != np.ndarray and type(bins) != int:
                raise TypeError('bins must be [int,int] or [np.ndarray,np.ndarray].')
            elif type(bins[0]) == np.ndarray: binaxes=None

        # Check that bin_lim is the right format
        if type(bin_lim) == str:
            if   bin_lim.strip().lower()[0] == 'e': bin_lim='e'
            elif bin_lim.strip().lower()[0] == 'l': bin_lim='l'
            else: raise ValueError('If bin_lim is of type str, it must be a recognized value.')
        elif type(bin_lim) == np.ndarray and bin_lim.shape != (2,2) and bin_lim.shape != (2,):
            raise ValueError('If bin_lim is of type np.ndarray, it must have shape (2,) or (2,2).')
        elif type(bin_lim) == list and not (
                ( type(bin_lim[0]) == list and (
                    isinstance(bin_lim[0][0],float) or isinstance(bin_lim[0][0],int) or
                    isinstance(bin_lim[0][0],np.floating) or isinstance(bin_lim[0][0],np.integer) )
                ) or isinstance(bin_lim[0],int) or isinstance(bin_lim[0],float) or
                isinstance(bin_lim[0],np.floating) or isinstance(bin_lim[0],np.integer) ):
            raise TypeError('Listen, that\'s not the way this works. You can\'t just expect me to b'
                           +'e able to parse a\nrandom list like that.')

        # Set bins to bin edges unless bins=[int,int] and binaxes!=['lin','lin']
        if binaxes and not ( bins and binaxes == ['lin','lin'] ):
            if bin_lim == None or bin_lim=='l' or bin_lim=='e':
                cenxyz = np.array([self.cenx,self.ceny,self.cenz])
                if bin_lim == None:
                    rng = [ [ max( self.halo_lim[j,0],cenxyz[j]-self.boxsize/2 ) ,
                              min( self.halo_lim[j,1],cenxyz[j]+self.boxsize/2 ) ] for j in index ]
                elif bin_lim == 'l':
                    rng = [[cenxyz[j]-self.boxsize/2,cenxyz[j]+self.boxsize/2] for j in index ]
                else:
                    rng = [ [ self.halo_lim[j,0],self.halo_lim[j,1] ] for j in index ]
            elif isinstance(bin_lim,list):
                if isinstance(bin_lim[0],list): rng = bin_lim
                else                          : rng = [bin_lim,bin_lim]
            elif isinstance(bin_lim,np.ndarray):
                if bin_lim.shape == (2,): rng = [ bin_lim[0] , bin_lim[1] ]*2
                else                    : rng = [[bin_lim[i,j]for j in range(2)]for i in range(2)]
            N = [np.ceil((rng[j][1]-rng[j][0])/self.cellsize).astype(int) for j in range(2)]

            # If no bins are given, use the lattice scale from ICs fields
            if (not bins) and binaxes == ['lin','lin'] : bins = N

            # If bin axes are not both linear, set bins to bin edges
            else:
                if bins: N = bins
                bins = [np.zeros(0),np.zeros(0)]
                for j in range(2):
                    if binaxes[j]=='lin': bins[j] = np.linspace( rng[j,0] , rng[j,1] , N+1 )
                    else: bins[j] = np.logspace( np.log10(rng[j,0]) , np.log10(rng[j,1]) , N+1 )

        # Make histogram
        hist2d = np.histogram2d( x, y, bins=bins, weights=weights )

        # Smoothing histogram
        if smoothing is not None:

            # Top-hat filter kernel
            if smoothing == 'top-hat':

                # If no smoothing scale specified, calculate automatically
                if smoothing_scale is None:

                    # Maximum halo top-hat filter radius
                    R_th_max = np.max(self.R_th)

                    # If plotting in Eulerian coordinates, smooth by max halo top-hat filter radius
                    # divided by the cube root of critical density factor for virialisation in this
                    # Peak Patch run
                    if coord=='e':
                        smoothing_scale = R_th_max*self.dcrit**(-1./3.)

                    # If plotting in Lagrangian coordinates, just smooth with max halo top-hat
                    # filter radius
                    else:
                        smoothing_scale = R_th_max

                raise ValueError('I haven\'t actually finished writing this function yet, use smoot'
                        +'hing=\'gauss\'. Eventually it should involve e.g. fourier_top_hat_filter_'
                        +'2D() defined below.')

            # Gaussian filter kernel
            else:
                from scipy.ndimage import gaussian_filter

                # If no smoothing scale specified, calculate automatically
                if smoothing_scale is None:

                    # Maximum halo top-hat filter radius
                    R_th_max = np.max(self.R_th)

                    # If plotting in Eulerian coordinates, smooth by max halo top-hat filter radius
                    # divided by the cube root of critical density factor for virialisation in this
                    # Peak Patch run
                    if coord=='e':
                        smoothing_scale = 2*R_th_max*self.dcrit**(-1./3.)

                    # If plotting in Lagrangian coordinates, just smooth with max halo top-hat
                    # filter radius
                    else:
                        smoothing_scale = 2*R_th_max

                print(smoothing_scale)

                # Convert smoothing scale to pixels
                smoothing_scale /= ( (   hist2d[1][1] - hist2d[1][0]
                                       + hist2d[2][1] - hist2d[2][0] )/2 )

                print(smoothing_scale)

                # Apply a Gaussian filter
                self.hist2d = [ gaussian_filter( hist2d[0] , smoothing_scale ),
                                hist2d[1],
                                hist2d[2] ]

        else:
            self.hist2d = [ hist2d[0] , hist2d[1] , hist2d[2] ]

        return self.hist2d



#    # Generate halo continuous field from added halo catalogue
#    def halo_hist3D(self, other stuff):
#    # note: include option to drop buffers



    def hmf(self, bins=None, bin_axis='log', coord='e', hmf_type='dn' ):
        # Make a 2D histogram of halos added to PeakPatch object with self.add_halos(...)
        #
        # Optional parameters:
        # 
        #     bins : int or NumPy array or [NumPy array,NumPy array] or None
        #         Number of bins if type(bins)==int, the number of bins along each axis if bins is a
        #         list of int, symmetric bin edges if bins is a NumPy array or assymetric bins if
        #         bins is a list of two NumPy arrays. If bins==None, it will try to automatically
        #         bin at the lattice scale of Peak Patch initial conditions fields. Default is None.
        #
        #     bin_axis : str
        #         'log' for log-spaced mass bins, 'lin' for linearly spaced mass bins.
        # 
        #     coord : str
        #         Whether to bin by Eulerian postions (set to 'e', 'euler', 'eulerian') or
        #         Lagrangian positions (set to 'l', 'lagrange', 'lagrangian'). Not case sensitive.
        #         Default is 'e'.
        # 
        #     hmf_type : str
        #         The type of halo mass function to output. We have lots of options because there
        #         are lots of different ways that people do these. The options are:
        #             'dndm'    or 'dn/dm'    - dN / dM
        #             'dndlnm'  or 'dn/dlnm'  - dN / d ln( M )
        #             'dndlogm' or 'dn/dlogm' - dN / d log( M )
        #             'dn'      or 'hist'     - dN
        #             'pdf'                   - dN / Ntot
        #             'srel'                  - relative entropy = pdf * some const.
        #         Where dN is the count of halos in a bin, and dM is the width of a mass bin.
        # 
        # Returns:
        #     H : ndarray, shape (n)
        #         Histogram of halos per mass bin.
        # 
        #     bin_edges : ndarray, shape(n+1)
        #         Array containing the bin edges of histogram H.

        # Mass range of halo catalogue
        Mmin , Mmax = np.min(self.M) , np.max(self.M)

        # Check that bin_axis is in the correct form
        if type(bin_axis) != str:
            raise TypeError('bin_axes must be a string.')
        else:
            bin_axis = bin_axis.strip().lower()
            if bin_axis != 'log' and bin_axis != 'lin':
                raise ValueError('bin_axes must be \'lin\' or \'log\'.')

        # Check that hmf_type is in the correct form
        if type(hmf_type) != str:
            raise TypeError('hmf_type must be a string.')
        else:
            hmf_type = hmf_type.strip().lower()
            if   ( hmf_type == 'dn'         or hmf_type == 'histogram'  or
                   hmf_type == 'histogramme'                                   ): hmf_type = 'dn'
            elif ( hmf_type == 'srelative'  or hmf_type == 's_rel'      or
                   hmf_type == 's_relative'                                    ): hmf_type = 'srel'
            elif ( hmf_type == 'dn/dm'      or hmf_type == 'dn_dm'             ): hmf_type = 'dndm'
            elif ( hmf_type == 'dn/dlnm'    or hmf_type == 'dn/dln(m)'  or
                   hmf_type == 'dn_dlnm'    or hmf_type == 'dn_dln(m)'         ): hmf_type = 'dndlnm'
            elif ( hmf_type == 'dn/dlogm'   or hmf_type == 'dn/dlog(m)' or
                   hmf_type == 'dn_dlogm'   or hmf_type == 'dn_dlog(m)'        ): hmf_type = 'dndlogm'
            elif ( hmf_type != 'dndm' and hmf_type != 'dndlnm' and hmf_type != 'dndlogm' and
                   hmf_type != 'pdf'  and hmf_type != 'dn'     and hmf_type != 'srel' ):
                raise ValueError("hmf_type must be one of 'dndm', 'dndlnm', 'dndlogm', 'pdf', 'dn'"+
                        ", and 'srel' (not case\nsensitive).")

        # Mass binning
        if bins is None:
            bins = int( (self.neff-1)/2 )

        # If bins is a floating point number just convert to an integer
        if isinstance( bins , float ):
            print('Warning: bins must be an integer or array, not a float. Converting to integer.')
            bins = int(bins)

        # If bins is an integer, set bins to bin edges with scaling based on bin_axis
        elif isinstance( bins, int ):
            if bin_axis == 'log':
                bins = np.logspace( np.log10(Mmin) , np.log10(Mmax) , bins+1 )
            else:
                bins = np.linspace( Mmin , Mmax , bins+1 )

        # Define bin centres
        self.bin_centres = ( bins[1:] + bins[:-1] )/2

        # Bin halos by mass where each bin j centred at M_j has dn[j] = dN(M_j) halos with a mass
        # range M \in [ bin_edges[j] , bin_edges[j+1] )
        dn , self.bin_edges = np.histogram( self.M , bins=bins )

        # For a probability density function dN(M_j) / N
        if hmf_type == 'pdf':
            self.pdf = dn / self.N_halos
            return self.pdf , self.bin_edges

        # For a regular histogram dN(M_j)
        if hmf_type == 'dn':
            self.dn = dn
            return self.dn , self.bin_edges

        # For halos per mass dN(M_j) / dM_j, where 2dM_j = bin_edges[j+1]-bin_edges[j]
        if hmf_type == 'dndm':
            dm = ( self.bin_edges[1:] - self.bin_edges[:-1] )/2
            self.dndm = dn / dm
            return self.dndm , self.bin_edges

        # For halos per ln mass dN(M_j) / dlnM_j, where dlnM_j = ln(bin_edges[j+1])-ln(bin_edges[j])
        if hmf_type == 'dndlnm':
            dlnm = ( np.log( self.bin_edges[1:] / self.bin_edges[:-1] ) )/2
            self.dndlnm = dn / dlnm
            return self.dndlnm , self.bin_edges

        # For halos per log mass dN(M_j) / dlogM_j, where dlogM_j defined similarly to dlnM_j
        if hmf_type == 'dndlogm':
            dlogm = ( np.log10( self.bin_edges[1:] / self.bin_edges[:-1] ) )/2
            self.dndlogm = dn / dlnm
            return self.dndlogm , self.bin_edges

        # For relative entropy S_rel = M dN(M_j) / dlnM_j
        if hmf_type == 'srel':
            dlnM = ( np.log( self.bin_edges[1:] / self.bin_edges[:-1] ) )/2
            self.srel = self.bin_centres * dn / dlnM
            return self.srel , self.bin_edges



    ################################################################################################
    ###   Halo correlation functions   #############################################################



    ################################################################################################
    ###   ICs fields   #############################################################################



    ################################################################################################
    ###   Power spectra   ##########################################################################

    def get_power_spectrum(self, field_file=None, field_type=None, method='total', overwrite=False,
                           verbose=True):
        # Run fortran scripts to read the power spectrum from a 3D field.
        # 
        # Optional parameters:
        # 
        #     field_file : str or None
        #         The path to a file containing a Peak-Patch-initital-conditions-formatted
        #         cosmological field or None to automatically select a field based on the field_type
        #         parameter.
        #
        #     field_type : str or None
        #         The type of Peak Patch initial conditions cosmological field described by the file
        #         field_file, or None to automatically select field type based on field_file.
        #         Allowed field_type values are
        #             rhog  : Gaussian energy density field (with mean subtracted out)
        #             rho   : energy density field (with mean subtracted out)
        #             zetag : Gaussian comoving scalar perturbation field $\zeta$
        #             zeta  : comoving scalar perturbation field $\zeta$
        #             chi   : other early universe field
        #
        #     method : str
        #         The method used to calculate the power spectrum. The current allowed methods are
        #             'total'   : to directly measure <|\psi(k)|^2> by taking the average over the
        #                         entire Fourier transform of the cosmological field's square
        #             'partial' : uses a subresolution method where we only look at every Nth k
        #                         value (currently with N fixed at 4, should update to be any input
        #                         subres value)
        #             'random'  : take a random sample of modes (not implemented yet)
        #         Default is 'total'.
        #
        #     overwrite : bool
        #         Reading power spectra is one of the most computationally expensive tasks run in
        #         peakpatchtools, which is why it's been written as a separate (hopefully soon to be
        #         fully parallelised) Fortran programme that is compiled and run by this python
        #         script. To save you from having to run and rerun this cript, the final power
        #         spectrum is saved in a common format, so that we can check if one has already been
        #         made. If overwrite is set to False, the Fortran programme is only run if no
        #         existing power spectrum file can be found. If overwrite is set to True, then any
        #         existing script is overwritten.
        # 
        #     verbose : bool
        #         If set to False, text output to screen is suppressed. Default is True.

        # To do:
        # - make script get_partial_power.f90 take a subres value instead of just using default 4
        # - make a script to make rho-rhog and zeta-zetag fields and then measure power spectra

        # Types of fields and corresponding file names
        field_types = ['rhog','rho' ,'zetag','zeta'  ,'chi']
        field_names = ['rhog','Fvec','zetag','zetang','chi']
        def fname(ftype): return '{0}/fields/{1}_{2}'.format( self.run_dir, ftype, self.run_name )
        def pname(ftype): return '{0}/tables/p_{1}_{2}'.format( self.run_dir, ftype, self.run_name )
        def field_t2n(ftype): return field_names[ field_types.index(ftype) ]

        ############################################################################################
        # Determine values for kwargs if none given and locate field files                         #
        ############################################################################################

        # Assume Gaussian overdeensity field if no field type or field file specified
        if not field_type and not field_file:
            field_type,field_name='rhog','rhog'

        # Try to read field type from file name if only field file given
        elif not field_type:
            for i in range(len(field_types)):
                if os.path.basename(field_file) == os.path.basename(fname(field_names[i])):
                    field_type = field_types[i] 
                    field_name = field_names[i]
                    break
                    
        # If field type given, determine which allowed type it refers to
        else:
            if type(field_type)==str:
                field_type = field_type.lower()
                if field_type not in field_types:
                    raise OSError(('Field type {0} not recognized. Allowed types are:\n\'{1}\', \''+
                        '{2}\', \'{3}\', \'{4}\', \'{5}\'').format( field_type,
                        field_types[0],field_types[1],field_types[2],field_types[3],field_types[4]))
                else:   
                    field_name = field_t2n(field_type)
            elif type(field_type)==int and field_type in range(5):
                field_type = field_types[field_type]
                field_name = field_t2n(field_type)
            else:
                raise OSError(('Field type {0} not recognized. Allowed types are:\n\'{1}\', \'{2}'+
                    '\', \'{3}\', \'{4}\', \'{5}\'').format( field_type,
                    field_types[0],field_types[1],field_types[2],field_types[3],field_types[4]))
                    
        # Auto-locate field file if none passed
        if not field_file:
            field_file = fname(field_name)
            if not os.path.exists(field_file):
                if field_type=='zeta' and os.path.exists(fname('zetag')): field_file=fname('zetag')
                elif os.path.exists(fname('Fvec')):                       field_file=fname('Fvec')
                else:
                    raise OSError(('file not found: {0}\nUnable to auto-locate field file.'
                                  ).format(field_file) )
        
        ############################################################################################
        # Determine if power spectrum has been read before, use it or re-calculate and overwrite   #
        ############################################################################################

        # Now we check if a power spectrum already exists for this field
        if os.path.exists( pname(field_type) ):
            if overwrite==False:
                if verbose==True:
                    print( ('It appears that a power spectrum has already been read from this fiel'+
                            'd:\n{0}\nThis spectrum will be used. To measure the power spectrum ag'+
                            'ain and overwrite the existing file, run\nthis function again with ar'+
                            'gument overwrite=True.').format(pname(field_type)) )
            else:
                os.system( 'rm -rf {0}'.format(pname(field_type)) )
        else:
            overwrite=True

        ############################################################################################
        # Calculate power spectrum                                                                 #
        ############################################################################################

        # If overwrite is set to True at this point, then we run the power sepctrum measuring script
        if overwrite == True:

            # Copy (and overwrite) the tools directory from Peak Patch to the run directory
            if os.path.exists( '{0}/tools'.format( self.run_dir ) ):
                os.system( 'rm -rf {0}/tools'.format( self.run_dir ) )
            os.system( 'cp -r {0}/tools {1}/tools'.format( peak_patch_dir, self.run_dir ) )

            # Check what machine you're using
            host = str(subprocess.check_output('hostname'))

            # Run fortran script for Niagara
            if host[:5]=='b\'nia' and host[-16:]=='.scinet.local\\n\'':
                #os.system( 'cd {0}/tools/;python get_powerspectrum_niagara.py {1} {2} {3}'
                #           .format( self.run_dir, field_file, self.nbuff, self.boxsize ) )
                # Haven't really built in Niagara support yet. My get_pwoerspectrum.f90 scripts
                # should really be parallelised and set up to run as jobs on Niagara or Sunnyvale.
                raise ValueError('Haven\'t written scripts for submitting Niagara jobs yet.')

            # Run fortran script for CITA machines
            elif host in [ 'b\'sheep\\n\''  , 'b\'mussel\\n\''  , 'b\'kingcrab\\n\'' ,
                           'b\'homard\\n\'' , 'b\'lobster\\n\'' , 'b\'calamari\\n\'' ,
                           'b\'shrimp\\n\'' , 'b\'prawn\\n\''   ]:
                
                # Fortran script and executable to read power spectrum from field file
                if method == 'total':
                    f90 = '{0}/tools/get_power_omp.f90'.format( self.run_dir )
                    exe = '{0}/tools/get_power_omp'.format( self.run_dir )
                elif method == 'partial':
                    f90 = '{0}/tools/get_partial_power.f90'.format( self.run_dir )
                    exe = '{0}/tools/get_partial_power'.format( self.run_dir )

                # Modules to load on CITA computers
                intel_module = 'intel'
                fftw_module  = 'fftw'

                # Compilers and flags for running
                fortrancompiler = 'ifort'
                fortranoptions  = ('-w'                  # suppress all warnings
                                  +' -fpp'               # use the Fortran preprocessor
                                  +' -qopenmp -parallel')# use OpenMP compiler directives
                                  #+' -traceback'       )# make debugging easier with traceback
                fftw_include    = '-I/cita/modules/fftw/3.3.10/include'
                fftw_lib_dir    = '-L/cita/modules/fftw/3.3.10/lib'
                fftw_libs       = '-lfftw3f_omp -lfftw3f_threads -mkl'# -lfftw3f_mpi -lfftw3f -lm'
                # above paths can be found by running e.g. `module show fftw` after loading modules
                # Note that if warnings are not suppressed with -w then you get annoying warnings
                # about formatting of write statements where P(k) is written to a file

                # Run scripts
                os.system(('chmod -x {7};module purge;module load {0} {1};'+
                           '{2} {3} {4} {5} {6} {7} -o {8};'+
                           '{8} << EOF\n"{9}"\n{10}\n{11}\nEOF')
                           #        {0}            {1}            {2}               {3}
                           .format( intel_module , fftw_module  , fortrancompiler , fortranoptions ,
                           #        {4}            {5}            {6}               {7}
                                    fftw_include , fftw_lib_dir , fftw_libs       , f90            ,
                           #        {8}            {9}            {10}              {11}
                                    exe          , field_file   , self.nbuff      , self.boxsize  ))

                # Cleanup
                os.system('rm -rf {0}'.format(exe))

                # Define the file containing the power spectrum so that it can be read in
                if method == 'total':
                    os.system('mv {0}_power.dat {1}'.format( field_file, pname(field_type) ))
                elif method == 'partial':
                    os.system('mv {0}_power.dat {1}_partial.dat'
                        .format( field_file, pname(field_type)[:-4] ))
                    warnings.warn('Note that pname currently won\'t point to files ending in _parti'
                            +'al.dat so the code can\'t see if\nyou\'ve calculated one of those alr'
                            +'eady and it\'s always gonna overwrite.',Warning)
                else:
                    raise ValueError('method {0} not understood'.format(method))

            # Run fortran script for Sunnyvale
            elif host[-13:] == '.sunnyvale\\n\'' :
                
                # Fortran script and executable to read power spectrum from field file
                if method == 'total':
                    f90 = '{0}/tools/get_power_omp.f90'.format( self.run_dir )
                    exe = '{0}/tools/get_power_omp'.format( self.run_dir )
                elif method == 'partial':
                    f90 = '{0}/tools/get_partial_power.f90'.format( self.run_dir )
                    exe = '{0}/tools/get_partial_power'.format( self.run_dir )

                # Modules to load on CITA computers
                openmpi_module = 'openmpi/4.1.6-gcc-ucx'
                fftw_module    = 'fftw/3.3.10-openmpi-ucx'
                gsl_module     = ''#'gsl/2.7.1'

                # Compilers and flags for running
                fortrancompiler = 'mpifort -DOMPI_SKIP_MPICXX -fallow-argument-mismatch'
                fortranoptions  = ('-Wextra'
                                  +' -mcmodel=large -fno-common -Wno-deprecated -fPIC'
                                  #+' -cpp'               # use the Fortran preprocessor
                                  +' -fopenmp'# -parallel' # use OpenMP compiler directives
                                  +' -O3 -ffast-math -march=znver3 -mtune=znver3 -funroll-loops'
                                  +' -finline-functions -march=native -flto -fno-fat-lto-objects'
                                  ) # + some starq flags
                                  #+' -traceback'       )# make debugging easier with traceback
                fftw_include    = '-I/cita/modules/fftw/3.3.10-openmpi-ucx/include'
                fftw_lib_dir    = '-L/cita/modules/fftw/3.3.10-openmpi-ucx/lib'
                fftw_libs       = '-lfftw3f_omp -lfftw3f_threads' # ilfftw3f -lstdc++ -mkl -lfftw3f_mpi -lfftw3f -lm'
                # above paths can be found by running e.g. `module show fftw` after loading modules
                # Note that if warnings are not suppressed with -w then you get annoying warnings
                # about formatting of write statements where P(k) is written to a file

                # Run scripts
                os.system(('chmod -x {8};module purge;module load {0} {1} {2};'+
                           '{3} {4} {5} {6} {7} {8} -o {9};'+
                           '{9} << EOF\n"{10}"\n{11}\n{12}\nEOF')
                           #        {0}              {1}            {2}            {3}
                           .format( openmpi_module , fftw_module  , gsl_module   , fortrancompiler ,
                           #        {4}              {5}            {6}            {7}
                                    fortranoptions , fftw_include , fftw_lib_dir , fftw_libs       ,
                           #        {8}             {9}             {10}           {11}
                                    f90            , exe          , field_file   , self.nbuff      ,
                           #        {12}
                                    self.boxsize  ))

                # Cleanup
                os.system('rm -rf {0}'.format(exe))

                # Define the file containing the power spectrum so that it can be read in
                if method == 'total':
                    os.system('mv {0}_power.dat {1}'.format( field_file, pname(field_type) ))
                elif method == 'partial':
                    os.system('mv {0}_power.dat {1}_partial.dat'
                        .format( field_file, pname(field_type)[:-4] ))
                    warnings.warn('Note that pname currently won\'t point to files ending in _parti'
                            +'al.dat so the code can\'t see if\nyou\'ve calculated one of those alr'
                            +'eady and it\'s always gonna overwrite.',Warning)
                else:
                    raise ValueError('method {0} not understood'.format(method))

            # Assume cita server and run fortran script for CITA machines
            else:
                raise ValueError('I don\'t recognize the machine you\'re trying to run on.')

        # Read in the power spectrum P(k)
        p = np.loadtxt( pname(field_type), delimiter=' ' )

        # Save power spectrum as intrinsic variable
        if   field_type == 'rhog' : self.k_from_rhog  , self.p_from_rhog  = p[:,0] , p[:,1]
        elif field_type == 'rho'  : self.k_from_rho   , self.p_from_rho   = p[:,0] , p[:,1]
        elif field_type == 'zetag': self.k_from_zetag , self.p_from_zetag = p[:,0] , p[:,1]
        elif field_type == 'zeta' : self.k_from_zeta  , self.p_from_zeta  = p[:,0] , p[:,1]
        else                      : self.k_from_chi   , self.p_from_chi   = p[:,0] , p[:,1]



    def get_k_pulse(self, overwrite=False, method='zeta', verbose=True):
        # Once you've run the get_power_spectrum() for field_type='zetag' and feild_type='zeta',
        # this script calculates the value of $k_\mathrm{pulse}$, the approximate scale of PING non-
        # Gaussianities. Note that $k_\mathrm{pulse}$ is intended to be used only with Peak Patch
        # non-Gaussianity models 6, 7 and 8. In this model, it is approximated as being the maximum
        # value of the PING power spectrum divided by the Gaussian power spectrum. The value is
        # calculated and stored in an attribute of class PeakPatch, k_pulse.
        # 
        # The standar deviation is also calculated from the power spectrum
        # $ \sigma^2 = \frac{1}{2\pi^2} \int k^2 P(k) dk $
        # 
        # Optional parameters:
        # 
        #     overwrite : bool
        #         If True, previous value of k_pulse is overwritten, if not, the calculation is not
        #         run.
        # 
        #     method : str
        #         Determines whether to calcualte k_pulse from overdensity or zeta fields. Options
        #         are 'zeta' to calculate based on the non-Gaussian zeta field or 'rho' to caluclate
        #         based on the nno-Gaussian density field. The default is 'zeta'.
        #
        #     verbose : bool
        #         If set to False, text output to screen is suppressed. Default is True.

        # Check if k_pulse has already been calculated
        if hasattr(self,'k_pulse'):
            if overwrite == False:
                if verbose == True:
                    print( justify_text( 'This run already has a specified value of k_pulse. To ov'+
                                         'erwrite that value, use overwrite=True. Exiting.' ))
                return

        # To caluclate using zeta method
        if method == 'zeta':

            # Check that PING and Gaussian zeta powers have been read from the IC fields
            if not hasattr(self,'p_from_zetag') or not hasattr(self,'p_from_zeta'):
                raise AttributeError(justify_text( ('The run {0} did not have attributes `p_from_z'+
                        'etag\' and/or `p_from_zeta\'. Make sure to run <this_run>.get_power_spect'+
                        'rum(field_type=\'zeta\') and <this_run>.get_power_spectrum(field_type=\'z'+
                        'etag\') before using this function.').format(str(self)),
                    width=100, first_line_buffer='AttributeError' ))

            # Calculate k_pulse as maximum value of P_{\zeta\zeta}(k) / P_{\zeta\zeta,Gaussian}(k)
            P = self.p_from_zeta / self.p_from_zetag
            j_max = np.argmax( P )
            self.k_pulse_from_zeta = self.k_from_zetag[ j_max ]

            # Calcualte zeta_pulse
            self.zeta_pulse = P[ j_max ]

            # Integrate power to calculate sigma_pulse
            P_interp = ( P[1:]                + P[:-1]                ) / 2
            k_interp = ( self.k_from_zeta[1:] + self.k_from_zeta[:-1] ) / 2
            dk       =   self.k_from_zeta[1:] - self.k_from_zeta[:-1]
            self.sigma_zeta_pulse = np.sqrt( np.sum( k_interp**2 * P_interp * dk ) /(2*np.pi**2) )

        # To caluclate using zeta method
        elif method == 'delta':

            # Check that PING and Gaussian zeta powers have been read from the IC fields
            if not hasattr(self,'p_from_deltag') or not hasattr(self,'p_from_delta'):
                raise AttributeError(justify_text( ('The run {0} did not have attributes `p_from_d'+
                        'eltag\' and/or `p_from_delta\'. Make sure to run <this_run>.get_power_spe'+
                        'ctrum(field_type=\'delta\') and <this_run>.get_power_spectrum(field_type='+
                        '\'deltag\') before using this function.').format(str(self)),
                    width=100, first_line_buffer='AttributeError' ) )

            # Calculate k_pulse as maximum value of P_{\delta\delta}(k)/P_{\delta\delta,Gaussian}(k)
            P = self.p_from_delta / self.p_from_deltag
            j_max = np.argmax( P )
            self.k_pulse_from_delta = self.k_from_deltag[ j_max ]

            # Calcualte zeta_pulse
            self.delta_pulse = P[ j_max ]

            # Integrate power to calculate sigma_pulse
            P_interp = ( P[1:]                 + P[:-1]                 ) / 2
            k_interp = ( self.k_from_delta[1:] + self.k_from_delta[:-1] ) / 2
            dk       =   self.k_from_delta[1:] - self.k_from_delta[:-1]
            self.sigma_delta_pulse = np.sqrt( np.sum( k_interp**2 * P_interp * dk ) /(2*np.pi**2) )

        # Calculate based on the chi power spectrum
        elif method == 'chi':

            # Check that chi power has been read from IC fields
            if not hasattr(self,'p_from_chi'):
                raise AttributeError(justify_text( ('The run {0} did not have attribute `p_from_ch'+
                        'i\'. Make sure to run <this_run>.get_power_spectrum(field_type=\'chi\') b'+
                        'efore using htis function.').format(str(self)),
                    width=100, first_line_buffer='AttributeError' ) )

            # Calculate k_pulse as the maximum value of P_{\chi\chi}(k)
            P = self.k_from_chi**2 * self.p_from_chi / (2*np.pi**2)
            j_max = np.argmax( P )
            self.k_pulse_from_chi = self.k_from_chi[ j_max ]
            
            # Calculate chi_pulse
            self.chi_pulse = P[ j_max ]

            # Calculate sigma_pulse
            P_interp = ( P[1:]               + P[:-1]               ) / 2
            k_interp = ( self.k_from_chi[1:] + self.k_from_chi[:-1] ) / 2
            dk       =   self.k_from_chi[1:] - self.k_from_chi[:-1]
            self.sigma_chi_pulse = np.sqrt( np.sum( k_interp**2 * P_interp * dk ) /(2*np.pi**2) )

        # Calcualte from the theoretical initial conditions chi power, not the measured power of chi
        elif method == 'chi_ic':
            
            # Check for p_chichi.dat

            self.add_IC_power_spectra(spectrum_type='P_ping')
            
            # Calculate k_pulse as the maximum value of P_{\chi\chi}(k)
            P = self.k_ping**2 * self.p_chichi_ping / (2*np.pi**2)
            j_max = np.argmax( P )
            self.k_pulse_from_chi_ic = self.k_ping[ j_max ]
            
            # Calculate chi_pulse
            self.chi_ic_pulse = P[ j_max ]

            # Calculate sigma_pulse
            P_interp = ( P[1:]           + P[:-1]           ) / 2
            k_interp = ( self.k_ping[1:] + self.k_ping[:-1] ) / 2
            dk       =   self.k_ping[1:] - self.k_ping[:-1]
            self.sigma_chi_ic_pulse = np.sqrt( np.sum( k_interp**2 * P_interp * dk ) /(2*np.pi**2) )



    # This function determines a measure of the PING non-Gaussianity strength by fitting the ratio
    # of the zeta powers in the non-Gaussian and Gaussian cases to a function. The function D(k) is
    # well modelled by a hyperbola in lnD vs lnk, which asymptotes to D ~ k^3 at low k and D ~ const
    # at some cutoff wavelength k_cutoff,
    #
    #     D(k) = D_0 n^3 k^3 k_puse^-3 exp[ 3 ln^2(n) / ln( k / n k_pulse ) ]
    #
    # where k_cutoff = n k_pulse, and the hyperbola ln[D(k)] reaches a maximum D(k_pulse) = D_0 at k=k_pulse.
    def fit_D0_kpulse_n( self, rg=None, fit_lims=None ):
        from scipy.optimize import curve_fit

        # If no Gaussian run is passed, set it to the same as the non-Gaussian run
        if rg is None:
            rg = self

        # Check that the PeakPatch objects have the appropriate power spectra read already
        if not ( hasattr(self,'k_from_zeta') and hasattr(self,'p_from_zeta') ):
            if not hasattr(self,'zeta'):
                self.add_field( field_type='zeta' )
            self.get_power_spectrum( field_type='zeta' )
        if not ( hasattr(rg,'k_from_zetag') and hasattr(rg,'p_from_zetag') ):
            if not hasattr(rg,'zetag'):
                rg.add_field( field_type='zetag' )
            rg.get_power_spectrum( field_type='zetag' )

        # Get D(k) from power spectra
        k  = self.k_from_zeta
        pg = k**3 / (2*np.pi**2) * rg.p_from_zetag
        p  = k**3 / (2*np.pi**2) * self.p_from_zeta
        D  = p/pg-1

        # Don't fit any part of the function that is equal to 0
        zero = np.array([0])
        _D_ = np.concatenate(( zero, D, zero ))
        i_pulse = np.argmax(_D_)
        fit_min = i_pulse + 1 - (_D_[i_pulse+1::-1]<=0).argmax()
        fit_max = i_pulse     + (_D_[i_pulse+1:   ]<=0).argmax()
        del(i_pulse)

        # Confine further if fit_lims specifies a stronger constraint on the k domain
        if not fit_lims is None:
            if fit_lims[0] > fit_min:
                fit_min = fit_lims[0]
            if fit_lims[1] < fit_max:
                fit_max = fit_lims[1]

        # Get ln D and ln k
        lnk = np.log( k[fit_min:fit_max] )
        lnD = np.log( D[fit_min:fit_max] )

        # Initial guesses for D_0 and k_pulse
        i_pulse   = np.argmax(lnD)
        lnk_pulse = lnk[i_pulse]
        lnD_0     = lnD[i_pulse]

        # Find k_cutoff, where D(k) goes to 0
        dlnD_dlnk_min = 0
        i_dlnD_dlnk   = i_pulse-1
        for j in range(i_pulse,len(lnk)):
            if D[j] <= 0:
                break
            dlnD_dlnk = np.log( D[j]/D[j-1] ) / np.log( k[j]/k[j-1] )
            if dlnD_dlnk < dlnD_dlnk_min:
                i_dlnD_dlnk   = j-1
                dlnD_dlnk_min = dlnD_dlnk

        # Initial guess for k_cutoff and therefore n
        i_cutoff   = i_dlnD_dlnk
        lnk_cutoff = lnk[i_cutoff]
        lnn        = lnk_cutoff - lnk_pulse

        # Trim any trailing zeros before fitting
        lnk = lnk[:i_cutoff]
        lnD = lnD[:i_cutoff]

        # Find fit parameters
        lnk_pulse, lnD_0, lnn = curve_fit( lnD_of_lnk_n,
                                           lnk,
                                           lnD,
                                           p0 = [ lnk_pulse, lnD_0, lnn ]
                                         )[0]

        # Save the fit as attributes of the PeakPatch class object
        self.D_ping_fit_params = {
                'D_0'     : np.exp(lnD_0),
                'k_pulse' : np.exp(lnk_pulse),
                'n'       : np.exp(lnn)
                }

        # Return fit values
        return ( self.D_ping_fit_params['k_pulse'],
                 self.D_ping_fit_params['D_0'    ],
                 self.D_ping_fit_params['n'      ]  )



    ################################################################################################
    ###   WebSky maps   ############################################################################



    ################################################################################################
    ###   WebSky power spectra   ###################################################################



    ################################################################################################
    ###                                                                                          ###
    ###                          PLOTTING Peak Patch AND WebSky RESULTS                          ###
    ###                                                                                          ###
    ################################################################################################
    # The functions in this section are used for plotting of the Peak Patch and WebSky data prodcuts
    # loaded into a PeakPatch object using the load functions above and optionally manipulated using
    # the post-processing functions above. Included
    # functions are:
    # 
    # Figure formatting functions:
    # 
    #     make_fig() : make a figure object
    # 
    #     make_ax() : make an axis object
    #
    # Halo mass functions and histograms:
    #
    #     plot_halos() : plot halos as circles
    #  
    #     plot_halo_hist() : make histograms of Dark Matter halos catalogued by Peak Patch
    #     
    #     plot_halo_hist2d() : make 2d histograms of DM halos, analagous to slices of ICs fields
    #     
    #     plot_halo_hist3d() : make 3d histograms of DM halos, analagous to ICs fields
    # 
    #     plot_hmf() : make a halo mass function of various forms
    # 
    #     plot_hmzf() : make a 2D halo function of mass and redshift
    # 
    # Halo correlation functions
    # 
    # ICs fields:
    # 
    #     plot_field_slice() : make a plot of a slice of 3D initial conditions fields
    # 
    # Power spectra:
    # 
    #     
    # 
    # WebSky maps:
    # 
    # WebSky power spectra:
    # 
    # See bolow for more information on their usage.

    ################################################################################################
    ###   Figure formatting functions   ############################################################

    def add_pcolormesh(fig,ax):
        ax.set_aspect(1)
        fig.colorbar( pcm, cax=make_axes_locatable(ax).append_axes('right', size='5%', pad=0.05) )
        return fig,ax
    
    # def make_fig( fig,axs,axobjs ):
    #     
    #     if type(axobjs) != dict and type(axobjs) != tuple:
    #         raise TypeError('axobjs must be dict or tuple of dicts.')
    #     elif ( repr(type(axs)) != "<class 'matplotlib.axes._subplots.AxesSubplot'>" and
    #            type(axs) != np.ndarray ):
    #         raise TypeError('axs must be an array or AxesSubplot.')
    #     elif type(axs) == np.ndarray:
    #         if ( repr(type(axs[0])) != "<class 'matplotlib.axes._subplots.AxesSubplot'>" and
    #                type(axs[0]) != np.ndarray ):
    #             raise TypeError('axs must be an array or AxesSubplot.')
    #         elif type(axs[0]) == np.ndarray:
    #             if repr(type(axs[0])) != "<class 'matplotlib.axes._subplots.AxesSubplot'>":
    #                 raise TypeError('axs must be an array or AxesSubplot.')
    #             elif int(axs.shape[0]*axs.shape[1]) != len(axobjs):
    #                 raise IndexError('axs and axobjs must have the same number of elements.')
    #         elif len(axs) != len(axobjs):
    #             raise IndexError('axs and axobjs must have the same number of elements.')
    #     elif type(axobjs) == dict:
    #         if repr(type(axs)) != "<class 'matplotlib.axes._subplots.AxesSubplot'>":
    #             raise IndexError('axs and axobjs must have the same number of elements.')
    #         else:
    #             axobjs = (axobjs,)
    #
    #     # For subplots
    #     elif type(axs) == tuple:
    #         for axobjs_j in axobjs:
    #             if type(axobjs_j) != dict: raise TypeError('each element of tuple must be a dict.')
    #
    #     # Loop over
    #     for i in range(len(
    #         # Returns pcolormesh and labels
    #         return axtype,pcm,xlabel,ylabel,cbarlabel
    # 
    # def make_ax(): return None



    ################################################################################################
    ###   Halo Mass Functions and Histograms   #####################################################

    def plot_halos_2D( self, ax, plane='xz' , coord='lagrangian' , selection='slice',
                       method='true size', intercept=None, thickness=None, cmap=None,
                       color='k', lw=0.5 ):
        # 
        # Arguments
        # 
        #     method
        #         'intersection'
        #         'true size'
        #
        #     selection
        #         'cross section'
        #         'slice'

        # 
        theta = np.linspace(0,2*np.pi,100)

        if thickness == None:
            thickness = 2*np.max(self.R_th)

        # To plot halos in the yz plane
        if plane == 'yz':

            # Automatically set the intercept if no intercept has been set
            if intercept == None:
                intercept = self.cenx

            # For plotting Lagrangian halos
            if coord.strip().lower()[0] == 'l':

                # Forth intersection selection criterion
                if selection == 'intersection':
                    which = np.argwhere( ( abs( intercept - self.xL ) <= self.R_th ) ).T[0]

                # For the slice selection criterion
                else: # selection == 'slice':
                    which = np.argwhere( ( abs( self.xL - intercept ) <= thickness/2 ) ).T[0]

                # Set variables for x-axis, y-axis of plot and halo size R_th
                x,y,z,R = self.yL[which] , self.zL[which] , self.xL[which] , self.R_th[which]

            # For plotting Eulerian halos
            else:

                # Forth intersection selection criterion
                if selection == 'intersection':
                    which = np.argwhere( ( abs( intercept - self.x ) <= self.R_th ) ).T[0]

                # For the slice selection criterion
                else: # selection == 'slice':
                    which = np.argwhere( ( abs( self.x - intercept ) <= thickness/2 ) ).T[0]

                # Set variables for x-axis, y-axis of plot and halo size R_th
                x,y,z,R = self.y[which]  , self.z[which]  , self.x[which] , self.R_th[which]

        # To plot halos in the yz plane
        elif plane == 'xz':

            # Automatically set the intercept if no intercept has been set
            if intercept == None:
                intercept = self.ceny

            # For plotting Lagrangian halos
            if coord.strip().lower()[0] == 'l':

                # Forth intersection selection criterion
                if selection == 'intersection':
                    which = np.argwhere( ( abs( intercept - self.yL ) <= self.R_th ) ).T[0]

                # For the slice selection criterion
                else: # if selection == 'slice':
                    which = np.argwhere( ( abs( self.yL - intercept ) <= thickness/2 ) ).T[0]

                # Set variables for x-axis, y-axis of plot and halo size R_th
                x,y,z,R = self.zL[which] , self.xL[which] , self.yL[which] , self.R_th[which]

            # For plotting Eulerian halos
            else:

                # Forth intersection selection criterion
                if selection == 'intersection':
                    which = np.argwhere( ( abs( intercept - self.y ) <= self.R_th ) ).T[0]

                # For the slice selection criterion
                else: # if selection == 'slice':
                    which = np.argwhere( ( abs( self.y - intercept ) <= thickness/2 ) ).T[0]

                # Set variables for x-axis, y-axis of plot and halo size R_th
                x,y,z,R = self.z[which]  , self.x[which]  , self.y[which] , self.R_th[which]

        # To plot halos in the yz plane
        else:

            # Automatically set the intercept if no intercept has been set
            if intercept == None:
                intercept = self.cenz

            # For plotting Lagrangian halos
            if coord.strip().lower()[0] == 'l':

                # Forth intersection selection criterion
                if selection == 'intersection':
                    which = np.argwhere( ( abs( intercept - self.zL ) <= self.R_th ) ).T[0]

                # For the slice selection criterion
                else: # if selection == 'slice':
                    which = np.argwhere( ( abs( self.zL - intercept ) <= thickness/2 ) ).T[0]

                # Set variables for x-axis, y-axis of plot and halo size R_th
                x,y,z,R = self.xL[which]  , self.yL[which]  , self.zL[which] , self.R_th[which]

            # For plotting Eulerian halos
            else:

                # Forth intersection selection criterion
                if selection == 'intersection':
                    which = np.argwhere( ( abs( intercept - self.z ) <= self.R_th ) ).T[0]

                # For the slice selection criterion
                else: # if selection == 'slice':
                    which = np.argwhere( ( abs( self.z - intercept ) <= thickness/2 ) ).T[0]

                # Set variables for x-axis, y-axis of plot and halo size R_th
                x,y,z,R = self.x[which]  , self.y[which]  , self.z[which] , self.R_th[which]
        
        # To plot cross sections of the halos
        if method == 'cross section':
            which      = np.argwhere( ( abs(intercept-z) < R ) ).T[0] 
            x,y,z,R_th = x[which], y[which], z[which], R[which]
            R          = np.sqrt( R_th**2 - abs(z-intercept)**2 )

        # If colourmap
        if cmap != None:

            if method == 'cross section':
                R_rel = ( R_th - np.min(self.R_th) ) / ( np.max(self.R_th) - np.min(self.R_th) )
            else:
                R_rel = ( R - np.min(self.R_th) ) / ( np.max(self.R_th) - np.min(self.R_th) )
            for j in range(len(x)):
                ax.plot( x[j]+R[j]*np.cos(theta) , y[j]+R[j]*np.sin(theta) , lw=lw ,
                         color=cmap(R_rel[j]) )

        # If all halos are the same colour
        else:
            for j in range(len(x)):
                ax.plot( x[j]+R[j]*np.cos(theta) , y[j]+R[j]*np.sin(theta) , lw=lw , c=color )

            

    def plot_halo_hist(self, ax ): return None
        # Function for plotting halo histograms


    
    def plot_halo_hist2d(self, fig, ax, plot_type='pcolormesh' , **kwargs ):
        # Function for plotting Peak Patch halos as 2d histograms (designed to accompany funtion
        # `halo_hist2d()` above).
        # 
        # Parameters:
        # 
        #     ax : matplotlib.axes._subplots.Axessubplot
        #         An axes subplot object (e.g. one made using `fig,ax=plt.subplots`) in which to
        #         plot the 2D histogram
        # 
        #     plot_type : str
        #         The type of plot to add to the axes ax. Options are:
        #         'pcolormesh' - pseudo-colour mesh plot
        #         'contour'    - plot contour lines to show topographical view of mass field
        #         Default is 'pcolormesh'.

        # Set the default colour map as greyscale
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'Greys'

        # For this to work you probably have to remove xlim and ylim from the dict after using them
        # # Set axis limits
        # if hasattr(kwargs,'xlim'): ax.set_xlim( kwargs['xlim'][0] , kwargs['xlim'][1] )
        # if hasattr(kwargs,'ylim'): ax.set_ylim( kwargs['ylim'][0] , kwargs['ylim'][1] )

        # Plot the halo 2D histogram as a matplotlib color mesh plot
        if plot_type == 'pcolormesh':
            X,Y = np.meshgrid( self.hist2d[1], self.hist2d[2] , indexing='ij' )

            if self.mass_weighted2d == False:
                unit_density = (   ( self.hist2d[1][1] - self.hist2d[1][0] )
                                 * ( self.hist2d[2][1] - self.hist2d[2][0] ) )
            else:
                unit_density = (   ( self.hist2d[1][1] - self.hist2d[1][0] )
                                 * ( self.hist2d[2][1] - self.hist2d[2][0] ) )

            pcm = ax.pcolormesh( X, Y, self.hist2d[0]/unit_density, **kwargs )
            ax.set_aspect(1)

            # Make colorbar
            cax  = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)
            cbar = fig.colorbar( pcm, cax=cax )
            if 'cbarlabel' in kwargs:
                cbarlabel = kwargs['cbarlabel']
            else:
                if self.mass_weighted2d == False:
                    cbarlabel=r'Halo numberdensity $[\mathrm{{Mpc}}^{{-2}}]$'
                else:
                    cbarlabel=r'Halo numberdensity $[M_\odot \mathrm{{Mpc}}^{{-2}}]$'
            if 'rotation' in kwargs: rotation = kwargs['rotation']
            else:                    rotation = 90
            cbar.set_label( cbarlabel, rotation=rotation )

        # Set axis labels
        ax.set_xlabel(r'$x ~ [\mathrm{Mpc}]$')
        ax.set_ylabel(r'$z ~ [\mathrm{Mpc}]$')

        return pcm,ax

# use self.mass_weighted2d and self.hist2d etc

# >>> import os, sys, matplotlib.pyplot as plt, numpy as np ; sys.path.insert( 0 , '<...>/peakpatch/python' ) ; from peakpatchtools import PeakPatch
# >>> r = PeakPatch('/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.04.29_LIM_production_runs/ng7_cenz7150Mpc_C_mlambda25/')
# >>> r.add_halos()
# >>> fig,ax = plt.subplots()
# >>> hist = r.halo_hist2d(1)
# >>> c,ax = plot_halo_hist2d(ax)
# >>> fig.colorbar( c, ax=ax )
# >>> fig.savefig('/cita/d/www/home/njcarlson/figures/thing.png',dpi=300,bbox_inches='tight')

#/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.04.29_LIM_production_runs/ng0_cenz7500Mpc_R/




    ################################################################################################
    ###   Halo correlation functions   #############################################################

    ################################################################################################
    ###   ICs fields   #############################################################################

    def field_label(field_type):
        field_types = ['rhog','rho','zetag','zeta','chi']
        labels      = [r'$\bar{\rho}_G\delta_G$' , r'$\bar{\rho}\delta$' , r'$\zeta_G$' , 
                       r'$\zeta$' , r'$\chi$' ] 
        for i in range(len(field_types)):
            if field_type.lower() == field_types[i] : return labels[i]

        # If the argument is not a recognized field type
        raise ValueError('Field type {0} not recognized.'.format(field_type))



    def plot_field_slice( self, fig, ax, field_type, field=None, plot_type='pcolormesh', plane='xz',
                          intercept=None, **kwargs ):
        # Function for plotting Peak Patch initial conditions fields from attributes rho, rhog,
        # zeta, zetag, and chi or from a field passed at command line. Note that these fields should
        # be formatted as PeakPatch field slice attributes are (see add_field above).
        # 
        # Parameters
        # 
        #     fig : matplotlib figure
        #         The figure in which to plot the field slice.
        # 
        #     ax : matplotlib.axes._subplots.Axessubplot
        #         An axes subplot object (e.g. one made using `fig,ax=plt.subplots`) in which to
        #         plot the 2D histogram, one of the axes contained within the figure, fig.
        #
        #     field_type : str
        #         One of ['rhog','rhong','rho','zetag','zetang','zeta','chi'] defining the type of
        #         field to plot on ax. If 'rhong' or 'zetang' and field=None, the field will be the
        #         rho-rhog or zeta-zetag respectively.
        # 
        #     field : None or [ np.ndarray, str, np.ndarray ], default None
        #         If None, slices are plotted from the field attributes of the PeakPatch object (see
        #         parameter field_type above for more on this). Otherwise it must be a list format-
        #         ted like a PeakPatch field attributes.
        # 
        #     plot_type : str, default 'pcolormesh'
        #         The type of plot to make: 'pcolormesh' for plt.pcolormesh, 'contour' for
        #         plt.contour.
        # 
        #     plane : str, default 'xz'
        #         If self has a full 3D field loaded instead of a single slice, then a plane of that
        #         full field is chosen. This parameter is only used if field[1]==None (meaning that
        #         a 3D field has been loaded). The default value is 'xz' to plot the x-z plane. The
        #         other allowed values are 'xy' and 'yz'.
        # 
        #     intercept : float or None, default None
        #         Similar to plane, this parameter is only used if a 3D field is passed. It is used
        #         to determine what plane of the 3D field to plot. The default value is None, in
        #         which case it will take the middle plane (e.g. if plane='xz' and intercept is set
        #         to None, the value of self.ceny will be used). The intercept represents a real-
        #         valued number in Mpc, if for instance the default is used and the field has an
        #         even number of cells `next', the cell that is plotted will be `next/2` otherwise
        #         it will be `(next+1)/2)'.
        # 
        #     **kwargs magic variables passed as variables and treated as a dict
        #         Allows for setting of key word arguemnets in the matplotlib pcolormesh method,
        #         such as vmin and vmax so taht you can keep colour bars consistent between plots.
        #         Supported kwargs for this function:
        #         
        #         - xlabel : str
        #               The label for the horizontal axis. If not set, use the default:
        #               "r'comoving distance ${0}~[\mathrm{{Mpc}}]$'"
        #               formatted according to the slice passed. This option is useful if you want
        #               to plot several panes with the same axis type to suppress labels and tidy up
        #               your plots.
        #         
        #         - ylabel : string
        #               The label for the vertical axis. If not set, use the default:
        #               "r'comoving distance ${0}~[\mathrm{{Mpc}}]$'"
        #               formatted according to the slice passed. This option is useful if you want
        #               to plot several panes with the same axis type to suppress labels and tidy up
        #               your plots.
        #         
        #         - xlim : tuple of floats, length 2
        #               Axis limits for the horizontal axis.
        #         
        #         - ylim : tuple of floats, length 2
        #               Axis limits for the vertical axis
        #         
        #         - vmin : float
        #               Lower limit for the colour bar.
        #         
        #         - vmax : float
        #               Upper limit for the colour bar.
        #         
        #         - cbarlabel : str
        #               The label for the colour bar. If none passed, a label is chosen based on the
        #               field passed. Like with the axis label kwargs, this option is usefull if you
        #               want to suppress automatic labels.

        # Check formatting for field type
        if type(field_type) != str: raise TypeError('field_type must be a string.')
        else                      : field_type.strip().lower()
        if field_type == 'delta' : field_type = 'rho'
        if field_type == 'deltag': field_type = 'rhog'
        if field_type not in ['rhog','rhong','rho','zetag','zetang','zeta','chi']:
            raise ValueError('Field type {0} not recognized.'.format(field_type))

        # Select field if none given as argument
        if field==None:

            # If only non-Gaussian compoenent of density field ( $\rho - \rho_G$ )
            if field_type == 'rhong':
                if not hasattr(self,'rho'):
                    raise AttributeError('{0} has no attribute \'rho\'.'.format(self))
                elif not hasattr(self,'rhog'):
                    raise AttributeError('{0} has no attribure \'rhog\'.'.format(self))
                else:
                    if ( self.rho[1] != self.rhog[1] #or self.rho[1] == None or self.rhog[1] == None
                            or (self.rho[2] != self.rhog[2]).any() ):
                        raise AttributeError('rho and rhog must be field slices along the same axi'+
                                             's of\nthe field with the same dimensions.')
                    field = [ self.rho[0]-self.rhog[0] ,
                              self.rho[1]              ,
                              self.rho[2]              ]

            # If only non-Gaussian component of zeta field ( $\zeta - \zeta_G$ )
            elif field_type == 'zetang':
                if not hasattr(self,'zeta'):
                    raise AttributeError('{0} has no attribute \'zeta\'.'.format(self))
                elif not hasattr(self,'zetag'):
                    raise AttributeError('{0} has no attribure \'zetag\'.'.format(self))
                else:
                    if ( self.zeta[1] != self.zetag[1] or self.zeta[1]==None or self.zetag[1]==None
                            or (self.zeta[2] != self.zetag[2]).any() ):
                        raise AttributeError('zeta and zetag must be field slices along the same axi'+
                                             's of\nthe field with the same dimensions.')
                    field = [ self.zeta[0]-self.zetag[0] ,
                              self.zeta[1]              ,
                              self.zeta[2]              ]

            # For other field types, save space by making a NumPy pointer
            elif field_type == 'rho'   : field = [ np.ndarray.view(self.rho  [0]), self.rho  [1],
                                                   np.ndarray.view(self.rho  [2])                 ]
            elif field_type == 'rhog'  : field = [ np.ndarray.view(self.rhog [0]), self.rhog [1],
                                                   np.ndarray.view(self.rhog [2])                 ]
            elif field_type == 'zeta'  : field = [ np.ndarray.view(self.zeta [0]), self.zeta [1],
                                                   np.ndarray.view(self.zeta [2])                 ]
            elif field_type == 'zetag' : field = [ np.ndarray.view(self.zetag[0]), self.zetag[1],
                                                   np.ndarray.view(self.zetag[2])                 ]
            elif field_type == 'chi'   : field = [ np.ndarray.view(self.chi  [0]), self.chi  [1],
                                                   np.ndarray.view(self.chi  [2])                 ]

        # Check that field passed as argument is formatted correctly
        else:
            if not type(field) == list: raise TypeError('field must have type list.')
            elif len(field) != 3: raise IndexError('field must be list of length 3.')
            elif type(field[0])!=np.ndarray or type(field[1])!=str or type(field[2])!=np.ndarray:
                raise TypeError('field indices must have type [ np.ndarray, str, np.ndarray ]')
            elif len(field[0].shape) != 2 and len(field[2].shape) != 2 and len(field[2]) != 2:
                raise IndexError('field[0].shape should be (n,n) where field[2].shape is (2,n).')
            elif field[0].shape != ( field[2].shape[1] , field[2].shape[1] ):
                raise IndexError('field[0].shape should be (n,n) where field[2].shape is (2,n).')
            elif field[1] not in ['yz','xz','xy']:
                if field[1].strip().lower() in ['yz','xz','xy']: field[1] = field[1].strip().lower()
                else: raise ValueError('field[1] must be in [\'yz\',\'xz\',\'xy\'].')

        # If a `field' is 3D not a plane
        if field[1] == None:
            
            # Make sure plane is properly formatted
            if type(plane) != str:
                raise TypeError('plane must be a string.')
            field[1] = plane.strip().lower()
            if field[1] not in ('yz','xz','xy'):
                raise ValueError('plane must be one of "yz", "xz", and "xy".')
            
            # Automatically set the intercept if it is set to None
            if intercept == None:
                if   field[1] == 'yz': intercept = self.cenx
                elif field[1] == 'xz': intercept = self.ceny
                else                 : intercept = self.cenz

            # Make sure intercept is formatted properly if not set to None
            if type(intercept) not in (float,int):
                raise TypeError('intercept must be a number or None.')

            # Setup for a yz plane
            if field[1] == 'yz':

                # If the intercept is not within the box, throw an error
                if intercept < field[2][0,0] or intercept > field[2][0,-1]:
                    raise ValueError(('the value of intercept xL={0} Mpc is not within the sim'+
                                          'ulation volume.').format(intercept))

                # Find the intercept as an `index' in the 3D field `field[0][index,:,:]'. If
                # intercept is equal to one of the bin edges, index will be the higher bin
                # except for the final bin. This is ensured by setting index initially to the
                # final bin index
                index = len( field[2][0] ) - 1
                for i in range( index ):
                    if intercept >= field[2][0,i] and intercept < field[2][0,i+1]:
                        index = i
                        break

                # set field[0] to just the relevant slice of the field
                field[0] = field[0][index,:,:].T
                # By convention, all field slices are viewed as if from "above", in other words,
                # the normal vector to the plane being plotted (in this case the yz plane)
                # oriented along the Cartesian basis vector (in this case x) points to the
                # observer. To achieve this, the transpose is taken above.  

            # Setup for an xz plane
            elif field[1] == 'xz':

                # If the intercept is not within the box, throw an error
                if intercept < field[2][1,0] or intercept > field[2][1,-1]:
                    raise ValueError(('the value of intercept yL={0} Mpc is not within the sim'+
                                      'ulation volume.').format(intercept))

                # Find the intercept as an `index' in the 3D field `field[0][index,:,:]'. If
                # intercept is equal to one of the bin edges, index will be the higher bin
                # except for the final bin. This is ensured by setting index initially to the
                # final bin index
                index = len( field[2][1] ) - 1
                for i in range( index ):
                    if intercept >= field[2][1,i] and intercept < field[2][1,i+1]:
                        index = i
                        break

                # set field[0] to just the relevant slice of the field
                field[0] = field[0][:,index,:]
                # By convention, all field slices are viewed as if from "above", in other words,
                # the normal vector to the plane being plotted (in this case the xz plane)
                # oriented along the Cartesian basis vector (in this case y) points to the
                # observer. To achieve this, the transpose is NOT taken above.  

            # Setop for an xy plane
            else:

                # If the intercept is not within the box, throw an error
                if intercept < field[2][2,0] or intercept > field[2][2,-1]:
                    raise ValueError(('the value of intercept zL={0} Mpc is not within the sim'+
                                      'ulation volume.').format(intercept))

                # Find the intercept as an `index' in the 3D field `field[0][index,:,:]'. If
                # intercept is equal to one of the bin edges, index will be the higher bin
                # except for the final bin. This is ensured by setting index initially to the
                # final bin index
                index = len( field[2][2] ) - 1
                for i in range( index ):
                    if intercept >= field[2][2,i] and intercept < field[2][2,i+1]:
                        index = i
                        break

                # set field[0] to just the relevant slice of the field
                field[0] = field[0][:,:,index].T
                # By convention, all field slices are viewed as if from "above", in other words,
                # the normal vector to the plane being plotted (in this case the xy plane)
                # oriented along the Cartesian basis vector (in this case z) points to the
                # observer. To achieve this, the transpose is taken above.  
                    
        # Used later to set default axis labels
        l = r'comoving distance ${0}~[\mathrm{{Mpc}}]$'

        # Set default x and y axis labels
        if field[1] == 'yz':
            xlabel , ylabel = l.format('y') , l.format('z')
            axes = ( field[2][1] , field[2][2] )
        elif field[1] == 'xz': 
            xlabel , ylabel = l.format('z') , l.format('x')
            axes = ( field[2][2] , field[2][0] )
        elif field[1] == 'xy':
            xlabel , ylabel = l.format('x') , l.format('y')
            axes = ( field[2][0] , field[2][1] )

        # Overwrite default axis labels if something else is passed in **kwargs
        if 'xlabel' in kwargs: xlabel = kwargs['xlabel']
        if 'ylabel' in kwargs: ylabel = kwargs['ylabel']

        # Set axis limits
        if 'xlim' in kwargs: ax.set_xlim( kwargs['xlim'][0] , kwargs['xlim'][1] )
        if 'ylim' in kwargs: ax.set_ylim( kwargs['ylim'][0] , kwargs['ylim'][1] )

        # Set the colour map
        if 'cmap' in kwargs:

            # Use Planck colour map
            if kwargs['cmap'] == 'planck': cmap = planck_cmap()

            # Use Earth Tones colour map
            elif kwargs['cmap'] == 'earth_tones': cmap = earth_tones_cmap()

            # Use matplotlib default colour maps
            else: cmap = kwargs['cmap']

        # Use the default colour map (viridis)
        else:
            cmap = 'viridis'

        # Make colour mesh object
        if plot_type == 'pcolormesh':
            if 'vmin' in kwargs and 'vmax' in kwargs:
                pcm = ax.pcolormesh( axes[0] , axes[1] , field[0] , cmap=cmap,
                                     vmin=kwargs['vmin'], vmax=kwargs['vmax'])
            else:
                pcm = ax.pcolormesh( axes[0] , axes[1] , field[0] , cmap=cmap)
            ax.set_aspect(1)
            cax  = make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05)
            cbar = fig.colorbar( pcm, cax=cax )

            # Colour bar labels
            if 'cbarlabel' in kwargs:
                cbarlabel = kwargs['cbarlabel']
            else:
                if   field_type == 'rho'   : cbarlabel=r'$\bar{\rho}\delta$'
                elif field_type == 'rhog'  : cbarlabel=r'$\bar{\rho}\delta_G$'
                elif field_type == 'rhong' : cbarlabel=r'$\bar{\rho}\delta-\bar{\rho}\delta_G$'
                elif field_type == 'zeta'  : cbarlabel=r'$\zeta$'
                elif field_type == 'zetag' : cbarlabel=r'$\zeta_G$'
                elif field_type == 'zetang': cbarlabel=r'$\zeta-\zeta_G$'
                elif field_type == 'chi'   : cbarlabel=r'$\chi$'
            if 'rotation' in kwargs: rotation = kwargs['rotation']
            else:                    rotation = 270
            cbar.set_label( cbarlabel, rotation=rotation )

            # Set axis labels
            if xlabel != '': ax.set_xlabel(xlabel)
            if ylabel != '': ax.set_ylabel(ylabel)

            # Define the type of axes object used by make_fig() to make figure files
            axtype = 'Peak Patch IC field slice'

            # Put colour bar axis in scientific notation with power mulitplying the scale above
            if 'cbar_sci_not' in kwargs: cbar_sci_not = kwargs['cbar_sci_not']
            else                       : cbar_sci_not = True
            if cbar_sci_not == True:
                from matplotlib.ticker import ScalarFormatter
                formatter = ScalarFormatter(useMathText=True)
                formatter.set_powerlimits((0,0))
                cbar.ax.yaxis.set_major_formatter(formatter)

            # Returns pcolormesh and labels
            return { 'fig'           : fig       ,
                     'ax'            : ax        ,
                     'axtype'        : axtype    ,
                     'pcolormesh'    : pcm       ,
                     'cax'           : cax       ,
                     'cbar'          : cbar      ,
                     'xlabel'        : xlabel    ,
                     'ylabel'        : ylabel    ,
                     'colorbarlabel' : cbarlabel }
        
        # Make contours
        elif plot_type == 'contour':
            axtype = 'Peak Patch IC field slice contour'
            topo = ax.contour( field[2][0] , field[2][1] , field[0] ) # , levels = ...

            return { 'ax'     : ax     ,
                     'axtype' : axtype ,
                     'contour': topo   }



    ################################################################################################
    ###   Power spectra   ##########################################################################

    ################################################################################################
    ###   WebSky maps   ############################################################################

    ################################################################################################
    ###   WebSky power spectra   ###################################################################


    def savefig( self, fig, filename=None, filetype=None, bbox_inches='tight', pad_inches=0,
            transparent=False, dpi=300 ):
        #
        # Typically, figures are saved to backed-up storage in `/mnt/raid-cita/<username>/`. For ease of use
        #

        if not filename:
            filename = '/cita/d/www/home/{0}/figures/{1}'.format( os.getlogin() , self.short_name )

            # Save to raid-cita figures
            raid_cita = '/mnt/raid-cita/{0}/figures/{1}'.format( os.getlogin() , 'CASCA23/'+self.short_name )

            # Save link to website
            www = '/cita/d/www/home/{0}/figures/{1}'.format( os.getlogin() , 'CASCA23/'+self.short_name )

        # Maybe check in the run directory to make stuff better so that it's not just dumping all the figures in raid-cita/njcarlson/figures/ so that it's easier to see. Maybe also add a thing that does a mkdir -p if the directory doesn't exist.

        # Save figure to raid-cita
        fig.savefig( filename, bbox_inches=bbox_inches, pad_inches=pad_inches,
                transparent=transparent, dpi=dpi )

        # Make a link to your raid-cita figure on your website because it's a pain to look at them otherwise
        os.system( 'cd {0}; ln -s {1} {2}'.format( raid_cita, www )  )



    # other functions

    def plot_stitched_lightcone( peakpatches ):
        # 
        # assumes that you've already made a series of 
        # 
        # PARAMETERS:
        # 
        #     peakpatches : tuple of PeakPatch objects
        #         Stitches together a tuple of PeakPatch objects into a single light cone.
        # 
        # USAGE:
        # 
        # Make halo catalogues:
        # 
        #     >>> c1 = PeakPatch('<run1>')
        #     >>> c2 = PeakPatch('<run2>')
        #     ...
        #     >>> lim = ['e',np.array( [[-np.inf,np.inf],[-25,25],[-np.inf,np.inf]] ) ]
        #     >>> c1.add_halos( mass_cut='h', lim=lim )
        #     >>> c2.add_halos( mass_cut='h', lim=lim )
        #     ...
        #     >>> fig,ax = plot_stitched_lightcone( (c1,c2,...) )

        # Check that peakpatches has length greater than 1
        if not len(peakpatches) > 1: raise IndexError('this script is for stitching together the h'+
                'alo catalogues from multiple runs. If you\nwant to plot the halo catalogue from a'+
                'single run, use the appropriate PeakPatch methods to do so.')

        # In the current version, stitched lightcone runs must have the same lattice spacing and
        # share boundaries

        # Check that lattice spacings are consistent
        if ( [peakpatches[0].cellsize]*len(peakpatches) 
                != [peakpatches[j].cellsize for j in range(len(peakpatches))] or 
                [peakpatches[0].boxsize]*len(peakpatches)
                != [peakpatches[j].boxsize for j in range(len(peakpatches))] ):
            raise ValueError('all runs must have the same spatial resolution and volume.')

        # Define shared params
        boxsize = peakpatches[0].boxsize

        # Check that edges of adjacent boxes match up (Note that this assumes peakpatches are listed
        # in increasing order of cenz)
        cenzs = [ peakpatches[j].cenz for j in range(len(peakpatches)) ]
        for j in range(1,len(cenzs)):
            if cenzs[j]-boxsize/2 != cenzs[j-1]+boxsize/2:
                raise ValueError('runs must share boundaries.')

        # Define lightcone boudnaries
        lc_lims = np.array([ [ peakpatches[0].cenx-boxsize/2 , peakpatches[0].cenx+boxsize/2 ],
                             [ peakpatches[0].ceny-boxsize/2 , peakpatches[0].ceny+boxsize/2 ],
                             [ cenzs[0]-boxsize/2            , cenzs[-1]+boxsize/2           ] ])

        # Loop over all runs
        for pkp in peakpatches:

            # Read halos for each run
            pkp.add_halos()# mass_cut='high',
            #    lim=['e', np.array([[-np.inf,np.inf],[-np.inf,np.inf],[-np.inf,np.inf]],dtype=np.float32)] )

            # Make histogram for each run
            pkp.halo_hist2d( axis=1, bins=None, coord='e', bin_lim=None, mass_weighted=False )

        # Make single histogram
        X,Y = np.meshgrid( peakpatches[0].hist2d[1] , 
                np.concatenate( [ peakpatches[j].hist2d[2][:-1] for j in range(len(cenzs)) ]+
                                [ [ peakpatches[-1].hist2d[2][-1] ] ]) )#, indexing='ij' )
        hist2d = np.concatenate( [ peakpatches[j].hist2d[0] for j in range(len(cenzs)) ], axis=1 )


        fig,ax = plt.subplots()
        pcm = ax.pcolormesh( X, Y, hist2d.T )#, cmap='Greys' ) # vmin = ? , vmax = ?
        fig.colorbar( pcm, cax=make_axes_locatable(ax).append_axes('right',size='5%',pad=0.05) )
        # Label axes
        ax.set_xlabel(r'$x$ [Mpc]')
        ax.set_ylabel(r'$z$ [Mpc]')
        ax.set_aspect(1)

        # Add redshift axis on second axis
        # ax2 = ax.twiny()
        # ax2.set_xticks( ax.get_xticks() )
        # ax2.set_xbound( ax.get_xbound() )
        # ax2.set_xticklabels( [0]+[ round(z_at_value(Planck18.comoving_distance, chi*u.Mpc ).value,2) 
        #     for chi in ax.get_xticks()[1:] ] )
        # ax2.set_xlabel(r'redshift $z$')
    
        # Save to raid-cita figures
        fig_name = 'lightcone.png'
        raid_cita = '/mnt/raid-cita/{0}/figures'.format( os.getlogin() )
        www = '/cita/d/www/home/{0}/figures/CASCA'.format( os.getlogin() )
        fig.savefig( raid_cita+'/'+fig_name, dpi=600, bbox_inches='tight' )

        # Check if link exists on homespace and write/overwrite a link to the figure on raid-cita
        if os.path.exists( '{0}/{1}'.format(www,fig_name) ):
            os.system('rm -f {0}/{1}'.format(www,fig_name))
        os.system( 'cd {0};ln -s {1}/{2} ./{2}'.format( www,raid_cita,fig_name ) )

        return fig,ax,pcm



    ################################################################################################
    ###   Scripts for adding information from logfiles to run   ####################################
    ################################################################################################

    def add_logfile(self, logfile=None ):

        # Automatically set logfile
        if logfile == None:
            self.logfile = '{0}/logfiles/{1}_{2}.stdout'.format( self.run_dir , self.run_name ,
                                                                 self.seed )

        # Check that logfile exists, if it can't be found, try to locate one
        if not os.path.isfile(self.logfile):
            warnings.warn('Could not locate log file:\n{0}\nWill attempt to locate a logfile.'
                          .format(self.logfile),Warning)
            if os.path.isdir('{0}/logfiles'.format(self.run_dir)):
                counter=0
                for filename in os.listdir('{0}/logfiles'.format(self.run_dir)):
                    if filename.endswith('.stdout'):
                        counter+=1
                        self.logfile='{0}/logfiles/{1}'.format(self.run_dir,filename)
                if counter == 0:
                    raise ValueError('no file with extension \'.stdout\' was found.')
                elif counter > 1:
                    raise ValueError('more than one file with extension \'.stdout\' was found.\nYo'+
                                     'u must specify which to use.\'')
            else:
                raise ValueError('logfiles directory does not exist.')

        # Read logfiles
        flag_bad_termination = False
        with open(self.logfile) as f:
            lines = f.readlines()
            for j in range(len(lines)):
                l = lines[j].strip().replace(' ','')

                # Add memory usage
                if l.split(':')[0] == 'FFTWmemory':
                    self.FFTW_memory = float( l.split(':')[1].split('GB')[0] )
                elif l.split(':')[0] == 'Gridmemory':
                    self.Grid_memory = float( l.split(':')[1].split('GB')[0] )
                elif l.split(':')[0] == 'S2Cmemory':
                    self.S2C_memory = float( l.split(':')[1].split('GB')[0] )
                elif l.split(':')[0] == 'Catalogmemory':
                    self.Catalog_memory = float( l.split(':')[1].split('GB')[0] )
                elif l.split(':')[0] == 'Latticehuntmemory':
                    self.Lattice_hunt_memory = float( l.split(':')[1].split('GB')[0] )
                elif l.split(':')[0] == 'Totalmemory':
                    self.Total_memory = float( l.split(':')[1].split('GB')[0] )

                # Peak Finding
                elif l.split(':')[0] == 'Total number of peaks kept'.replace(' ',''):
                    self.N_peaks_kept = float( l.split(':')[1] )
                elif l.split(':')[0] == 'Total number of peaks dumped'.replace(' ',''):
                    self.N_peaks_dumped = float( l.split(':')[1] )

                # Errors
                elif l=='=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES'.replace(' ',''):
                    flag_bad_termination = True

                # Add field means and standard deviations
                elif j+2 < len(lines):
                    lp1 = lines[j+1].strip().replace(' ','').split('=')
                    lp2 = lines[j+2].strip().replace(' ','').split('=')
                    if l=='whitenoise:' and lp1[0]=='mean' and lp2[0]=='sigma':
                        if not hasattr(self,'white_noise_gaussian_mean'):
                            self.white_noise_gaussian_mean  = float(lp1[1])
                            self.white_noise_gaussian_sigma = float(lp2[1])
                        else:
                            self.white_noise_mean  = float(lp1[1])
                            self.white_noise_sigma = float(lp2[1])
                    elif l=='zeta:' and lp1[0]=='mean' and lp2[0]=='sigma':
                        if not hasattr(self,'zeta_gaussian_mean'):
                            self.zeta_gaussian_mean  = float(lp1[1])
                            self.zeta_gaussian_sigma = float(lp2[1])
                        else:
                            self.zeta_mean  = float(lp1[1])
                            self.zeta_sigma = float(lp2[1])
                    elif l=='density:' and lp1[0]=='mean' and lp2[0]=='sigma':
                        if not hasattr(self,'density_gaussian_mean'):
                            self.density_gaussian_mean  = float(lp1[1])
                            self.density_gaussian_sigma = float(lp2[1])
                        else:
                            self.density_mean  = float(lp1[1])
                            self.density_sigma = float(lp2[1])
                    elif l=='chi:' and lp1[0]=='mean' and lp2[0]=='sigma':
                        self.chi_mean  = float(lp1[1])
                        self.chi_sigma = float(lp2[1])
                    elif '-displacement:' in l and lp1[0]=='mean' and lp2[0]=='sigma':
                        if '2lpt' in l:
                            if not hasattr(self,'displacement_2LPT_mean'):
                                displacement_2LPT_mean  = [0.,0.,0.]
                                displacement_2LPT_sigma = [0.,0.,0.]
                            if 'x' in l:
                                displacement_2LPT_mean[0]  += float(lp1[1])
                                displacement_2LPT_sigma[0] += float(lp2[1])
                            elif 'y' in l:
                                displacement_2LPT_mean[1]  += float(lp1[1])
                                displacement_2LPT_sigma[1] += float(lp2[1])
                            else:
                                displacement_2LPT_mean[2]  += float(lp1[1])
                                displacement_2LPT_sigma[2] += float(lp2[1])
                        else:
                            if not hasattr(self,'displacement_1LPT_mean'):
                                displacement_1LPT_mean  = [0.,0.,0.]
                                displacement_1LPT_sigma = [0.,0.,0.]
                            if 'x' in l:
                                displacement_1LPT_mean[0]  += float(lp1[1])
                                displacement_1LPT_sigma[0] += float(lp2[1])
                            elif 'y' in l:
                                displacement_1LPT_mean[1]  += float(lp1[1])
                                displacement_1LPT_sigma[1] += float(lp2[1])
                            else:
                                displacement_1LPT_mean[2]  += float(lp1[1])
                                displacement_1LPT_sigma[2] += float(lp2[1])

        # Get total number
        if hasattr(self,'N_peaks_kept') and hasattr(self,'N_peaks_dumped'):
            self.N_peaks_premerge = self.N_peaks_kept + self.N_peaks_dumped
            if self.N_peaks_premerge == 0:
                warnings.warn( justify_text( 'The logfiles indicate that this run does not contain'+
                        'any peaks. Typically this means that the initial conditions fields includ'+
                        'e unphysical features, like if the power spectrum gets blown out or if th'+
                        'ere are NaNs in the fields.', first_line_buffer='Warning' ), Warning )
        else:
            warnings.warn( justify_text( 'Could not find the number of peaks kept and dumped in lo'+
                    'gfiles. This typically means that the run has failed at some point before hom'+
                    'ogeneous ellipsoidal collapse.', first_line_buffer='Warning' ), Warning )

        # Flag
        if flag_bad_termination == True:
            warnings.warn( justify_text( 'The logfiles indicate that this run encountered errors. '+
                    'You can look at where it indicates a "BAD TERMINATION OF ONE OF YOUR APPLICAT'+
                    'ION PROCESSES" in the logfiles or run with traceback to see where it failed.',
                    first_line_buffer='Warning' ), Warning )



    ################################################################################################
    ###   Homogeneous Ellipsoidal Collapse Table   #################################################
    ################################################################################################

    def add_homogeneous_ellipsoidal_collapse_table(self):
        # Read the interpollation table for Homogeneous Ellipsoidal collapse calculations for this
        # run.
        
        homeltab_file = self.run_dir + '/' + self.TabInterpFile

        # Read in the raw data
        homeltab_file_in = open( homeltab_file, 'rb' )
        
        # The first 9 values are a header
        self.homeltab_Nf     =   np.fromfile( homeltab_file_in, dtype=np.int32  , count=1 )[0]
        self.homeltab_Nev    =   np.fromfile( homeltab_file_in, dtype=np.int32  , count=1 )[0]
        self.homeltab_Npv    =   np.fromfile( homeltab_file_in, dtype=np.int32  , count=1 )[0]
        self.homeltab_f_min  =   10**np.fromfile( homeltab_file_in, dtype=np.float32, count=1 )[0]
        self.homeltab_f_max  =   10**np.fromfile( homeltab_file_in, dtype=np.float32, count=1 )[0]
        self.homeltab_ev_min =   np.fromfile( homeltab_file_in, dtype=np.float32, count=1 )[0]
        self.homeltab_ev_max =   np.fromfile( homeltab_file_in, dtype=np.float32, count=1 )[0]
        self.homeltab_pv_min = ( self.homeltab_ev_max * 
                                 np.fromfile( homeltab_file_in, dtype=np.float32, count=1 )[0] )
        self.homeltab_pv_max = ( self.homeltab_ev_max *
                                 np.fromfile( homeltab_file_in, dtype=np.float32, count=1 )[0] )

        # Read in all tabulated values
        self.homeltab = np.reshape(
                np.fromfile( homeltab_file_in, dtype=np.float32, count=-1 ),
                ( self.homeltab_Npv , self.homeltab_Nev , self.homeltab_Nf ) )
        # This still needs to be checked to make sure that the shape is correct



    def add_homeltab(self):
        # Convenience function for add_homogeneous_ellipsoidal_collapse_table().
        self.add_homogeneous_ellipsoidal_collapse_table()



####################################################################################################
###   Scripts for making spreadsheets of runs   ####################################################
####################################################################################################

def make_run_table( tablefile, runs, delimiter=',', file_extension='.csv', params=None ):
    # Pass a file path `tablefile' and a tuple of PeakPatch objects and the run characteristics are
    # tabulated in a `separator'-separated file of type file_extension.
    # 
    # Arguments are
    #     tablefile : str or None
    #         The path at which to save the table of run parameters.
    # 
    #     runs : list of PeakPatch
    #         A list of objects of class `PeakPatch', the python class for Peak Patch/WebSky runs. 
    #         The parameters from these runs are tabulated in the csv file.
    #
    #     delimiter : str (default is ',')
    #         The delimiter for the table of parameters. Default is ',' for a typical csv file.
    #
    #     

    default_params = [ 'seed', 'run_name', 'run_dir', 'boxsize', 'nmesh', 'nbuff', 'ntile', 'next',
        'neff', 'tlimit', 'nnodes', 'tpnode', 'ntasks', 'ncpus', 'nopmth', 'ievol', 'cenx', 'ceny',
        'cenz' ] 

    cosmology_params = [ 'Omx', 'OmB', 'Omvac', 'h', 'ns', 'As', 'sigma8', 'tau', 'mnu' ]

    map_params = [ 'maps', 'nside_map', 'npix_map', 'fov_map', 'zmin_map', 'zmax_map', 'ellmax' ]

    fNL_params  = [ 'NonGauss', 'fNL' ]

    PING_params = [ 'NonGauss', 'm_phi', 'm_chi', 'phi_w', 'phi_p', 'vev', 'm_tach', 'a_e' ]

    params = []
    if params == None:
        prams = default_params + cosmology_params
    elif 'fNL' in params:
        params += default_params + cosmology_params + fNL_params
    elif 'PING' in params:
        params += default_params + cosmology_params + PING_params
        
    if 'map' in params:
        params += map_params

    import csv
    with open( tablefile,'w',newline='' ) as csv_file:
        csvwriter = csv.writer( csv_file, delimiter=delimiter )

        # Add header row
        row = []
        for key in kwargs:
            row += [key]

        # Cycle through runs in the tuple `runs'
        for run in runs:
            
            row = []
            for key in kwargs:
                if hasattr(run,key):
                    row += [ vars(run)[key] ] 
            
    # As well as including columns for PeakPatch attributes like boxsize, nmesh, various instability
    # parameters, etc., it would also be good to include some boolean values that determine if halos
    # were found, if the homel calculation succeeded, etc. and maybe some assessment of why the run
    # might not have succeeded or where it failed if that is the case.
    return 0



def add_to_run_table( tablefile, runs ):
    # Adds a run to a preexisting table file made with `make_run_table()'.
    return 0



def sort_run_table( tablefile, sortby ):
    # Sorts a table file created using `make_run_table()'.
    return 0



####################################################################################################
###   Scripts for running components of the Peak Patch calculation   ###############################
####################################################################################################

# Script for checking if a run will work
def setup_run( n, boxsize, r_buff=36.0, cluster='niagara', n_cpu_per_node=None, iLPT=2, cat_len=11):
    # This script will do a quick assessment of possible run geometries on Niagara given a desired
    #  box sidelength `boxsize' in Mpc and resolution `n' (the sidelength of the box in discrete
    #  cells).
    # 
    # Note that as of 2024, Niagara is typically RAM limitted for Peak Patch runs, meaning that
    # there is not enough RAM/node to use every CPU without running out of RAM. Typically this means
    # that we will request about 32 CPUs per node instead of the maximum 40. This is not strictly
    # speaking allowed, as SciNet encourages jobs that leave no CPUs idle (to encourage people to
    # parallelise their jobs, so that you can run e.g. one parallel job on one node instead of 40
    # serial jobs on 40 nodes), but there's nothing we can do about it, the 202 GB (188GiB) per node
    # is simply not enough for all the Fourier space stuff we do in Peak Patch.
    #
    # Note that the Niagara supercomputer has the following characteristics:
    # Number of nodes: 2024 (80960 cores total)
    # RAM/node:        188 GiB (202 GB)                               <--- 1 GiB = 2^10 MiB = 2^30 B
    # CPU cores/node:  40 (80 hyperthreads)
    # Our allocation allows for no more than 1000 nodes per job, and all jobs are limited at 24
    # hours of runtime.
    # 
    # Parameters:
    #
    #     n : int
    #         Run resolution, defined as the number of discrete cells per side of the cubic run, so
    #         the total number of cells is N=n^3. Note that this INCLIDES BUFFERS.
    # 
    #     boxsize : float
    #         The sidelength of the box (EXCLUDING BUFFERS) in Mpc, so the total volume of the box
    #         excluding buffers is boxsize^3 Mpc^3 (or with buffers (boxsize+2*r_buff)^3 Mpc^3).
    #         Note that this is not in Mpc/h, so we assume a value of h, although there isn't really
    #         any cosmology calculation done here, the size is just relative to the buffer thickness
    #         `r_buff' which must be greater than about 32 Mpc to avoid serious edge effects from
    #         large halos.
    # 
    #     r_buff : float
    #         The buffer thickness in Mpc. See schematic below for full summary of geometries.
    #         default : 36.0 (value shown to avoid edge effects).
    # 
    #     cluster : str
    #         Which machine you're running on. 
    #         default : 'niagara' for the SciNet Niagara Supercomputer
    # 
    #     n_cpu_per_node : int
    #         The number of CPUs per node requested for your job. Niagara can use up to 40, but that
    #         often doesn't leave enough RAM per process and can result in failed runs. It is
    #         recommended to use 32 or fewer for large runs. Making this a power of 2 also
    #         simplifies things as the number of FFTW slabs must be evenly divided across processors.
    #         default : None (will automatically determine a value based on what cluster you're
    #                   using).
    # 
    #     iLPT : int
    #         Whether the run will use 1st or 2nd order Lagrangian perturbation theory (this
    #         determines the number of early universe fields we need to generate, affecting total
    #         RAM requirements for the run).
    #         default : 2
    # 
    #     cat_len : int
    #         The number of parameters used to describe each halo in the halo catalogue. Each
    #         parameter is a single 32-bit floating point number, thus the total will have an effect
    #         on the RAM requirements to save the catalogues in temporary memory during peak finding
    #         and merging.
    #         default : 11, to output full strain information and get antisymmetric halos requires
    #                   cat_len=33
    # 
    # For a run with n_tile = 2, n_buff=1, n=10, boxsize=8.0 Mpc and therefore r_buff=1 Mpc, here is
    # an axample of what a slice of the box would look like.
    # 
    # |<------------- n (voxels) ------------>|
    #     |<------- boxsize (Mpc) ------->|
    #  --- --- --- --- --- --- --- --- --- ---   -
    # | b | b | b | b | b | b | b | b | b | b |  ^
    #  --- --- --- --- --- --- --- --- --- ---   |  ---------
    # | b | 1 | 1 | 1 | 1 |2b1|   |   |   | b |  |          ^
    #  --- --- --- --- --- --- --- --- --- ---   |          |
    # | b | 1 | 1 | 1 | 1 |2b1|   |   |   | b |  |          |
    #  --- --- --- --- --- --- --- --- --- ---   |          |
    # | b | 1 | 1 | 1 | 1 |2b1|   |   |   | b |  |          |
    #  --- --- --- --- --- --- --- --- --- ---   |          |
    # | b | 1 | 1 | 1 | 1 |2b1|   |   |   | b | n (voxels)  |
    #  --- --- --- --- --- --- --- --- --- ---   |          |
    # | b |3b1|3b1|3b1|3b1|4b1|   |   |   | b |  |       boxsize (Mpc)
    #  --- --- --- --- --- --- --- --- --- ---   |          |
    # | b |   |   |   |   |   |   |   |   | b |  |          |
    #  --- --- --- --- --- --- --- --- --- ---   |          |
    # | b |   |   |   |   |   |   |   |   | b |  |          |
    #  --- --- --- --- --- --- --- --- --- ---   |          |
    # | b |   |   |   |   |   |   |   |   | b |  |          v
    #  --- --- --- --- --- --- --- --- --- ---   |  ---------
    # | b | b | b | b | b | b | b | b | b | b |  v 
    #  --- --- --- --- --- --- --- --- --- ---   - 
    # |<-- n_mesh (voxels) -->|         ->|   |<- n_buff (voxels)
    #                                   ->|   |<- r_buff (Mpc)
    # Legend:
    # b   - cells in the buffer
    # 1   - cells in tile 1
    # ib1 - cells in tile i that are also in the buffer for tile 1, note that some cells labelled 1
    # are also in the buffer for the adjaceent tiles, but I've just labelled them 1 for simplicity
    # here.
    
    ################################################################################################
    # Automatically set machine characteristics                                                    #
    ################################################################################################

    # For UofT/SciNet's Niagara Supercomputer
    if cluster.lower() == 'niagara':
        ram_per_node = 188.0 # GiB, the amount of RAM per node on Niagara
        N_nodes_max  = 1000  # number of available nodes per job on Niagara
        
        # Set the number of CPUs to use per node if set to None
        if n_cpu_per_node == None:
            n_cpu_per_node = 32 # Intel Skylake cores
                                # Niagara allows for up to 40, for large runs fewer may be necessary

    # Different queues on CITA's Sunnyvale supercomputer
    elif cluster.lower() == 'hpq':
        ram_per_node = 29.802 # GiB (32 GB)
        N_nodes_max  = 32
        if n_cpu_per_node == None:
            n_cpu_per_node = 12 # Intel Sandybridge cores
    elif cluster.lower() == 'sandyq':
        ram_per_node = 59.604 # GiB (64 GB)
        N_nodes_max  = 32
        if n_cpu_per_node == None:
            n_cpu_per_node = 16 # Intel Sandybridge cores
    elif cluster.lower() == 'greenq':
        ram_per_node = 119.209 # GiB (128 GB)
        N_nodes_max  = 8
        if n_cpu_per_node == None:
            n_cpu_per_node = 32 # Intel Skylake cores
    elif cluster.lower() == 'starq':
        ram_per_node = 953.674 # GiB (1 TB)
        N_nodes_max  = 12
        if n_cpu_per_node == None:
            n_cpu_per_node = 128 # AMD cores

    # Need to add support for other systems
    # like maybe mac and linux personal computers

    else:
        raise ValueError( justify_text('I don\'t understand cluster "{0}".'.format(cluster),
                          first_line_buffer='ValueError' ) )
    
    # Determine the number of 3D early universe fields that will be sotred (note that this is
    # independent of the non-Gaussianity model because we never have more than one early-universe
    # scalar field in RAM when doing non-Gaussian initial conditions calculations)
    if   iLPT == 1: N_fields = 4
    elif iLPT == 2: N_fields = 7
    else:
        raise ValueError( justify_text( ('Only 1st and 2nd order Lagrangian perturbation theory ar'+
                    'e supported. iLPT must be either 1 or 2, not {0}.').format(iLPT),
                first_line_buffer='ValueError' ) )

    # Set some geometry variables
    r_latt     = boxsize/n                  # the lattice spacing in Mpc
    min_n_buff = int( r_buff/r_latt + 0.5 ) # minumum buffer thickness (in voxels)

    # Determine m_mesh and n_buff based on other parameters
    get_n_mesh = lambda n, n_buff, n_tile: (n+2*n_buff*(n_tile-1))/n_tile
    get_n_buff = lambda n_mesh, n_buff, n_tile, r_buff, boxsize: int(
        r_buff * (n_mesh-2*n_buff) * n_tile / boxsize + 0.5 )
    
    # Check if a number is prime (FFTW is very slow for primes so primes are to be avoided)
    def is_prime(a): return all( a % i for i in range( 2 , int(a/2)+1 ) )

    # Good values of n_mesh, n_tile and n_buff
    n_mesh = []
    n_tile = []
    n_buff = []
    found  = 0

    # Calculate roughly how many 32-bit floats fit per CPU on the computer architecture
    nmax = int( ( 0.5*ram_per_node  # 1/2 RAM in GiB per node
                  /n_cpu_per_node   # / number of CPUs per node
                  *2**30            # * 2^30 B/GiB
                  /4                # / 4 B per 32-b float
                  /(N_fields+1.25)  # / ( 7 fields + 1 working field + 1/4 extra )
                )**(1./3)         # ]^{1/3} to convert from volume to length
                +0.5)             # rounded to the nearest integer

    # Not much point running tiles smaller than nmesh = 199
    nmin = 199
    if n < nmin:
        raise ValueError( justify_text( ('Potential parallelization configurations for a box with '+
                    'n={0}, boxsize={1} Mpc and buffersize={2} Mpc require that integers n_mesh, n'+
                    '_buff and n_tile be found such that\n    n = (nmesh-2nbuff)*ntile+2nbuff\nAnd'+
                    ' for the {3} cluster {4} < nmesh <= {5}, which means that n must be less than'+
                    ' {4}.').format( n, boxsize, r_buff, cluster, nmin, nmax ),
                first_line_buffer='ValueError' ) )
    # nmin and nmax are used to find tiling in which the simulation volume is divided into cubic
    # tiles of side length n_mesh voxels (including buffers) with n_mesh <= nmax

    ################################################################################################
    # Find possible values for discretisation and parallelisation                                  #
    ################################################################################################

    # Find possible values for n_tile, n_buff, n_mesh for specified n and boxsize, starting with
    # n_buff = minumum n_buff value, increment up to this value plus 5
    try_n_buff = min_n_buff
    for j in range(5):

        # start by trying a single tile, increment up to 50 tiles
        try_n_tile = 1
        for i in range(50):

            # estimate n_mesh from n, current guess at n_buff and current guess at n_tile
            try_n_mesh = get_n_mesh( n, try_n_buff, try_n_tile )

            # if these values of n_mesh, n_buff are valid, add them to the list of possible values
            if ( try_n_mesh.is_integer() == True and try_n_mesh < nmax and try_n_mesh > nmin and 
                 is_prime(int(try_n_mesh)) == False ):
                if try_n_buff > get_n_buff( try_n_mesh, try_n_buff, try_n_tile, r_buff, boxsize ):

                    # round n_mesh down to an integer
                    try_n_mesh = int(try_n_mesh)

                    # Add set of values to list of possible run configurations
                    n_mesh.append( try_n_mesh )
                    n_tile.append( try_n_tile )
                    n_buff.append( try_n_buff )
                    found += 1
            try_n_tile += 1
        try_n_buff += 1

    # If there were no values found, raise a value error
    if found == 0:
        raise ValueError( justify_text( ('No potential parallelization configurations could be fou'+
                    'nd with n={0}, boxsize={1} Mpc and buffersize={2} Mpc. Parallelization requir'+
                    'es that integers nmesh, nbuff and ntile be found such that\n    n = (nmesh-2n'+
                    'buff)*ntile+2nbuff\nwhere\n    nbuff >= (nmesh-2nbuff)*ntile*buffersize/boxsi'+
                    'ze\nand for the {3} cluster, {4} < nmesh < {5}. Any of these conditions not b'+
                    'eing met could cause this error.').format( 
                        n, boxsize, r_buff, cluster, nmin, nmax ),
                first_line_buffer='ValueError' ) )

    # Print possible values
    print( justify_text( ('\nTable 1: tile and buffer geometries that parallelize a Peak Patch run'+
                ' with n={0}, boxsize={0} Mpc, and r_buff={1} Mpc. Note that this is purely geomet'+
                'rical, not accounting for computer limitations (see table 2 for more).\n').format(
                n, boxsize, r_buff ) ) )
    print(' n_tile | n_mesh | n_buff \n--------+--------+--------')
    for i in range(len(n_mesh)):
        print( '{0:7} |{1:7} |{2:7}  '.format( n_tile[i], n_mesh[i], n_buff[i] ) )

    ################################################################################################
    # Check if cluster architecture can support candidate discretisation/parallelisation schema    #
    ################################################################################################

    # Assume a halo catalogue has at most 1 million halos in it
    N_halos_max = 1.0e6
    
    # Memory of each code section:
    def fftw_mem( n, N_nodes ):
        # Returns memory in GiB/node needed to perform FFTs on N_fields (=7 for 2LPT) fields of
        # dimension nxnxn 32-bit floating point values spread over nnodes computational nodes.
        return(   N_fields  # 7 fields per tile
                * 4         # * 4 B per 32-b float
                / 2**30     # / 2^{30} B per GiB
                * n**3      # * n^3 32-b floats per field file
                / N_nodes ) # / nnodes nodes on cluster
    
    def grid_mem( n_mesh, n_cpu_per_node ):
        # Returns the memory in GiB needed to store and work on fields (usually 7 for 2LPT + 1 
        # working field + 1/4 extra space) for a serial Peak Patch run with ngrid^3 32-bit floating
        # point values times the number of CPUs per node.
        return(   n_cpu_per_node # number of CPUs per node
                * n_mesh**3      # * number of 32-b floats per parallelization "tile"
                * (   N_fields   # *( 7 fields
                    + 1          #    + 1 working field
                    + .25 )      #    + 1/4 field for extra space )
                * 4              # * 4 B per 32-b float
                / 2**30 )        # / 2^{30} B per GiB
    
    def cat_mem( N_halos, n_cpu_per_node ):
        # Returns the memory in GiB needed to store catalogue files (which will have maximum size 
        # N_halos x cat_len where cat_len is 33 for runs where we output shear data and 11 for runs
        # without shear) times the number of CPUs per node.
        return (   n_cpu_per_node # number of CPUs per node
                 * cat_len        # * number of 32-b floats per peak in catalogue file
                 * N_halos        # * max number of peaks per catalogue file
                 * 4              # * 4 B per 32-b float
                 / 2**30 )        # / 2^{30} of B per GiB
        # Note that in some legacy versions of Peak Patch, cat_len was either 11 or 23, so you may
        # see 23 pop up randomly in some places.
    
    def lat_mem( n_buff, n_cpu_per_node ):
        # Returns the memory in GiB needed to perform peak-finding on lattice.
        return (   n_cpu_per_node      # number of CPUs per node
                 * 24                  # * 24 parameters per voxel in peak finding
                 * 4/3*np.pi*n_buff**3 # * volume (in voxels) of largest peaks
                 * 4                   # * 4 B per 32-b float
                 / 2**30 )             # / 2^{30} of B per GiB
    
    def s2c_mem( n, n_mesh, N_nodes, n_cpu_per_node ):
        # Returns the memory in GiB used to do the slab-to-cube transformation
        return ( ( # number of voxels in a tile:
                     n_mesh**3      # n_mesh^3
                   # number of voxels in a slab:
                   + n              # sidelength of full cube /w buff in voxels
                   / N_nodes        # / number of nodes                         
                   / n_cpu_per_node # / number of CPUs per node                 
                   * n_mesh**2 )    # n_mesh^2                                         
                 * n_cpu_per_node   # * number of CPUs per node
                 * 4                # * 4 B per 32-b float
                 / 2**30 )          # / 2^{30} B per GiB

    # try all possible n_tile, n_buff, n_mesh values, see which combination splits up best for
    # cluster architecture
    if   iLPT==1: iLPT_string='1st order Lagrangian perturbation theory'
    elif iLPT==2: iLPT_string='2nd order Lagrangian perturbation theory'
    print( justify_text( ('\nTable 2: Peak Patch run parameters with n={0}, boxsize={1} Mpc, r_buf'+
                'f={2} Mpc producing a halo catalogue with {3} 32-bit numbers per halo, using {4} '+
                'with {5} processors/node on {6}.\n').format( n, boxsize, r_buff, cat_len, 
                iLPT_string, n_cpu_per_node, cluster ) ) )
    print( ('     Peak Patch run parameters     |               memory in GiB for each step\n-----'+
            '------------------------------+------------------------------------------------------'+
            '---\n N_nodes| n_tile | n_mesh | n_buff |  FFTW  | initial|  halo  |  peak  |  slab  '+
            '| total\n        |        |        |        |        | fields |catalogs| finding| to '+
            'cube| /{0} GiB\n--------+--------+--------+--------+--------+--------+--------+------'+
            '--+--------+------------').format(ram_per_node) )
    N_nodes = np.empty( (0,4), dtype=np.int32 )
    any_check_ram  = False # Switched to True if a configuration with memory less than noderam is found
    any_check_proc = False # Switched to True if a configuration with an integer n/nnodes/nproci found
    for i in range(len(n_tile)):

        try_n_tile = n_tile[i] # number of tiles
        try_n_buff = n_buff[i] # n_buff
        try_n_mesh = n_mesh[i] # n_mesh
    
        # Number of processes is whichever is less, the number of tiles n_tile^3 or the number of
        # CPUs per node nproc
        try_n_cpu_per_node = min( try_n_tile**3, n_cpu_per_node )
    
        try_N_nodes = 1     # start with one node and count up
        min_reached = False # flag
        while not min_reached and try_N_nodes <= N_nodes_max:
            fftw = fftw_mem( n, try_N_nodes )
            grid = grid_mem( try_n_mesh , try_n_cpu_per_node )
            cat = cat_mem( N_halos_max, try_n_cpu_per_node )
            lat = lat_mem( try_n_buff, try_n_cpu_per_node )
            s2c = s2c_mem( n, try_n_mesh, try_N_nodes, try_n_cpu_per_node )
            total_mem = fftw+grid+cat+lat+s2c # Total memory in GiB
 
            # Check that the memory requirement is less than the total memory per node of cluster
            if total_mem < ram_per_node:
                any_check_ram = True  # Will remain set to True if any run passes
                check_ram     = True  # True if this run passes
            else:
                check_ram     = False # False if this run doesn't pass

            # Check that the number of parallel processes for the FFTW slab decomposition can be
            # evenly distributed across the allotment of parallel processors
            if ( n / try_N_nodes / try_n_cpu_per_node ).is_integer():
                any_check_proc = True  # Will remain set to True if any run passes
                check_proc     = True  # True if this run passes
            else:
                check_proc     = False # False if this run doesn't pass
 
            # If the previous two conditions are met, the run is allowed, print details to screen
            if ( check_ram and check_proc ):

                print( ('{0:7} |{1:7} |{2:7} |{3:7} |{4:7.5g} |{5:7.5g} |{6:7.5g} |{7:7.5g} |{8:7.'+
                        '5g} |{9:11.9g}').format(try_N_nodes,try_n_tile,try_n_mesh,try_n_buff,
                                                fftw,grid,cat,lat,s2c,total_mem) )
                min_reached = True
                
                # Add to list of geometries that does not overload RAM and CPU limits
                N_nodes = np.vstack([ N_nodes, 
                                      np.array([ try_N_nodes, try_n_tile, try_n_buff, try_n_mesh ]) ])
            try_N_nodes+=1
    
    # If no configuration that works could be found throw a Value Error
    if not any_check_ram and any_check_proc:
        raise ValueError( ('No configuration was found for which the memory re'
            +'quired was less than the\n{0} GiB/node allowed by {1}. You lik'
            +'ely need to choose a different value\nof n with more factors so '
            +'that it will be more likely that we can find a\nresoluiton ntile'
            +' for the parallelization of the volume.'
            ).format(ram_per_node,cluster) )
    elif any_check_ram and not any_check_proc:
        raise ValueError( ('No configuration was found for which n/nnodes/npro'
            +'ci was an integer, which\nis required for the parallelization sc'
            +'heme used by FFTW. The values used\nwere\nn     ={0}\nnnodes={1}'
            +'\nnproci=min(ntile**3,nproc)\nIf n or nnodes are prime numbers, '
            +'that is likely the culprit, as the\nalgorithm tries to find ntil'
            +'e that divides the total volume evenly (note\nthat nproc for {2}'
            +' is {3}).').format( n, N_nodes_max, cluster, n_cpu_per_node ) )
    
    elif not any_check_ram and not any_check_proc:
        raise ValueError( ('No configuration was found for which the memory re'
            +'quired was less than the\n{0} GiB/node allowed by {1}, or for '
            +'which n/nnodes/nproci was an integer.\nYou likely need to choose'
            +' a different value of n with more factors so that\nit will be mo'
            +'re likely that we can find a resoluiton ntile for the paral-\nle'
            +'lization of the volume.').format( ram_per_node, cluster ) )
     
    else:
        # set list of n_mesh, n_buff, n_tile, and N_nodes
        n_mesh  = N_nodes[:,3]
        n_buff  = N_nodes[:,2]
        n_tile  = N_nodes[:,1]
        N_nodes = N_nodes[:,0]

        # use setup that requires minimum number of nodes, if multiple have the same number of
        # nodes, pick the one with the fewest tiles
        N_nodes_min = min(N_nodes)
        n_tile_min  = 999
        for i in range(len(N_nodes)):
            if N_nodes[i] <= N_nodes_min:
                N_nodes_min = N_nodes[i]
                if n_tile[i] <= n_tile_min:
                    n_tile_min = n_tile[i]
                    ind = i
        
        # Print optimum run
        print('\n===================================='            )
        print('Values optimized for number of nodes  '            )
        print('ntile  = ' , n_tile[ind]                           )
        print('nbuff  = ' , n_buff[ind]                           )
        print('nmesh  = ' , n_mesh[ind]                           )
        print('nnodes = ' , N_nodes[ind]                          )
        print('tpnode = ' , min(n_tile[ind]**3 , n_cpu_per_node ) )
        print('====================================\n'            )



# This script runs the early universe instability codes
def run_instability( ins_dir    = os.getcwd() , # the directory containing all the instability code
                     a_e        = 1.0e-50     , # scale factor at end of instability (a=1 now)
                     boxsize    = 10000.0     , # Peak Patch box size (for determining k range)
                     m_phi      = 1.0         , # m^2_\phi in units of 10^-5(8\pi)^-0.5 m_P
                     m_chi      = 1.0         , # m^2_\chi in units of 10^-5(8\pi)^-0.5 m_P
                     m2_phi     = None        , # convenience variable to specify m^2 instead of m
                     m2_chi     = None        , # convenience variable to specify m^2 instead of m
                     phi_w      = 0.12547     , # \phi_w width of instability in \phi
                     phi_p      = 8.49953     , # \phi_p value of \phi at centre of instability
                     m_tach     = 25.0        , # m_\lambda instability strength
                     vev        = 0.1         , # vev vacuum expectation value
                     nstep      = 2**15       , # maximum number of steps to converge
                     stepadapt  = 1           , # number of steps between each time step
                     dt0        = 2**-5       , # base time step in integrator 
                     nk         = 100         , # number of k modes in chi power
                     k0         = None        , # fundamental k mode
                     kn         = None        , # nyquist k mode (largest k in chi power spectrum)
                     save_out   = False       , # save raw output from instability codes
                     overwrite  = False       , # True to overwrite previous instability calculation
                     ng_model   = 10          , # Non-Gaussianity model
                     strength   = -1          ):

    # The default formatted instability directory in Peak Patch
    pp_ins_dir = peak_patch_dir+'/src/instability/'

    # Check if the directory ins_dir exists, if not, make it, if it does exist and overwrite is True
    # then it will be overwritten
    if not os.path.exists(ins_dir):
        os.system( 'mkdir -p {0}'.format( ins_dir ) )
    if os.path.exists( ins_dir ):
        if overwrite == True:
            os.system( 'rm -rf {0};mkdir -p {0}'.format( ins_dir ) )
    table_dir = '{0}/../../tables'.format( ins_dir )
    if not os.path.exists( table_dir ):
        table_dir = '{0}/tables/'.format( ins_dir )
        os.system( 'mkdir -p {0}'.format( table_dir ) )

    # Copy the default formatted instability directory from Peak Patch without overwriting
    os.system( 'cp -n -r {0}/. {1}/.'.format( pp_ins_dir, ins_dir ) )

    # Set scalar field masses (note that squared masses override masses)
    if m2_phi == None:
        m2_phi = m_phi**2
    else:
        m_phi  = m2_phi**.5
    if m2_chi == None:
        m2_chi = m_chi**2
    else:
        m_chi  = m2_chi**.5

    # Unit conversions from Early Universe to LSS codes
    def kj(k0,kn,j,n): return k0*(kn/k0)**(j/n)
    lcode        = 2.6259e-52 # Mpc/(10^5 reduced planck masses)
    acode_approx = 3          # the approximate early u scale factor
    if k0  == None: k0   = 2    * np.pi / a_e * acode_approx * lcode / boxsize / 10 
    if kn  == None: kn   = 2001 * np.pi / a_e * acode_approx * lcode / boxsize * 10 
    dk   = kj(k0,kn,-1,nk)
    kmax = kj(k0,kn,nk+1,nk)

    # Print the parameters for ease of interpretting
    out_params = [ a_e, m2_phi, m2_chi, phi_w, phi_p, m_tach, vev, nstep, stepadapt, dt0, nk,k0,kn ]
    with open(ins_dir+'/instability_params.txt','w') as ins_param_file:
        for j in out_params:
            ins_param_file.write( '{0}\n'.format(j) )

    ################################################################################################
    ### Original PING Model                                                                      ###
    ################################################################################################
    if ng_model <= 10:

        # Use scipy's curve fitting routine
        from scipy.optimize import curve_fit

        ############################################################################################
        ### RUN THE V_0 CASE                                                                     ###
        ############################################################################################

        # Make a copy of the Peak Patch macro
        os.system('cd {0};cp peak-patch_macros.h temp.h'.format(ins_dir))

        # Calculate k and P_{\chi\chi}(k) in natural units \hbar=c=1, save the V_0 case power
        p_chichi_V_0 = run_evolve_corr(
                ins_dir    = ins_dir   ,   m2_phi    = m2_phi      ,   m2_chi  = m2_chi   ,
                phi_p      = phi_p     ,   phi_w     = phi_w       ,   vev     = vev      ,
                nk         = nk        ,   dk        = dk          ,   kmax    = kmax     ,
                nstep      = nstep     ,   stepadapt = stepadapt   ,   dt0     = dt0      ,
                lambda_chi = 0.0       ,   phi_init  = None        ,   phi_fin = None     ,
                set_macros = True      ,   save_out  = save_out                            )
        np.savetxt(table_dir+'/p_chichi_V0.txt',p_chichi_V_0)

        ############################################################################################
        ### RUN THE V_0 + \Delta V CASE                                                          ###
        ############################################################################################

        # Replace the copy of peak-patch_macros.h with lambda_chi = 0 used in the V_0 case
        os.system('cd {0};rm -f peak-patch_macros.h;mv temp.h peak-patch_macros.h'.format(ins_dir))

        # Calculate chi to phi quadratic fit
        lambda_chi = m_tach**2 * vev**-2
        dphi_of_chi, alpha_e, H_e = run_ballistic_ensemble(
                ins_dir         = ins_dir      ,   m2_phi     = m2_phi     ,
                m2_chi          = m2_chi       ,   phi_p      = phi_p      ,
                phi_w           = phi_w        ,   vev        = vev        ,
                nk              = nk           ,   dk         = dk         ,
                kmax            = kmax         ,   nstep      = nstep      ,
                stepadapt       = stepadapt    ,   dt0        = dt0        ,
                lambda_chi      = lambda_chi   ,   phi_init   = None       ,
                phi_fin         = None         ,   set_macros = True       ,
                save_ins_params = True         ,   save_out   = save_out    )

        # Fit to \Delta\phi = A \chi^2 + C
        def dphi_fit(chi,a,c): return a*chi**2+c
        popt,cov = curve_fit(dphi_fit, dphi_of_chi[:,0], dphi_of_chi[:,1])

        # Calculate k and P_{\chi\chi}(k) in natural units \hbar=c=1, save the V_0+DeltaV power
        p_chichi = run_evolve_corr(
                ins_dir    = ins_dir      ,   m2_phi    = m2_phi      ,   m2_chi  = m2_chi   ,
                phi_p      = phi_p        ,   phi_w     = phi_w       ,   vev     = vev      ,
                nk         = nk           ,   dk        = dk          ,   kmax    = kmax     ,
                nstep      = nstep        ,   stepadapt = stepadapt   ,   dt0     = dt0      ,
                lambda_chi = lambda_chi   ,   phi_init  = None        ,   phi_fin = None     ,
                set_macros = False        ,   save_out  = save_out                            )
        np.savetxt( table_dir+'/p_chichi_DeltaV.txt',p_chichi )

        ############################################################################################
        ### SUBTRACT V_0 TO GET APPROXIMATE \chi POWER                                           ###
        ############################################################################################

        # Make a subtracted power spectrum to account for V_0 case
        p_chichi[:,1] -= p_chichi_V_0[:,1]
        np.savetxt( table_dir+'/p_chichi.txt',p_chichi )

        # set k to units of Mpc^-1 at present
        p_chichi[:,0] *=   a_e*np.exp(-alpha_e)/lcode
        p_chichi[:,1] *= ( a_e*np.exp(-alpha_e)/lcode ) ** -3.

        # Take primordial <|chi|^2> to Peak-Patch formated power P(k)
        p_chichi[:,1] *= (2*np.pi)**-3

        # In this approximation, the differenced chi power spectrum should asymptote to zero, and
        # negative values can lead to NaNs, so we just set P to 0 for any value after the dropoff
        after_cutoff = False
        for j in range(len(p_chichi[:,1])):
            if after_cutoff == True:
                p_chichi[j,1] = 0.0
            else:
                if p_chichi[j,1] < 0.0:
                    p_chichi[j,1] = 0.0
                    after_cutoff = True

        # Save chi power spectrum and wave numbers as Fortran-ordered 32-bit floats in an
        # unformatted binary file to be read in by pktable.f90
        p_chichi.astype(np.float32).tofile( table_dir+'/p_chichi.dat' )

        return alpha_e , H_e , popt[0]

    ################################################################################################
    ### Original PING Model                                                                      ###
    ################################################################################################
    if ng_model == 11:

        # Set the coupling constant
        lambda_chi = m_tach**2 * vev**-2

        # Calculate k and P_{\chi\chi}(k) in natural units \hbar=c=1 at the beginning of the
        # instability, k=k_i, save the results to a file
        p_chichi_i = run_evolve_corr(
                ins_dir    = ins_dir      ,   m2_phi    = m2_phi      ,   m2_chi  = m2_chi        ,
                phi_p      = phi_p        ,   phi_w     = phi_w       ,   vev     = vev           ,
                nk         = nk           ,   dk        = dk          ,   kmax    = kmax          ,
                nstep      = nstep        ,   stepadapt = stepadapt   ,   dt0     = dt0           ,
                lambda_chi = lambda_chi   ,   phi_init  = None        ,   phi_fin = phi_p+phi_w   ,
                set_macros = True         ,   save_out  = True                                     )
        np.savetxt( table_dir+'/p_chichi_p.txt', p_chichi_i )
        os.system( 'cd {0};mv evolve_correlation_DeltaV.txt evolve_correlation_i.txt'
                   .format(table_dir) )

        # Calculate k and P_{\chi\chi}(k) in natural units \hbar=c=1 at the middle of the
        # instability, k=k_p, save the results to a file
        p_chichi_p = run_evolve_corr(
                ins_dir    = ins_dir      ,   m2_phi    = m2_phi      ,   m2_chi  = m2_chi   ,
                phi_p      = phi_p        ,   phi_w     = phi_w       ,   vev     = vev      ,
                nk         = nk           ,   dk        = dk          ,   kmax    = kmax     ,
                nstep      = nstep        ,   stepadapt = stepadapt   ,   dt0     = dt0      ,
                lambda_chi = lambda_chi   ,   phi_init  = None        ,   phi_fin = phi_p    ,
                set_macros = True         ,   save_out  = True                                )
        np.savetxt( table_dir+'/p_chichi_p.txt', p_chichi_p )
        os.system( 'cd {0};mv evolve_correlation_DeltaV.txt evolve_correlation_p.txt'
                   .format(table_dir) )

        # Calculate k and P_{\chi\chi}(k) in natural units \hbar=c=1 at the end of the instability,
        # k=k_e, save the results to a file
        p_chichi_e = run_evolve_corr( 
                ins_dir    = ins_dir      ,   m2_phi    = m2_phi      ,   m2_chi  = m2_chi   ,
                phi_p      = phi_p        ,   phi_w     = phi_w       ,   vev     = vev      ,
                nk         = nk           ,   dk        = dk          ,   kmax    = kmax     ,
                nstep      = nstep        ,   stepadapt = stepadapt   ,   dt0     = dt0      ,
                lambda_chi = lambda_chi   ,   phi_init  = None        ,   phi_fin = None     ,
                set_macros = True         ,   save_out  = True                                )
        np.savetxt( table_dir+'/p_chichi_e.txt', p_chichi_e )
        os.system( 'cd {0};mv evolve_correlation_DeltaV.txt evolve_correlation_e.txt'
                   .format(table_dir) )

        # Read in alpha and H
        alpha_i, H_i = np.loadtxt('{0}/evolve_correlation_i.txt'.format(table_dir))[0,:2]
        alpha_p, H_p = np.loadtxt('{0}/evolve_correlation_p.txt'.format(table_dir))[0,:2]
        alpha_e, H_e = np.loadtxt('{0}/evolve_correlation_e.txt'.format(table_dir))[0,:2]

        # alpha_e should be 0, but may be +/- a very small number, so we set it to 0
        alpha_i = 0.0

        # Interpollate the chi power at H_i and H_e respectively with k = a_i H_i and k = a_e H_e
        k_i = np.exp(alpha_i) * H_i
        k_e = np.exp(alpha_e) * H_e
        P_i = loglog_interp( k_i, p_chichi_i[:,0], p_chichi_i[:,1] )
        P_e = loglog_interp( k_e, p_chichi_e[:,0], p_chichi_e[:,1] )
        
        # Set the instability strength
        Q_code = np.sqrt( P_e/P_i )
        if strength != -1:
            Q      = strength
            C = Q**2 / Q_code**2
            p_chichi_e[:,1] *= C
        else:
            strength = Q_code

        # set k to units of Mpc^-1 at present
        p_chichi_i[:,0] *=   a_e*np.exp( -alpha_i ) / lcode
        p_chichi_i[:,1] *= ( a_e*np.exp( -alpha_i ) / lcode ) ** -3.
        p_chichi_p[:,0] *=   a_e*np.exp( -alpha_p ) / lcode
        p_chichi_p[:,1] *= ( a_e*np.exp( -alpha_p ) / lcode ) ** -3.
        p_chichi_e[:,0] *=   a_e*np.exp( -alpha_e ) / lcode
        p_chichi_e[:,1] *= ( a_e*np.exp( -alpha_e ) / lcode ) ** -3.

        # Take primordial <|chi|^2> to Peak-Patch formated power P(k)
        p_chichi_i[:,1] *= (2*np.pi)**-3
        p_chichi_p[:,1] *= (2*np.pi)**-3
        p_chichi_e[:,1] *= (2*np.pi)**-3

        # Find the horison scale k_p at the middle of the instability in LSS units
        k_instability  = np.exp(alpha_p) * np.sqrt( m_tach**2 - m_chi**2 )
        k_instability *= a_e*np.exp( -alpha_e ) / lcode

        # Make the chi smoothing kernel W(k)
        W_chi = np.heaviside( 2 * k_instability - p_chichi_e[:,0] , 0 )

        # The power of the smoothed chi field is P_{(W*\chi)^2}(k) = W^2(k) P_{\chi\chi}(k)
        p_chichi_e[:,1] *= W_chi**2

        # Save chi power spectrum and wave numbers as Fortran-ordered 32-bit floats in an
        # unformatted binary file to be read in by pktable.f90
        p_chichi_i.astype(np.float32).tofile(table_dir+'/p_chichi_i.dat')
        p_chichi_p.astype(np.float32).tofile(table_dir+'/p_chichi_p.dat')
        p_chichi_e.astype(np.float32).tofile(table_dir+'/p_chichi_e.dat')

        # Make the (W*\chi)^2 to \Delta\phi transfer function
        T_Wchi2phi       = np.zeros( ( len(p_chichi_e[:,0]), 2 ) )
        T_Wchi2phi[:,0]  = p_chichi_e[:,0]
        T_Wchi2phi[:,1]  = np.heaviside( 4 * k_instability - p_chichi_e[:,0] , 0 )

        # Save filter and transfer function
        T_Wchi2phi.astype(np.float32).tofile(table_dir+'/T_Wchi2phi.dat')
    
        # Return
        return (alpha_i,alpha_p,alpha_e) , (H_i,H_p,H_e) , strength



def run_evolve_corr( ins_dir, m2_phi=1.0, m2_chi=1.0, phi_p=8.0, phi_w=0.1, vev=0.1, nk=1000,
        dk=3.434760, kmax=3434.760, nstep=32768, stepadapt=1, dt0=0.03125, lambda_chi=0.0,
        phi_init=None, phi_fin=None,  m_lambda=None, m_tach=None, set_macros=False,
        save_out=False ):

    # Both "m_lambda" and the older term "m_tach" are used in different places in Peak Patch, so for
    # convenience, you can set the effective mass term using either (m_lambda overrides m_tach)
    if m_lambda == None and not m_tach == None:
        m_lambda = m_tach

    # Set lambda =_chi from m_lambda and vev if m_lambda is presesnt
    if not m_lambda == None:
        lambda_chi = m_lambda**2 * vev**-2

    # Set parameters in macros
    if set_macros == True:
        set_instability_macros( ins_dir=ins_dir, m2_phi=m2_phi, m2_chi=m2_chi, phi_p=phi_p,
                phi_w=phi_w, vev=vev, nk=nk, dk=dk, kmax=kmax, nstep=nstep, stepadapt=stepadapt,
                dt0=dt0, lambda_chi=lambda_chi, phi_init=phi_init, phi_fin=phi_fin,
                m_lambda=m_lambda, m_tach=m_tach, replace_macro_file=True )

    # Run correlation matrix code to get P_{\chi\chi}(k) in the V_0 case
    makeinstlog = ins_dir+'/log.log'
    os.system( ('cd {0};module load intel;make clean -f Makefile_corr >> {1};make -f Makefile_'+
                    'corr >> {1};./corr_test').format(ins_dir,makeinstlog) )

    # Save output files
    if save_out:

        # Check that tables directory exists
        if os.path.isdir( '{0}/../../tables'.format(ins_dir) ):
            save_dir = '{0}/../../tables'.format(ins_dir)
        else:
            save_dir = ins_dir

        # Save evolve correlation raw output for the V0 case
        if lambda_chi == 0:
            os.system( 'cp {0}/corr.out {1}/evolve_correlation_V0.txt'.format(ins_dir,save_dir) )

        # Save evolve correlation raw output for the V0 + DeltaV case
        else:
            os.system( 'cp {0}/corr.out {1}/evolve_correlation_DeltaV.txt'.format(ins_dir,save_dir))

    # Read in k and P_{\chi\chi}(k) in natural units \hbar=c=1
    p_chichi = np.loadtxt(ins_dir+'/corr.out',usecols=(-17,-6))
    os.system('cd {0};make clean -f Makefile_corr &> /dev/null'.format(ins_dir))

    return p_chichi



def run_ballistic_ensemble( ins_dir, m2_phi=1.0, m2_chi=1.0, phi_p=8.0, phi_w=0.1, vev=0.1, nk=1000,
        dk=3.434760, kmax=3434.760, nstep=32768, stepadapt=1, dt0=0.03125, lambda_chi=0.0,
        phi_init=None, phi_fin=None,  m_lambda=None, m_tach=None, set_macros=False,
        save_ins_params=True, save_out=False ):

    # Both "m_lambda" and the older term "m_tach" are used in different places in Peak Patch, so for
    # convenience, you can set the effective mass term using either (m_lambda overrides m_tach)
    if m_lambda == None and not m_tach == None:
        m_lambda = m_tach

    # Set lambda =_chi from m_lambda and vev if m_lambda is presesnt
    if not m_lambda == None:
        lambda_chi = m_lambda**2 * vev**-2

    # Set parameters in macros
    if set_macros == True:
        set_instability_macros( ins_dir=ins_dir, m2_phi=m2_phi, m2_chi=m2_chi, phi_p=phi_p,
                phi_w=phi_w, vev=vev, nk=nk, dk=dk, kmax=kmax, nstep=nstep, stepadapt=stepadapt,
                dt0=dt0, lambda_chi=lambda_chi, phi_init=phi_init, phi_fin=phi_fin,
                m_lambda=m_lambda, m_tach=m_tach, replace_macro_file=True )

    # Run ballistic ensemble code to get \chi->\Delta\phi transfer function and read it in
    makeinstlog = ins_dir+'/log.log'
    os.system( ('cd {0};module load intel;make clean -f Makefile_ens &> /dev/null;make -f '+
                'Makefile_ens > {1}').format(ins_dir,makeinstlog) )
    dphi_of_chi = np.loadtxt(ins_dir+'/ball_ens.out',usecols=(4,3,1,2))

    # Get lattice scale factor at end of instability for \chi=0 (typically \alpha(\chi=0)~1.1 where
    # a=e^\alpha=1 at the start of the instability)
    alpha_e = dphi_of_chi[0,2]
    print('Primordial scale factor a = ',np.exp(alpha_e))

    # Get lattice Hubble parameter at the end of the instability
    H_e = dphi_of_chi[0,3]
    print('Primordial Hubble parameter H_e = ',H_e)

    # Write scale factor and Hubble parameter to the parameter file
    if save_ins_params == True:
        with open(ins_dir+'/instability_params.txt','a') as ins_param_file:
            ins_param_file.write( '{0}\n{1}\n'.format( np.exp(alpha_e) , H_e ) )

    # Save ballistic ensemble raw output
    if save_out:

        # Check that tables directory exists
        if os.path.isdir( '{0}/../../tables'.format(ins_dir) ):
            save_dir = '{0}/../../tables'.format(ins_dir)
        else:
            save_dir = ins_dir

        os.system( 'cd {0};cp ball_ens.out {1}/ballistic_ensemble.txt'.format(ins_dir,save_dir) )

    # Cleanup
    os.system( 'cd {0};make clean -f Makefile_ens &> /dev/null'.format(ins_dir) )

    return dphi_of_chi, alpha_e, H_e



def set_instability_macros( ins_dir, m2_phi=1.0, m2_chi=1.0, phi_p=8.0, phi_w=0.1, vev=0.1,
        nk=1000, dk=3.434760, kmax=3434.760, nstep=32768, stepadapt=1, dt0=0.03125,
        lambda_chi=0.0, phi_init=None, phi_fin=None, m_lambda=None, m_tach=None,
        replace_macro_file=False ):

    # Both "m_lambda" and the older term "m_tach" are used in different places in Peak Patch, so for
    # convenience, you can set the effective mass term using either (m_lambda overrides m_tach)
    if m_lambda == None and not m_tach == None:
        m_lambda = m_tach

    # Set lambda =_chi from m_lambda and vev if m_lambda is presesnt
    if not m_lambda == None:
        lambda_chi = m_lambda**2 * vev**-2

    # Automatically set initial and final inflaton mean if not specified
    if phi_init == None:
        phi_init = phi_p + phi_w
    if phi_fin == None:
        phi_fin = phi_p - phi_w

    # Replace macro file with default from Peak Patch directory
    macro_file = 'peak-patch_macros.h'
    if replace_macro_file == True:
        os.system( 'cd {0};rm -f {1};cp {2}/src/instability/{1} ./{1}'
                   .format( ins_dir, macro_file, peak_patch_dir ) )

    # Write run parameters from Peak Patch parameter file into instability macro scripts
    s = 'cd '+ins_dir+";sed 's/{0} {1}/{0} {2}/g' {3} > temp;mv temp {3}"
    os.system( s.format( '__M2_PHI__'     , '1.0'      , float( m2_phi     ), macro_file ) )
    os.system( s.format( '__M2_CHI__'     , '1.0'      , float( m2_chi     ), macro_file ) )
    os.system( s.format( '__PHI_P__'      , '8.49953'  , float( phi_p      ), macro_file ) )
    os.system( s.format( '__PHI_W__'      , '0.12547'  , float( phi_w      ), macro_file ) )
    os.system( s.format( '__VEV__'        , '0.1'      , float( vev        ), macro_file ) )
    os.system( s.format( '__PHI_INIT__'   , '8.625'    , float( phi_init   ), macro_file ) )
    os.system( s.format( '__PHI_FIN__'    , '8.37406'  , float( phi_fin    ), macro_file ) )
    os.system( s.format( '__NK__'         , '1000'     , int(   nk         ), macro_file ) )
    os.system( s.format( '__DK__'         , '3.434760' , float( dk         ), macro_file ) )
    os.system( s.format( '__KMAX__'       , '3434.769' , float( kmax       ), macro_file ) )
    os.system( s.format( '__NSTEP__'      , '32768'    , int(   nstep      ), macro_file ) )
    os.system( s.format( '__STEPADAPT__'  , '1'        , int(   stepadapt  ), macro_file ) )
    os.system( s.format( '__DT0__'        , '0.03125'  , float( dt0        ), macro_file ) )
    os.system( s.format( '__LAMBDA_CHI__' , '62500.0'  , float( lambda_chi ), macro_file ) )



def powerspectrum_create( outfile=None, z=0.0, k_min=5.0e-6, k_max=5.0e3, nkpoints=1000, **kwargs ):
    # This function uses the python inerface for CAMB (Code for Anisotropies in the Microwave
    # Background) to create a matter power spectrum and transfer functions for a Peak Patch run. For
    # more information on CAMB, read the CAMB docs: https://camb.readthedocs.io/en/latest/.
    #
    # Arguments are
    #     outfile : str or None
    #         The path at which to save the power spectra. The default is None, in this case a file
    #         called power<i>.dat in the the tables subdirectory of the Peak Patch home direcotry is
    #         used (where <i> is nothing if no file power.dat exists or otherwise the lowest integer
    #         for which the file does not already exist.
    # 
    #     z : float
    #         The redshift at which to calculate the power spectrum. The default is 0.0 because that
    #         is what Peak Patch is designed to use (applying the linear growth factor D(t) to 
    #         evolve the density field back in time to some time after the CMB emission to cacluate
    #         halo collapse.
    # 
    #     k_min : float
    #         the minimum wavenumber to calculate the ... in Mpc^{1}
    # 
    #     k_max : float
    #         the maximum wavenumber to calculate the ... in Mpc^{1}
    # 
    #     nkpoints : int
    #         the number of elements in the wavenumber and power spectrum arrays.
    # 
    # Supported key word arguments and their default values:
    # 
    #     Density fractions:
    #     Omega_b      = 0.0493 : density fraction of baryons
    #     Omega_CDM    = 0.2645 : density fraction of cold dark matter
    #     Omega_m      = 0.3138 : density fraction of all nonrelativistic species that collapse
    #                             under gravity
    #     Omega_Lambda = 0.6862 : density fraction of cosmological-constant dark energy
    #     Omega_k      = 0.0    : equivalent density fraction of curvature
    #
    #     Current characteristics
    #     h        = 0.6736    : dimensionless Hubble constant "little h", H_0/(100 km/s/Mpc)
    #     sigma_8  = 0.8111    : standard deviation of field smoothed at 8 Mpc/h
    #     rho_crit = 2.7754e11 : critical energy density in units of h^2 M_sol Mpc^-3
    #     tau      = 0.0544    : optical depth of reionization
    #     m_nu     = 0.0       : sum of massive neutrino masses
    #
    #     Characteristics of the primordial power spectrum
    #     $\mathcal{P}_{\zeta\zeta}(k) = A_s (k/k_pivot)**(n_s-1)$
    #     n_s = 0.9649 : scalar spectral index
    #     A_s = 2.100e-9 : primordial power spectrum amplitude
    # 
    # Default values for cosmological model parameters are taken from to "Planck 2018 results. VI
    # Cosmological Parameters", Table 2, column TT,TE,EE+lowE+lensing.
    # 
    # Notes on formatting
    # Peak Patch initial conditions fields are generated from power spectrum and transfer
    # function files that are formatted in a slightly peculiar way. Here are a few important
    # notes on the formats used:
    # 
    # - Units
    #   Peak Patch always uses distance units of Mpc, not Mpc/h. This is a source for confusion
    #   because the LSS community very often uses Mpc/h. You can always divide out the h in the
    #   parameter file, but keep in mind that all calculations of dynamics are assuming this
    #   value of h so by just dividing h out, you're still not truly agnostic of its value.
    # 
    # - Power spectra files
    #   Files tabulating power spectra and transfer functions have some formatting querks as
    #   well, these files are generated using CAMB and the function powerspectrum_create()
    #   below, which saves them to the `peakpatch/tables/' directory. They are formatted with 4
    #   columns as follows
    #       column 1: wave number, k, [Mpc^{-1}]
    #       column 2: linear matter power spectrum at z=0 in units of Mpc^3
    #                 P_peakpatch = (2\pi)^{-3} P_{ff}(k)
    #       column 3: transfer function to primordial zeta power in units of Mpc^3
    #                 T_peakpatch = (2\pi)^{-3/2} T_{\zeta \to f}(k)
    #       column 4: linear theory power spectrum of primordial field sourcing non-
    #                 Gaussianities, P_{\chi\chi}(k), [Mpc^3] in the same format as the density
    #                 power.
    # where
    #     f = \delta \bar{\rho} = \rho - \bar{\rho}
    # is the density field minus it's average (or you can look at this as the dimensionless
    # overdensity paramter $\delta$ times the average density $\bar{\rho}$), 
    #     P_{ff}(k) = \langle | \tilde{f}(\mathbf{k}) |^2 \rangle
    # is the linear matter power spectrum in the typical form, taken as the Fourier transform of
    # $f$ averaged over spherical shells of each wavenumber $k$ (note that the Peak Patch linear
    # matter power spectrum has an extra factor of $(2\pi)^{-3}$ multiplying it), and
    #     T_{\zeta to f} = \sqrt{ P_{ff}(k) / P_{\zeta\zeta}(k) }
    # is the transfer function that relates the linear matter and primordial $\zeta$ power
    # spectra. Finally, the field $P_{\chi\chi}$ is the power spectrum for a second early 
    # universe field, either a transverse field during multi-field inflation or an effective
    # field during reheating.

    # Import CAMB
    try:
        import camb
        from camb import model, initialpower
    except:
        raise ImportError('Could not import CAMB. Exiting.')
        return

    # If no outfile is passed, automatically select one and make sure not to overwrite existing file
    if outfile == None:
        outfile = peak_patch_dir+'/tables/power{0}.dat'
        if os.path.isfile( outfile.format('') ):
            counter = 1
            while os.path.isfile( outfile.format(counter) ):
                counter += 1
            outfile = outfile.format(counter)
        else:
            outfile = outfile.format('')

    # Make sure path to outfile exists
    elif not os.path.isfile( outfile ):
        os.system( 'mkdir -p '+os.path.dirname(outfile) )

    # Check that kwargs is a dictionary
    if type(kwargs) != dict: raise TypeError('kwargs must be of type dict.')

    # Check that only two of the matter-component density fractions are provided in kwargs
    if 'Omega_b' in kwargs and 'Omega_CDM' in kwargs and 'Omega_m' in kwargs:
        raise KeyError('Since Omega_m = Omega_b + Omega_CDM, you can only specify two of these.')

    # Check that both Omega_Lambda and Omega_k are not specified
    if 'Omega_Lambda' in kwargs and 'Omega_k' in kwargs:
        raise KeyError('1 = Omega_m + Omega_k + Omega_Lambda, you cannot specify both Omega_k and '+
                        'Omega_Lambda')

    # Assign default values of kwargs
    Omega_b      = 0.0493
    Omega_CDM    = 0.2645
    Omega_m      = 0.3138
    Omega_Lambda = 0.6862
    Omega_k      = 0.0
    Omega_total  = 1.0
    h            = 0.6736
    sigma_8      = 0.8111
    rho_crit     = 2.7754e11 # h^2 M_sol Mpc^-3
    tau          = 0.0544
    m_nu         = 0.0 # eV
    n_s          = 0.9649
    A_s          = 2.100e-9
    k_pivot      = 0.05 # h / Mpc

    # Assign any values speficied as kwargs
    for key in kwargs:
        exec( '{0}={1}'.format(key,kwargs[key]) )

    # Automatically adjust matter energy density fractions
    if   'Omega_b' in kwargs and 'Omega_CDM' in kwargs: Omega_m   = Omega_b + Omega_CDM
    elif 'Omega_b' in kwargs and 'Omega_m'   in kwargs: Omega_CDM = Omega_m - Omega_b
    elif 'Omega_m' in kwargs and 'Omega_CDM' in kwargs: Omega_b   = Omega_m - Omega_CDM

    if   'Omega_lambda' in kwargs: Omega_k      = 1 - Omega_m - Omega_Lambda
    elif 'Omega_k'      in kwargs: Omega_Lambda = 1 - Omega_m - Omega_k

    # Setup CAMB linear matter power spectrum calculation by setting cosmology parameters
    H_0 = h*100
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H_0,ombh2=Omega_b*h**2,omch2=Omega_CDM*h**2,mnu=m_nu,omk=Omega_k,tau=tau)
    pars.set_dark_energy()                                 # set DE equation of state
    pars.InitPower.set_params( ns=n_s, As=A_s )            # set primordial power spectrum params
    pars.set_matter_power( redshifts = [z], kmax = k_max ) # set z and k ranges for calculation
    pars.NonLinear = model.NonLinear_none                  # linear matter power spectrum only
    pars.PK_WantTransfer = 1     # include transfer functions in power spectrum calculations
    pars.WantTransfer    = 1     # calculate and output transfer functions
    pars.Transfer.kmax   = k_max # set maximum wavenumber for transfer function calculations
  
    # Do calculations needed for power spectrum calculation, store in object of type CAMBdata
    results = camb.get_results(pars)
     
    # Calcualte matter power spectrum for CAMBdata object `results' where:
    # - kh[j] is an array of wavenumber k/h
    # - zz[i] is an array of redshift z
    # - pk[i,j] is a matrix of powerspectrum as a function of k/h and z
    kh, zz, pk = results.get_matter_power_spectrum( minkh = k_min/h, maxkh = k_max/h,
                                                    npoints = nkpoints )
  
    # Get sigma_8 (standard deviation for overdensity smoothed at 8 Mpc/h)
    # value at specified redshift z. Must be called after get_matter_power_spectrum()
    s8 = np.array(results.get_sigma8())
 
    # Calculate matter transfer function
    transfer = results.get_matter_transfer_data()

    # Store transfer function as array Tf[n,i,j], n=6 is the transfer function between linear matter
    #  power and the primordial zeta power at redshift z[j] and wavenumber k[i]
    Tf = transfer.transfer_data
    kk = transfer.q/h # comoving wavenumber k[i]/h
    Tf_m = Tf[6,:,0]  # transfer function T[i] we only calculated for one redshift so j=0
 
    # Normalize matter power spectrum by setting $\sigma_8$ equal to the specified value `sigma_8'
    norm = (sigma_8/s8)**2 # normalization constant
    k    = kh * h         # exact wavenumber k, not scaled by h
    pk   = norm * pk[0,:] / ( 2.*np.pi * h)**3 # normalized P_m(z=0,k)
    # Note here that we divide by h^3 to get pk in units of (Mpc/h)^3
 
    # Primordial zeta power spectrum
    pkzeta = 2*np.pi**2 * A_s / k**3 * (k/k_pivot)**(n_s-1)
 
    #Get transfer function
    Trans = np.sqrt(pk/pkzeta)

    # Light field power spectra for Preheating work (e.g. see arXiv:0903.3407
    Achi = (5.e-7)**2
    pkchi = 2*np.pi**2*Achi/k**3 #in units of sigmas
    pkchi = pkchi/(2*np.pi)**3 #for pp power spectra

    # Save the calculated matter power spectrum, transfer function, and preheating power spectrum
    np.savetxt( outfile , np.transpose([ k, pk, Trans, pkchi ]) , fmt='%1.4e' )



# Script for running homogeneous ellipsoidal collapse
def run_HEC( f, e_v, p_v, filename,
             fortrancompiler=None, fortrancompileroptions=None, **kwargs ):
    # This script runs the Homogeneous Ellipsoidal Collapse fortran code which is used in the Peak
    # Patch calculation to determine if a given peak in the smoothed density field is expected to
    # form a virialised structure. 
    # 
    # The strain along the three axes are determined by the overdensity $f$ of the density peak and
    # two parameters called the ellipticity $e_v$ and prolateness $p_v$. Together, these indicate
    # the anisotropy of the collapse. You might think of them as representing the difference in
    # gravitational force that a single particle placed at the same distance from the centre of the
    # overdensity along each of the three axes would feel:
    #     $ F_x = e^{ -f/3 (-3e_v+p_v) } F_g $
    #     $ F_y = e^{ +f/3 (     2p_v) } F_g $
    #     $ F_z = e^{ -f/3 (+3e_v+p_v) } F_g $
    # These three parameters along with a file at which to save the HEC results are the mandatory
    # arguments of this function. Additionally there are optional arguments that specify the fortran
    # compiler and compiler options, as well as many keyword argments that can be used to set
    # variables related to the numerical integration, cosmology, input/output and the HEC physics.
    # 
    # Arguments
    # 
    #     f : float
    #     The overdensity. The HEC calculation determines the collapse of an initially spherical
    #     overdensity in a linear density field. f is the magniude of the overdensity at its peak
    #     and relates to the prominence of the overdense feature.
    #
    #     e_v : float
    #         The ellipticity, which along with f and p_v define the strain field.
    #         e_v >= 0
    # 
    #     p_v : float
    #         The prolateness, which together with f and e_v define the strain field.
    #         -e_v <= p_v <= e_v
    # 
    #     filename : str
    #         The file at which to save the HEC output.
    # 
    # Optional arguments
    # 
    #     fortrancompiler : str or None, default None
    #         The fortran compiler to use, allowed values are
    #         - 'gfortran' for Gnu Fortran compiler
    #         - 'gcc' for Gnu Fortran installed as part of the Gnu C compiler
    #         - 'ifort' for Intel Fortran comipler
    #         - None to try and predict based on the machine you're using.
    # 
    #     fortrancompileroptions : str or None, default None
    #         This is a list of compiler options (flags) for running Fortran. The default is None,
    #         if None is passed it will try to decide what flags to use baesd on your fortran
    #         compiler. Common flags are:
    #             -0i is the optimisation with -01 being minimum and -04 being maximum
    #             -w supresses all warnings from the compiler
    # 
    # Key word arguments
    # 
    #     There are a lot of key word arguments allowed for this function, so I break them down into
    #     parameters related to the numerical integration calculation, to the parameters defining
    #     the homogeneous ellipsoidal collapse and virialisation, to the fundamental cosmological
    #     parameters, and to input and output for the function. 
    # 
    #     Numerical Integration kwargs:
    # 
    #         The numerical integration parameters are declared in the module 
    #         src/modules/GlobalVariables/params.f90.
    # 
    #         iwant_evmap : int, default 4
    #             To make additional tables of
    #             1) Virialisation redshift z_vir vs ellipticity e_v
    #             2) Virialisation redshift z_vir vs prolateness p_v for fixed e_v
    #             3) table of e_v and p_v
    #             4) Make no additional tables
    #
    #         nstepmax : int, default 1e4
    #             The maximum number of steps to take in the numerical integration calculation.
    # 
    #         iwant_rd : int, default 1
    #             If iwant_rd-1, Carlson's elliptic integrals are used, otherwise a simplified model
    #             is used.
    #
    #         tfac : float, default 0.01
    #             Factor multiplying the local 1-axis Hubble time to get a time step dt.
    # 
    #     Homogeneous Ellipsoidal Collapse and virialisation kwargs
    # 
    #         zinit : float, default 20 * f
    #             Initial redshift of the HEC calculation. This redshift should correspond to a time
    #             at which the local universe around the ellipsoid can be modelled by linear theory
    #             evolution from the post-inflation field. The default is to set it to 20 times the
    #             overdensity of the field at the point where the ellipsoid is found.
    # 
    #         dcrit : float, default 200.0
    #             The critical overdensity an ellipse must reach to be considered to form a
    #             virialised dark matter halo. Theoretically, virialisation occurs when
    #             $\rho/\rho_{crit} = 18 \pi^2 = 177.886... $. It is common convention to instead
    #             use the even value of 200, which is what is done by default here.
    # 
    #         e_vmax : float, default 0.0
    #             The maximum value of e_v in the e_v map if iwant_evmap==1.
    # 
    #         de_v : float, default 0.0
    #             The difference between e_v values in the e_v map if iwant_evmap==1.
    # 
    #         p_vbar : float, default 0.0
    #             Mean value of p_v in p_v map if iwant_evmap==2.
    # 
    #         dp_v : float, default 0.0
    #             Difference between p_v values in p_v map if iwant_evmap==2.
    # 
    #         e_vbar : float, default 0.0
    #             Mean value of e_v in p_v map if iwant_evmap==2.
    # 
    #         Fbar : float, default 0.0
    #             Mean value of f, only used if iwant_evmap==1 or 2
    # 
    #         fcoll_1 : float, default 0.01
    #         fcoll_2 : float, default 0.171
    #         fcoll_3 : float, default 0.171
    #             Radial freeze-out factors for each axis of the ellipsoid. If you start with a
    #             sphere of radius R, in the end you get an ellipsoid with semi-axes
    #             a,b,c = fcoll_1 * R , fcoll_2 * R , fcoll_3 * R. The default values come from two
    #             physically informed values, the virialisation collapse factor 200^{1/3} ~ 0.178
    #             (meaning that if all three axes collapsed to this factor, you get a virialised
    #             halo) and the total collapse factor 0.01 (chosen so that the value is sufficiently
    #             close to zero without causing unwanted numerical effects in ODE solvers). The
    #             defaults were chosen to best recreate the results halo mass functions (the number
    #             of halos as a function of mass) of N-body simulations.
    # 
    #         ivir_strat : int, default 2
    #             1) a_jeq=fcoll_3 a_b3, 2) a_jeq=fcoll_j a_bj
    # 
    #         iforce_strat : int, default 4
    #             0) no bg
    #             1) sbg
    #             3) bg+NLstrain
    #             4) stbg+Lstrain
    #             5) Lstrain
    #             6) SW b_i
    # 
    #     Cosmological kwargs
    #     
    #         All of the cosmological key word arguments are declared in cosmoparams.f90.
    # 
    #         The default values of all of the cosmological key word arguemnts are the Planck
    #         Collaboration's 2018 data release. These can be found in "Planck 2018 results. VI
    #         Cosmological parameters", Table 2, column TT,TE,EE+lowE+lensing.
    # 
    #         Omb : float, default 0.0493
    #             Omega_baryon, baryonic matter density fraction
    # 
    #         Omx : float, default 0.2645
    #             Omega_x, cold DM density fraction
    # 
    #         Omvac : float, default 0.6862
    #             Omega_Lambda, DE density fraction
    # 
    #         Omcur : float, default 0.0
    #             Omega_k, curvature parameter
    # 
    #         h : float, default 0.6736
    #             little h (dimensionless Hubble constant), h = H_0/(100 km s^-1 Mpc^-1)
    # 
    #         ns : float, default 0.9649
    #             n_s, spectral index
    # 
    #         sigma8 : float, default 0.8111
    #             sigma(z=0,r=8 Mpc)
    # 
    #     Input/Output parameters
    # 
    #         There is relevant one i/o parameter, declared in input_parameters.f90. It is used in
    #         formatting the output and assigned 1 in hpkvd.f90, so I did the same.
    # 
    #         ihard : int, default 1

    # The Peak Patch environment variable must be set for this script to work
    if not peak_patch_dir:
        raise SystemError('peak_patch_dir must be defined.')

    # Check that ellipticity and prolateness values make sense
    if e_v < 0:
        raise ValueError('ellipticity e_v must be greater than or equal to 0.')
    elif p_v < -e_v or p_v > e_v:
        raise ValueError('prolateness p_v must be in the range [-e_v,e_v].')

    # Try to determine what Fortran compiler to use if None passed
    if fortrancompiler == None:
        machine = subprocess.check_output('hostname').decode('utf-8')[:-1]

        # list of recognised hostnames that use Gnu Fortran compiler
        recognised_gfortran_machines = [ 'velociraptor', 'homes-MacBook-Air.local' ]

        # list of recognised hostnames that use Gnu Fortran as part of the Gnu C compiler
        recognised_gcc_machines = [ 'Nates-MacBook-Pro.local' ]

        # list of recognised hostnames that use Intel Fortran compiler
        recognised_intel_machines = [ 'mussel', 'kingcrab', 'homard', 'lobster', 'calamari',
                'shrimp', 'prawn'] + [ 'nia-login0{0}.scinet.local'.format(j) for j in range(1,9) ]

        # Set fortran compiler
        if   machine in recognised_gfortran_machines: fortrancompiler = 'gfortran'
        elif machine in recognised_gcc_machines     : fortrancompiler = 'gcc'
        elif machine in recognised_intel_machines   : fortrancompiler = 'ifort'
        else:
            raise ValueError( justify_text( ( 'Your machine is not recognised, please either speci'+
                                              'fy fortrancompiler or add your machine to the list '+
                                              'of regonised machines.'),
                                            first_line_buffer='ValueError' ) )

    # Make sure compiler is recognised
    if fortrancompiler not in ['gfortran', 'gcc', 'ifort']:
        raise SystemError('unsupported Fortran compiler '+fortrancompiler+'.')

    # Set Fortran compiler options if None specified
    if fortrancompileroptions == None:
        if   fortrancompiler == 'gfortran':
            fortranmodule          = 'gfortran'
            fortrancompileroptions = '-w'
        elif fortrancompiler == 'gcc':
            fortranmodule          = 'gcc'
            fortrancompileroptions = '-w'
        elif fortrancompiler == 'ifort':
            fortranmodule          = 'intel'
            fortrancompileroptions = '-04 -w'

    # Directories used by HomogeneousEllipsoid.f90
    dir_ = { 'homel'           : peak_patch_dir+'/src/modules/HomogeneousEllipsoid/' ,
             'src'             : peak_patch_dir+'/src/'                              ,
             'External'        : peak_patch_dir+'/src/modules/External/'             ,
             'GlobalVariables' : peak_patch_dir+'/src/modules/GlobalVariables/'      ,
             'Solvers'         : peak_patch_dir+'/src/modules/Solvers/'              ,
             'cosmology'       : peak_patch_dir+'/src/cosmology/'                     }

    # All Fortran programs that are dependencies of HomogeneousEllipsoid.f90
    f90s = { 'intreal_types'        : dir_['External']        ,
             'params'               : dir_['GlobalVariables'] ,
             'cosmoparams'          : dir_['GlobalVariables'] ,
             'input_parameters'     : dir_['GlobalVariables'] ,
             'Solvers'              : dir_['Solvers']         ,
             'Dlin_params'          : dir_['cosmology']       ,
             'psubs_Dlinear'        : dir_['cosmology']       ,
             'HomogeneousEllipsoid' : dir_['homel']           ,
             'runhomel'             : dir_['homel']            }

    # All fortran scripts to compile and link for the Homogeneous Ellipsoidal collapse calculation
    scripts_f90 = ''
    for i in f90s:
        scripts_f90 += f90s[i]+i+'.f90 '
 
    # Comment out calls to mpi_finalize(ierr) in HomogeneousEllipsoid.f90
    os.system('''
        cd {0}
        touch temp.f90
        sed 's/call mpi_finalize(ierr)/!call mpi_finalize(ierr)/' \\
            HomogeneousEllipsoid.f90 > temp.f90
        rm -f HomogeneousEllipsoid.f90
        mv temp.f90 HomogeneousEllipsoid.f90
        '''.format( dir_['homel'] ))
        # This feature is used for parallel Peak Patch runs on Niagara, but we
        # don't need it for single homel runs.
   
    # Change mode to avoid issues with Peak Patch .f90 scripts being executable
    os.system('cd {0};chmod -x {1}'.format( dir_['homel'], scripts_f90 ))
     
    # Make logfile
    os.system('rm -f {0};touch {0}'.format( filename ))
     
    # Compile and Link .f90 files
    exe_homel = 'runhomel' # Executable file name
    os.system(f'''
        cd {dir_['homel']}
        ml {fortranmodule}
        {fortrancompiler} {fortrancompileroptions} {scripts_f90} -o {exe_homel}
        ''')
   
    # HEC parameters, see documentation for full explanation of each.
    params = { 'iwant_evmap' : 4,      'nstepmax'    : 1e4,    'iwant_rd'    : 1,
               'tfac'        : 0.01,   'zinit'       : 20*f,   'dcrit'       : 200.,
               'e_vmax'      : 0.,     'de_v'        : 0.,     'p_vbar'      : 0.,
               'dp_v'        : 0.,     'e_vbar'      : 0.,     'Fbar'        : 0.,
               'fcoll_3'     : .171,   'fcoll_2'     : .171,   'fcoll_1'     : .01,
               'ivir_strat'  : 2,      'iforce_strat': 4,      'Omb'         : 0.0493,
               'Omx'         : 0.2645, 'Omvac'       : 0.6862, 'Omcur'       : 0.0,
               'h'           : 0.6736, #'ns'          : 0.9649, 'sigma8'      : 0.8111,
               'ihard'       : 1       }

    # overwrite any params in kwargs
    for param in params:
        if param in kwargs:
            params[param] = kwargs[param]

    # Make params into a list of command line arguments for the fortran executable
    params_to_read = str(len(filename)) + '\n"' + filename + '"\n'
    for i in params:
        params_to_read += str(params[i]) + '\n'

    # Run the executable runhomel
    os.system('cd {0}\n./{1} << EOF\n{2}{3} {4} {5}\nEOF'
          .format( dir_['homel'], exe_homel, params_to_read, f, e_v, p_v ))
     
    # Once execution is complete, clean up stray .mod files
    scripts_mod = ''
    for i in f90s:
        scripts_mod += i+'.mod '
    os.system('''
        cd {0}
        rm -f {1}
        '''.format( dir_['homel'], scripts_mod ))
    
    # Unfomment calls to mpi_finalize(ierr) in HomogeneousEllipsoid.f90
    os.system('''
        cd {0}
        touch temp.f90
        sed 's/!call mpi_finalize(ierr)/call mpi_finalize(ierr)/' \\
            HomogeneousEllipsoid.f90 > temp.f90
        rm -f HomogeneousEllipsoid.f90
        mv temp.f90 HomogeneousEllipsoid.f90
        '''.format( dir_['homel'] )) # So that we don't mess up Peak Patch

    # Return all the run parameters
    params['f'] , params['e_v'] , params['p_v'] = f , e_v , p_v
    return params



def plot_HEC_principle_axis_scale_factors( filename, params, ax, color='tab:blue',
        labels=[r'$a_1(z)$',r'$a_2(z)$',r'$a_3(z)$'], comoving=False,
        invert_x_axis=True, include_a_average=False, a_average_color='tab:purple',
        include_overdensity=False, overdensity_color='tab:red',
        find_z_collapse=True, z_collapse_color='black',
        z_collapse_labels=[r'$z={0}$',r'$z={0}\times10^{{{1}}}$']  ):
    # A script for plotting the output of the run_HEC script, which is a python runner for the
    # fortran Homogeneous Ellipsoidal Collapse calculation, a core piece of the Peak Patch
    # simulation.
    # 
    # Arguments
    # 
    #     filename : str
    #         The file at which to save the HEC output. This contains tabulated data describing the
    #         ellipsoidal collapse of an over-dense feature. The file has columns
    #             0. the scale factor $a$ (since the scale factor is uniformly increasing, this is 
    #                like a time coordinate). It can be4 converted into a redshift by doing
    #                $z = a^{-1} - 1$.
    #             1. a_1 principle axis scale factor for third and final axis to collapse
    #             2. a_2 principle axis scale factor for second axis to collapse
    #             3. a_3 principle axis scale factor for first axis to collapse
    #             4. 
    # 
    #     params : dict
    #         A dictionary object formatted as the output of pkp.run_HEC().
    # 
    #     ax : matplotlib axis
    # 
    # Optional arguments
    # 
    #     color : str, default 'tab:blue'
    # 
    #     labels : list of str, default [r'$a_1(z)$',r'$a_2(z)$',r'$a_3(z)$']
    # 
    #     comoving : bool, default False
    # 
    #     invert_x_axis : bool, default True
    # 
    #     include_a_average : bool, default False
    # 
    #     a_average_color : str, default 'tab:purple'
    # 
    #     include_overdensity : bool, default False
    # 
    #     overdensity_color : str, default 'tab:red'
    # 
    #     find_z_collapse : bool, default True

    # Read runhomel.log (note the call to expanduser, this is just in case the path includes '~/'
    # instead of '/User/home/' because the former will cause a crash
    raw = np.loadtxt( os.path.expanduser(filename), delimiter=',')

    # Define the redshift
    z = raw[:,0]**-1 - 1
   
    # Plot the average scale factor for the universe
    if include_a_average:
        if comoving:
            ax.plot( z, raw[:,0]*0+1, label=r'$\bar{a}$', ls='-', color=a_average_color )
        else:
            ax.plot( z, raw[:,0], label=r'$\bar{a}$', ls='-', color=a_average_color )

    # Plot of a_i(z) and delta(z)
    if comoving:
        ax.plot( z, raw[:,1]/raw[:,0], label=labels[0], ls='-',  color=color)
        ax.plot( z, raw[:,2]/raw[:,0], label=labels[1], ls='--', color=color)
        ax.plot( z, raw[:,3]/raw[:,0], label=labels[2], ls=':',  color=color)
    else:
        ax.plot( z, raw[:,1], label=labels[0], ls='-',  color=color)
        ax.plot( z, raw[:,2], label=labels[1], ls='--', color=color)
        ax.plot( z, raw[:,3], label=labels[2], ls=':',  color=color)
    if invert_x_axis: ax.invert_xaxis()
    ax.set_xlabel(r'redshift $z$')

    # Make regular y label
    if not include_overdensity:
        if comoving:
            ax.set_ylabel(r'principle axis scale factor $a_i/\bar{a}$')
        else:
            ax.set_ylabel(r'principle axis scale factor $a_i$')

    # Make coloured y label if plotting secondary axis with overdensity
    else:
        if comoving:
            ax.set_ylabel(r'principle axis scale factor $a_i/\bar{a}$', color=color )
        else:
            ax.set_ylabel(r'principle axis scale factor $a_i$', color=color )
        ax.tick_params(axis='y',colors=color)
        ax2 = ax.twinx()
        ax2.plot( z, raw[:,4]-1, label=r'$1+\delta_c(z)$', color='r')
        ax2.set_ylabel(r'overdensity $\delta$', color='r')
        ax2.set_yscale('log')
        ax2.tick_params(axis='y',colors='r')

    # Find the redshift at which each axis reaches it's collased state
    if find_z_collapse:
        z_coll_1, z_coll_2, z_coll_3 = 1100,1100,1100
        flag_1  , flag_2  , flag_3   = False,False,False
        for j in range(len(raw[:,1])-1):
            if flag_1 == False and raw[j,1] == raw[j+1,1]:
                z_coll_1 = z[j]
                flag_1   = True
            if flag_2 == False and raw[j,2] == raw[j+1,2]:
                z_coll_2 = z[j]
                flag_2   = True
            if flag_3 == False and raw[j,3] == raw[j+1,3]:
                z_coll_3 = z[j]
                flag_3   = True

        # Plot vertical line showing collapse redshift for a_1
        if flag_1 == True:
            print('axis 1 is collapsed at z = ',z_coll_1,sep='')
            base,exponent = '{0:.3e}'.format(z_coll_1).split('e')
            base,exponent = float(base), int(exponent)
            if exponent == 0:
                label = z_collapse_labels[0].format(base)
            else:
                label = z_collapse_labels[1].format(base,exponent)
            if comoving:
                ax.vlines( z_coll_1, 0, np.max(raw[:,1]/raw[:,0]), colors=z_collapse_color, ls='-',
                        label=label )
            else:
                ax.vlines( z_coll_1, 0, np.max(raw[:,1]), colors=z_collapse_color, ls='-',
                        label=label )
        else:
            print('axis 1 does not collapse by z=0')

        # Plot vertical line showing collapse redshift for a_2
        if flag_2 == True:
            print('axis 2 is collapsed at z = ',z_coll_2,sep='')
            base,exponent = '{0:.3e}'.format(z_coll_2).split('e')
            base,exponent = float(base), int(exponent)
            if exponent == 0:
                label = z_collapse_labels[0].format(base)
            else:
                label = z_collapse_labels[1].format(base,exponent)
            if comoving:
                ax.vlines( z_coll_2, 0, np.max(raw[:,1]/raw[:,0]), colors=z_collapse_color, ls='--',
                        label=label )
            else:
                ax.vlines( z_coll_2, 0, np.max(raw[:,1]), colors=z_collapse_color, ls='--',
                        label=label )
        else:
            print('axis 2 does not collapse by z=0')

        # Plot vertical line showing collapse redshift for a_3
        if flag_3 == True:
            print('axis 3 is collapsed at z = ',z_coll_3,sep='')
            base,exponent = '{0:.3e}'.format(z_coll_3).split('e')
            base,exponent = float(base), int(exponent)
            if exponent == 0:
                label = z_collapse_labels[0].format(base)
            else:
                label = z_collapse_labels[1].format(base,exponent)
            if comoving:
                ax.vlines( z_coll_3, 0, np.max(raw[:,1]/raw[:,0]), colors=z_collapse_color, ls=':',
                        label=label )
            else:
                ax.vlines( z_coll_3, 0, np.max(raw[:,1]), colors=z_collapse_color, ls=':',
                        label=label )
        else:
            print('axis 3 does not collapse by z=0')

    # Labels etc.
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper left')
    ax.set_title(r'Inititial conditions $\delta=$'+str(params['f'])+r', $e_v=$'+str(params['e_v'])+r', $p_v=$'+str(params['p_v']))

    return ax



def plot_HEC_f_ev_pv( filename, params, ax ):


    # Read runhomel.log (note the call to expanduser, this is just in case the path includes '~/'
    # instead of '/User/home/' because the former will cause a crash
    raw = np.loadtxt( os.path.expanduser(filename), delimiter=',')

    # Define the redshift
    z = raw[:,0]**-1 - 1
   
    # 
    ax.plot( z,  raw[:,4], label=r'$f(z) = \rho(z)/\bar{\rho}(z)$')
    ax.plot( z,  raw[:,5], label=r'$e_v(z)$')
    ax.plot( z,  raw[:,6], label=r'$p_v(z)$')
    ax.plot( z, -raw[:,6], label=r'$-p_v(z)$')
    
    #ax.plot( z, np.log( raw[:,1]/raw[:,3] ) / (2*raw[:,4]), label='test e')
    #ax.plot( z, np.log( raw[:,2]**2/raw[:,1]/raw[:,3] ) / (2*raw[:,4]), label='test p')
   
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.invert_xaxis()
    ax.set_xlabel(r'redshift $z$')
    ax.set_ylabel(r'eigenvalue')
    ax.legend()

    return ax
    


def plot_HEC_eigenvalues( filename, params, ax ):

    # Read runhomel.log (note the call to expanduser, this is just in case the path includes '~/'
    # instead of '/User/home/' because the former will cause a crash
    raw = np.loadtxt( os.path.expanduser(filename), delimiter=',')

    # Define the redshift
    z = raw[:,0]**-1 - 1
   
    # Eigenvalues of traceless strain tensor
    # lam1prime = -np.log( raw[:,1]*(raw[:,1]*raw[:,2]*raw[:,3])**(-1/3) )
    # lam2prime = -np.log( raw[:,2]*(raw[:,1]*raw[:,2]*raw[:,3])**(-1/3) )
    # lam3prime = -np.log( raw[:,3]*(raw[:,1]*raw[:,2]*raw[:,3])**(-1/3) )
    
    ax.plot( z, raw[:,7], label=r'$\lambda_1(z)$')
    ax.plot( z, raw[:,8], label=r'$\lambda_2(z)$')
    ax.plot( z, raw[:,9], label=r'$\lambda_3(z)$')
    # ax.plot(redshift, raw[:,7]+raw[:,8]+raw[:,9], label=r'$\delta(z)=\sum_i \lambda_i(z)$')
    # ax.plot(redshift, raw[:,4], label=r'$\delta(z)$', ls='--')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.invert_xaxis()
    ax.legend()
    plt.show()
    
    #plt.figure(7)
    #plt.plot( raw[1:,0], raw[1:,0]-raw[:-1,0], ls='none', marker='.' )
    #plt.show()

    return ax



####################################################################################################
###   Parallelization scripts   ####################################################################
####################################################################################################

def load_halos_pksc_to_npz( pksc_file, mass_cutoff, rho_m, z_range, z_chi_tab ):

    # Open halo catalogue file and load header
    c_in     = open( pksc_file , 'rb' )
    N_halos  = np.fromfile( c_in, dtype=np.int32  , count=1 )[0]
    R_th_max = np.fromfile( c_in, dtype=np.float32, count=1 )[0]
    z_obs    = np.fromfile( c_in, dtype=np.float32, count=1 )[0]

    # Load halos
    cols = int( (os.path.getsize(pksc_file)-12)/4/N_halos )
    c    = np.reshape( np.fromfile( c_in, dtype=np.float32, count=N_halos*cols ), (N_halos,cols) )

    # Perform halo mass cut
    c = c[ np.where( c[:,6] > mass_cutoff )[0] , : ] 

    # Save only the relenat data from the halo catalogue
    x ,y ,z  = c[ : ,  :3 ].T 
    dx,dy,dz = c[ : , 3:6 ].T 
    M        = c[ : ,  6  ]**3 * 4./.3*np.pi * rho_m
    zform    = c[ : , 10  ]*0 # Note: currrently not outputing the formation redshift, so it is set to zero

    # Get comoving distance, peculiar velocity of halos
    chi      = np.sqrt( x**2 + y**2 + z**2 )
    chi_inv  = np.where( chi!=0, 1/np.where( chi==0, np.nan, chi ), 0 ) 
    vpec     = ( x*dx + y*dy + z*dz ) * chi_inv
    zcos     = z_of_r_comoving( chi, z_chi_table=z_chi_tab )

    # The halo redshift is given by the sum of the redshift due to the Hubble flow and the redshift
    # due to the halo's peculiar velocity, e.g. https://arxiv.org/abs/1405.0105.
    ckms     = 299792.458 # km/s
    redshift = (1+zcos) * (1+vpec/ckms) - 1 

    # Perform redshift cut
    if z_range != (-1,-1):
        if z_range[0] == -1 and z_range[1] != -1: 
            where = np.where( ( redshift < z_range[1] ) )[0]
        elif z_range[0] != -1 and z_range[1] == -1: 
            where = np.where( ( redshift > z_range[0] ) )[0]
        else:
            where = np.where( ( redshift > z_range[0] ) & ( redshift < z_range[1] ) )[0]

        return ( x[where], y[where], z[where], dx[where], dy[where], dz[where], M[where],
                redshift[where], zform[where], zcos[where] )

    # Return the dark matter halo catalogue 
    else:
        return x, y, z, dx, dy, dz, M, redshift, zform, zcos



####################################################################################################
###   Parallelization scripts   ####################################################################
####################################################################################################

# Functions parallelized using the class Pool from the module multiprocessing must be global, not
# local, so we define them in this section instead of within the functions that they are used.

def poolfunc_pksc_to_npz( fin , fout , mass_cutoff , cosmo_header , rho_m , z_range , z_chi_tab ):
    # This function is used by PeakPatch.pksc_to_npz() to parallelize the process of making NumPy
    # formatted dark matter halo catalogues used by LIMLAM mocker from the unformatted binary halo
    # catalogue files output by Peak Patch. See that method for more information on what the inputs
    # for this function are.

    # Get halo data
    x,y,z,dx,dy,dz,M,redshift,zform,zcos = add_halos_pksc_to_npz( fin, mass_cutoff , rho_m ,
            z_range , z_chi_tab )

    # Save the halos in the NumPy format .npz
    np.savez( fout, cosmo_header=cosmo_header, x=x, y=y, z=z, vx=dx, vy=dy, vz=dz, M=M,
              zhalo=redshift, zform=zform, zcos=zcos )



# Add halos (must use a function not self.add_halos so it can be parallelised)
# To do: add parallelisation to self.add_halos() and add pksc2npz directly into it)
def add_halos_pksc_to_npz( catalogue_file , mass_cutoff , rho_m , z_range , z_chi_tab ):

    # Open halo catalogue file, load header
    c_in     = open( catalogue_file , 'rb' )
    N_halos  = np.fromfile( c_in, dtype=np.int32  , count=1 )[0]
    R_th_max = np.fromfile( c_in, dtype=np.float32, count=1 )[0]
    z_obs    = np.fromfile( c_in, dtype=np.float32, count=1 )[0]

    # Load halos
    cols     = int( (os.path.getsize(catalogue_file)-12)/4/N_halos )
    c        = np.reshape( np.fromfile( c_in, dtype=np.float32,
                           count=N_halos*cols ), (N_halos,cols) )

    # Perform mass cut
    c = c[ np.where( c[:,6] > mass_cutoff )[0] , : ]

    # Save relevant data from halo catalogue
    x ,y ,z  = np.ndarray.view( c[ : ,  :3 ].T )
    dx,dy,dz = np.ndarray.view( c[ : , 3:6 ].T )
    M        = np.ndarray.view( c[ : ,  6  ]   )**3 * 4./.3*np.pi * rho_m
    zform    = np.ndarray.view( c[ : , 10  ]   ) * 0.0
    # Note we're not currrently outputing the formation redshift, but it isn't needed so it
    # is set to zero
    chi      = np.sqrt( x**2 + y**2 + z**2 )
    chi_inv  = np.where( chi!=0, 1/np.where( chi==0, np.nan, chi ), 0 )
    vpec     = ( x*dx + y*dy + z*dz ) * chi_inv
    zcos     = z_of_r_comoving( chi, z_chi_table=z_chi_tab )
    # The halo redshift is given by the sum of the redshift due to the Hubble flow and the
    # redshift due to the halo's peculiar velocity, e.g. https://arxiv.org/abs/1405.0105 
    ckms     = 299792.458 # km/s
    redshift = (1+zcos) * (1+vpec/ckms) - 1

    # if no redshift range is specified, return all halos
    if z_range is None:
        return x, y, z, dx, dy, dz, M, redshift, zform, zcos

    # otherwise, perform redshift cut
    else:
        where = np.where( ( redshift > z_range[0] ) & ( redshift < z_range[1] ) )[0]
        return ( x[where] , y[where] , z[where] , dx[where] , dy[where] , dz[where] , M[where] ,
                 redshift[where] , zform[where] , zcos[where] )



####################################################################################################
###   Cosmology scripts   ##########################################################################
####################################################################################################

# Convenience functions for converting between scale factor and redshift
def a_of_z(z): return 1/(1+z)
def z_of_a(a): return 1/a-1



# Comoving distance as a function of redshift
def r_comoving_of_z( z, Omega_r0=0.0, Omega_m0=0.3138, Omega_k0=0.0, Omega_Lambda=0.6862, H_0=67.36,
        c=299792458.0 ):
    # This function numerically integrates the 1st Friedmann equation to calculate the comoving
    # radial distance to a cosmological observable given its redshift. The default values from all
    # cosmological parameters are taken from the Planck 2018 results. The default value for the
    # speed of light is the NIST definition.
    # 
    # The first Friedmann equation can be expressed in the form
    #     H^2 / H_0^2 = \Omega_{r,0} a^{-4} + \Omega_{m,0} a^{-3} + \Omega_{k,0} a^{-2}
    #                    + \Omega_\Lambda
    # This can be rearranged to get da/dt. Next, a light-like geodesic has a zero spacetime
    # interval, ds=0, so the metric for light moving in the x direction is
    #     0 = c^2 dt^2 - a^2(t) dx^2.
    # This gives us dx/dt. Integrating over dx, we can then use dx/dt and da/dt to get an integral
    # over da for x. This integrand is given by the function dchi_by_da below.
    # 
    # Arguments
    # 
    #     z : float
    #         The redshift of the cosmological observable.
    # 
    #     Omega_r0 : float, default 0.0
    #         The fraction of the total energy density that is relativistic at the time that it is
    #         observed (e.g. radiation).
    # 
    #     Omega_m0 : float, default 0.3138
    #         The fraction of the total energy density that clusters under gravity at the time it is
    #         observed (e.g. matter + dark matter).
    # 
    #     Omega_k0 : float, default 0.0
    #         The effective energy density fraction due to spatial curvature at present.
    # 
    #     Omega_Lambda : float, default 0.6862
    #         The effective energy density fraction due to cosmological constant (or vacuum energy).
    # 
    #     H_0 : float, default 67.36
    #         The Hubble constant at present in km/s/Mpc.
    # 
    #     c : float, default 299792458.0
    #         The speed of light in m/s.
    #
    # Returns
    # 
    #     chi : float or np.ndarray (same type as z)
    # 
    # Examples
    # >>> import numpy as np
    # >>> import peakpatchtools as pkp
    # >>> z = np.logspace( np.log10(1e-3) , np.log10(1e3) , 1000 )
    # >>> r_comoving_of_z( z )

    # Import the basic numerical integration function from scipy
    from scipy.integrate import quad

    # The integrand of our numerical integration
    def dchi_by_da( a, Omega_r0, Omega_m0, Omega_k0, Omega_Lambda, H_0, c_km_s ):
        return c_km_s / H_0 / np.sqrt( Omega_r0 + Omega_m0*a + Omega_k0*a**2 + Omega_Lambda*a**4 )

    # Speed of light in km/s
    c_km_s = c * 1e-3

    # Check that z is a NumPy array, if it isn't, make it into one
    type_z = None
    if not isinstance(z,np.ndarray):
        type_z = type(z)

        # Python can be dumb about integers vs floats so we make sure here that the output type is
        # not an integer
        if 'int' in repr(type_z):
            type_z = float

        # Make z a NumPy array
        z = np.array([z])

    # Perform the integral and return the comoving distance in Mpc
    chi = np.zeros( z.shape )
    for j in range(len(z)):
        chi[j] = quad( dchi_by_da, a=a_of_z(z[j]), b=1,
                       args=(Omega_r0,Omega_m0,Omega_k0,Omega_Lambda,H_0,c_km_s) )[0]

    # Return comoving distance in the same form as redshift
    if type_z is None or type_z is list:
        return chi
    else:
        return type_z(chi)



# Interpollation table for z(r_comoving) calculations
def z_of_r_comoving_table( zmin=1.0e-5, zmax=339.3, dlog10z=None, N=1000, Omega_r0=0.0, Omega_m0=0.3138,
        Omega_k0=0.0, Omega_Lambda=0.6862, H_0=67.36, c=299792458.0 ):
    # Generate a lookup table for redshift as a function of comoving distance
    # 
    # Arguments
    # 
    #     zmin : float, default 1.0e-5
    #         Minimum redshift of the lookup table. To keep things simple, we want to set this so
    #         that Peak Patch halos are going to be at redshifts greater than zmin. Typical Peak
    #         Patch runs have resolutions of hundreds of kpc to a few Mpc. The redshift of an object
    #         moving with no peculiar velocity relative to the Hubble flow at 100 kpc is about
    #         z = 2.25e-5 (or z=2.25e-4 for a distance of 1 Mpc), so to be safe, we choose
    #         zmin = 1.125e-5 to allow for comoving distances as low as 50 kpc. Typical peculiar
    #         velocities are on the order of tens to hundreds of km/s (for instance, Andromeda is
    #         moving toward the Milky Way at about 300 km/s or about 1e-3 c). Similarly, Tully-
    #         Fisher relations show that the asymptotic rotational velocity of most galaxies is less
    #         than 100 km/s, so this scale is reasonable. Redshift due to peculiar velocity is
    #             1+z = \gamma(1+v_p/c),
    #         so for small peculiar velocity relative to the speed of light v_p << c
    #             z ~ v_p/c
    #         For Andromeda-like velocities, that means z ~ 1e-3. This means that the Hubble flow
    #         redshift gets blurred out so that it is only accurate to about +/- 1e-3 or so, so
    #         there's not much point going mroe than an order of magnitude or so less than this for
    #         our minimum redshift.
    # 
    #     zmax : float, default 339.3
    #         Maximum redshift of the lookup table. Obviously, you want this to be higher than the
    #         highest value you will be interpollating. Since the Peak Patch equations of motion
    #         assume the radiation fraction is 0, we could then choose this value to be some well
    #         informed choice of a redshift above which this assumption breaks down. The matter
    #         radiation equality occurs at a redshift of z_eq=3403 +/- 26 (see 1807.06209 Table 2).
    #         Since 
    #             \Omega_m(z) = \Omega_{m,0}(1+z)^3
    #         and
    #             \Omega_r(z)=\Omega_{r,0}(1+z)^4,
    #         we can therefore set some threashold, I've chosen where the radiation density fraction
    #         \Omega_r(zmax) is 10% of the matter fraction \Omega_m(zmax) and solve for z. This
    #         gives zmax=339.3. If you wanted to be more conservative, you could go with 5% 
    #         (zmax=169.15) or 1% (zmax=33.03). First stars have been predicted to form at redshifts
    #         between 30 and 70 or so, so for most astrophysical observables, this is not an issue
    #         given Peak Patch cataogues are only accurate to the percent level when compared with
    #         N-body codes.
    #
    #     dlog10z : float, default None
    #         Log difference between redshift values in the table. The default value is None, in
    #         which case a value is chosen so that the table will have 1000 values in it.
    #             N * dlog10z = log(zmax) - log(zmin)
    #         where log() is the logarithm with base 10, so therefore
    #             dlog10z = N^-1 log(zmax/zmin).
    #         Setting N=1000, the default parmeters give us dlog10z = 7.5306e-3.
    # 
    #     N : integer, default 1000
    #         The number of values in the lookup table. This is overrided if dlog10z is set
    #         explicitly.
    # 
    #     Omega_r0 : float, default 0.0
    #         The fraction of the total energy density that is relativistic at the time that it is
    #         observed (e.g. radiation). In Peak Patch this is set to 0, so we do the same. In fact
    #         if you consider the matter radiation equality redshift z_eq = 3402 (see table 2 of
    #         1807.06209), this tells us that the radiation density fraction at present is about
    #         (9.22 +/- 0.14) * 10^{-5}. As discussed in above in the `zmax` section, setting this
    #         to 0 is find for redshifts less than a hundred or so.
    # 
    #     Omega_m0 : float, default 0.3138
    #         The fraction of the total energy density that clusters under gravity at the time it is
    #         observed (e.g. matter + dark matter).
    # 
    #     Omega_k0 : float, default 0.0
    #         The effective energy density fraction due to spatial curvature at present.
    # 
    #     Omega_Lambda : float, default 0.6862
    #         The effective energy density fraction due to cosmological constant (or vacuum energy).
    # 
    #     H_0 : float, default 67.36
    #         The Hubble constant at present in km/s/Mpc.
    # 
    #     c : float, default 299792458.0
    #         The speed of light in m/s.
    # 
    # Returns
    # 
    #     Two-dimensional array tabulating redshift as a function of comoving distance. The table
    #     indexes as follows:
    #         redshift values:           output[0,:]
    #         comoving distance values:  output[1,:]

    # If dlog10z is set to None, calculate it automatically from N, zmax and zmin
    if dlog10z is None:
        dlog10z = np.log10(zmax/zmin) / N

    else:
        # Number of values in lookup table
        N = np.ceil( np.log10(zmax/zmin) / dlog10z ).astype(np.int32)

    # Redshift values in the lookup table
    z_table = np.logspace( np.log10(zmin) , np.log10(zmax) , N )

    # Comoving distance table
    chi_table = r_comoving_of_z( z_table, Omega_r0=Omega_r0, Omega_m0=Omega_m0, Omega_k0=Omega_k0,
            Omega_Lambda=Omega_Lambda, H_0=H_0, c=c )

    # Return z and chi as one 2xN array
    return np.vstack(( z_table , chi_table ))



# Redshift as a fuction of comoving distance
def z_of_r_comoving( chi, z_chi_table=None, **kwargs ):
    # This function interpolates the redshift of an object moving with the Hubble flow as a function
    # of its comoving radial distance to the observer by calculating an interpollation table
    # generated by z_of_r_comoving_table().
    # 
    # Arguments
    # 
    #     chi : float or np.ndarray
    #         The comoving distance to a cosmological observable moving with the Hubble flow.
    # 
    #     z_chi_table : np.ndarray or None, default None
    #         An interpolation table relating redshift z and comoving distance chi. This can be
    #         generated by calling z_of_r_comoving_table().
    # 
    # Key Word Arguments
    # 
    #     **kwargs are only used if z_chi_table==None, in which case they are interpreted as
    #     arguments for the function z_of_r_comoving_table().
    #
    # Returns
    # 
    #     z : float or np.ndarray (same type as chi)
    # 
    # Examples
    # 
    # >>> import numpy as np
    # >>> import peakpatchtools as pkp
    # >>> chi = np.linspace( 0, 10000.0, 1000 )
    # >>> table = z_of_r_comoving_table()
    # >>> z_of_r_comoving( chi, z_chi_table=table )
    
    # If no table is passed, automatically make one. Key word arguments can be used to set the
    # values of the interpollation table
    if z_chi_table is None:
        z_chi_table = z_of_r_comoving_table( **kwargs )

    # Check that z is a NumPy array, if it isn't, make it into one
    type_chi = None
    if not isinstance( chi, np.ndarray ):
        type_chi = type(chi)

        # Python can be dumb about integers vs floats so we make sure here that the output type is
        # not an integer
        if 'int' in repr(type_chi):
            type_chi = float

        # Make chi a NumPy array
        chi = np.array( [ chi ] )

    # Optoinal key word argument for setting which interpollator is used
    if not hasattr(kwargs,'method'):
        kwargs['method'] = 'fast'

    # Fast logarithmic interpollator using scipy, this method has been tested against NASA
    # Extragalactic Database (NED) cosmology calculators and has less than 1% errors for redshifts
    # greater than 1e-1. Errors increase at lower redshift, but I think that is my code out-
    # performing NASA's codes because they also saturate to comoving distance of 0 below z ~ 2e-3.
    # NED calculators: https://ned.ipac.caltech.edu/Documents/References/CosmologyCalculators
    if kwargs['method'] == 'fast':
        from scipy.interpolate import interp1d

        # To avoid divide by zero errors, any value of chi less than the minimum of the
        # interpollation table is set to 10 times the table minimum, these will later be equated to
        # an interpollated redshift z=0
        chi[ chi == 0 ] = z_chi_table[1,0] / 10

        # Linearly interpolate log(z) from log(chi) using our chi(z) table, any points outside the
        # table bounds are set the minimum redshift of the table divided by 10, in the next step,
        # these are set to zero (but not here to avoid division by zero errors)
        z = 10**( interp1d( np.log10(z_chi_table[1,:]) , np.log10(z_chi_table[0,:]) ,
                bounds_error=False , fill_value=np.log10(z_chi_table[0,0]/10) )( np.log10(chi) ) )

        # Set any values outside the interpollation bounds to z=0
        z[ z < z_chi_table[0,0] ] = 0.0

    # Slow as shit interpollator that also extrapollates as power law if the comoving distance
    # values lie outside the bounds of the interpollation table
    else: # if kwargs['method'] == 'safe'

        # For values z lower than the minimum value of the table z_chi_table, we interpolate as a linear
        # function between the minimum value and z,chi=0,0 with the following slope
        m_low = z_chi_table[1,0]/z_chi_table[0,0]

        # For values z greater than the maximum value of the table z_chi_table, we extrapolate as a
        # power law so that log chi = m log z + b where
        m_high = ( np.log10( z_chi_table[1,-1] / z_chi_table[1,-2] ) / 
                   np.log10( z_chi_table[0,-1] / z_chi_table[0,-2] )   )

        b_high = ( ( np.log10( z_chi_table[0,-1] ) * np.log10( z_chi_table[1,-2] ) -
                     np.log10( z_chi_table[0,-2] ) * np.log10( z_chi_table[1,-1] )   ) /
                   np.log10( z_chi_table[0,-1] / z_chi_table[0,-2] )                     )

        # Find chi for each value of z
        z = np.zeros( chi.shape )
        outofboundsflag = False
        for j in range(len(z)):

            # Interpolate values below the minimum z in z_chi_table as a linear function from z,chi=0,0
            # to z,chi=z_chi_table[0,0],z_chi_table[0,1]
            if chi[j] < z_chi_table[1,0]:
                z[j] = chi[j] / m_low

            # Extrapolate values greater than the maximum z in z_chi_table as a power law
            elif chi[j] >= z_chi_table[1,-1]:
                z[j] = 10**( ( np.log10( chi[j] ) - b_high ) / m_high )

                # Throw a flag to warn if extrapolating beyond the table
                if outofboundsflag == False:
                    warnings.warn( ('comoving distance {0} out of bounds for table with '+
                                    '{1} < chi < {2}. Extrapolating as power law.').format( chi[j],
                                    z_chi_table[1,0] , z_chi_table[1,-1] ) , Warning )
                    outofboundsflag = True

            # Interpolate as a power law between values in the table
            else:

                for j_high in range(1,len(z_chi_table[0,:])):
                    if z_chi_table[1,j_high] > chi[j] and z_chi_table[1,j_high-1] <= chi[j]:
                        break
                j_low = j_high-1
    
                #j_low,j_high = j_low-1,j_high-1
    
                # Slope and intercept interpollation
                m_interp = ( np.log10( z_chi_table[1,j_high] / z_chi_table[1,j_low] ) /
                             np.log10( z_chi_table[0,j_high] / z_chi_table[0,j_low] )   )

                b_interp = ( ( np.log10(z_chi_table[0,j_high]) * np.log10(z_chi_table[1,j_low ]) -
                               np.log10(z_chi_table[0,j_low ]) * np.log10(z_chi_table[1,j_high])   ) /
                             np.log10( z_chi_table[0,j_high] / z_chi_table[0,j_low] )                   )

                # Interpolate chi
                z[j] = 10**( ( np.log10( chi[j] ) - b_interp ) / m_interp )

    # Return z in the same type as chi
    #if type_chi is None or type_chi is list:
    #    return z
    #else:
    #    return type_chi(z)

    return z


# Set the cosmology used in redshift-to-comoving distance calculations (next two functions
def setup_cosmology( h, N_eff, Omega_b, Omega_dm, Omega_Lambda, m_nu ):

    # Set Omega_Lambda such that Omega_total=1 if it is set to None
    if Omega_Lambda == None: 
        Omega_Lambda = 1 - Omega_b - Omega_dm

    # Load LambdaCDM cosmology class from astropy
    try:
        from astropy.cosmology import LambdaCDM
    except:
        raise ImportError( justify_text('Could not import the class `LambdaCDM\' from the module `'+
                +'astropy.cosmology\', which is necessary for the functioning of this process. Abo'+
                +'rting.'),first_line_buffer='ImportError' )
    
    # Return the LambdaCDM cosmology object
    return LambdaCDM( H0=h*100, Neff=N_eff, m_nu=m_nu, 
                      Om0=Omega_b+Omega_dm, Ob0=Omega_b, Ode0=Omega_Lambda )



# Redshift as a function of comoving distance using astropy
def z_of_r_comoving_astropy( r, z_max=3402, h=0.6736, N_eff=3.046, m_nu=0.0, Omega_b=0.0493,
        Omega_dm=0.2645, Omega_Lambda=None ):
    # Calculate the redshift as a function of comoving distance given a set of cosmological
    # parameters.
    # 
    # Arguments
    #
    #     r : float
    #         A comoving distance in Mpc at which to calculate reshift.
    # 
    #     z_max : float
    #         The maximum value of redshift for the calculation of z. The default is the Planck 2018
    #         result for the redshift of matter-radiation equality.
    # 
    # Optional arguments are:
    #     
    #     h : float
    #         The dimensionless Hubble constant $h = H_0 / (100 km/s/Mpc)$, the default value is
    #         from the Planck 2018 results.
    # 
    #     N_eff : float
    #         The effective extra relativistic degrees of freedom, (this is often referred to as
    #         effective number of neutrino species becaues if there were no realtivistic particles
    #         besides photons and neutrinos, N_eff would be the number of neutrino species) 
    #         $N_\mathrm{eff} = 8/7 (11/4)^{4/3} (\rho_r-\rho_\gamma)/\rho_\gamma$
    #         The default value is the standard model prediction of 3.046, Planck 2018 found
    #         N_eff = 2.99±0.17 so that would be another option.
    #
    #     m_nu : float
    #         The sum of the neutrino masses, we approximate this to 0 in most cases, so it is the
    #         default value.
    # 
    #     Omega_b : float
    #         Density fraction of baryonic matter, the default value is from the Planck 2018
    #         results.
    # 
    #     Omega_dm : float
    #         Density fraction of dark matter, the default value is from the Planck 2018 results.
    # 
    #     Omega_Lambda : float
    #         Density fraction of cosmological constant dark energy. If Omega_Lamda is set to None,
    #         it will be determined such that Omega_total = 1. The default for Omega_Lambda=None.
    
    # Try to load the z_at_value function from astropy
    try:
        machine = subprocess.check_output('hostname').decode('utf-8')[:-1]
        if machine[:3]+machine[-13:] == 'nia.scinet.local':
            from astropy.config import paths
            import tempfile
            # Create a temporary directory to store config
            temp_config_dir = tempfile.mkdtemp()
            paths.set_temp_config(temp_config_dir)
        from astropy.cosmology import z_at_value
    except:
        raise ImportError('Could not import function `z_at_value\' from astropy. Aborting.')

    # Set the csomology using astropy
    cosmology = setup_cosmology( h, N_eff, m_nu, Omega_b, Omega_dm, Omega_Lambda )
   
    # If r is a numpy array
    if type(r) == np.ndarray:
        z = np.zeros( len(r) )
        for j in range(len(z)):
            z[j] = z_at_value( cosmology.comoving_distance, r[j]*u.Mpc, zmax=z_max )

        # return the redshift for the specified comoving distance r in Mpc
        return z

    # If r is a number
    else:

        # return the redshift for the specified comoving distance r in Mpc
        return z_at_value( cosmology.comoving_distance, r*u.Mpc, zmax=z_max )



# Comoving distance as a function of redshift using astropy
def r_comoving_of_z_astropy( z, h=0.6736, N_eff=3.046, m_nu=0.0, Omega_b=0.0493, Omega_dm=0.2645,
        Omega_Lambda=None ):
    # Calculate the comoving distance in Mpc as a function of redshift given a set of cosmological
    # parameters.
    # 
    # Arguments
    #
    #     z : float
    #         Redshift at which to calculate the comoving distance in Mpc.
    # 
    # Optional arguments are:
    #     
    #     h : float
    #         The dimensionless Hubble constant $h = H_0 / (100 km/s/Mpc)$, the default value is
    #         from the Planck 2018 results.
    # 
    #     N_eff : float
    #         The effective extra relativistic degrees of freedom, (this is often referred to as
    #         effective number of neutrino species becaues if there were no realtivistic particles
    #         besides photons and neutrinos, N_eff would be the number of neutrino species) 
    #         $N_\mathrm{eff} = 8/7 (11/4)^{4/3} (\rho_r-\rho_\gamma)/\rho_\gamma$
    #         The default value is the standard model prediction of 3.046, Planck 2018 found
    #         N_eff = 2.99±0.17 so that would be another option.
    #
    #     m_nu : float
    #         The sum of the neutrino masses, we approximate this to 0 in most cases, so it is the
    #         default value.
    # 
    #     Omega_b : float
    #         Density fraction of baryonic matter, the default value is from the Planck 2018
    #         results.
    # 
    #     Omega_dm : float
    #         Density fraction of dark matter, the default value is from the Planck 2018 results.
    # 
    #     Omega_Lambda : float
    #         Density fraction of cosmological constant dark energy. If one of Omega_Lambda is set
    #         to None, it will be determined such that Omega_total = 1. The default for
    #         Omega_Lambda=None.
    
    # Try to load the z_at_value function from astropy
    try:
        from astropy.cosmology import z_at_value
    except:
        raise ImportError('Could not import function `z_at_value\' from astropy. Aborting.')

    # Set the csomology using astropy
    cosmology = setup_cosmology( h, N_eff, m_nu, Omega_b, Omega_dm, Omega_Lambda )
    
    # return the redshift for the specified comoving distance r in Mpc
    return cosmology.comoving_distance(z).value



# Zeldovich Linear growth factor as a fucntion of scale factor
def D_of_a( a, Omega_m=0.3138, Omega_k=0.0, Omega_Lambda=0.6862 ):
    # The Zeldovich approximation can be accurately employed in phsyical cosmology when in the
    # relatively early universe when photon pressure balances out gravitational collapse
    # sufficiently that dynamics can be modelled with linear theory in this regime, we can calculate
    # the globally averaged linear growth factor D as a function of the globally averaged FRW-metric
    # scale factor a. The growth factor D quantifies the amount that linear perturbtions from the
    # uniform fluid have grown since the end of inflation.
    # 
    # A derivation for the integral form for D(a) used below is given in sections 10-13 of "The
    # Large Scale Structure of the Universe" by P.J.E. Peebles (1980). A more succinct version can 
    # be found in the paper arXiv:astro-ph/0006089.
    # 
    # Arguments
    # 
    #     a : NumPy array or number
    #         The scale factor for the expansion of the universe with a=1 at present.
    # 
    #     Omega_m : float
    #         Density fraction of energy that clusters under gravity (matter and dark matter).
    # 
    #     Omega_k : float
    #         Effective density fractio of curvature.
    # 
    #     Omega_Lambda : float
    #         Density fraction of cosmological constant-like energy.
    # 
    # Because of the linear theory approximation, we must assume \Omega_r = 0.
    
    # Import scipy's numerical integration library
    import scipy.integrate as integrate

    # The Hubble parameter
    def H_over_H0(a, Omega_m, Omega_r, Omega_k, Omega_Lambda ):
        return ( Omega_m*a**-3 + Omega_r*a**-4 + Omega_k*a**-2 + Omega_Lambda )**(1/2)

    # The integrand
    def integrand(a, Omega_m, Omega_k, Omega_Lambda ):
        return 5*Omega_m/2 * a**-3 * H_over_H0(a,Omega_m, 0.0, Omega_k, Omega_Lambda)**-3

    # Determine D(a) using numerical integration
    D = a * integrate.quad( func=integrand, a=0.0, b=1.0, args=(Omega_m,Omega_k,Omega_Lambda) )[0]

    # Return D(a)
    return D



# Conveneience function for calcualting Zeldovich linear growth factor as a function of redshift
def D_of_z( z, Omega_m=0.3138, Omega_k=0.0, Omega_Lambda=0.6862 ):
    return D_of_a( 1/(1+z) , Omega_m, Omega_k, Omega_Lambda )



# HEC table maximum overdensity parameter as a function of redshift
def HEC_Frho_of_z( z, Omega_m=0.3138, Omega_k=0.0, Omega_Lambda=0.6862 ):
    # The Peak Patch calculation is built on the theory of homogeneous ellipsoidal collapse (HEC).
    # Overdense regions are modelled as ellipsoids of uniform density collapsing in a uniform
    # background under their own gravity. These ellipsoidal regions have no density gradient, which
    # means that the equations of motion can actually be fully solved. We use the length of time it
    # takes for such regions to collapse to estimate which regions of the initial linear density
    # field will eventually virialise within the duration of our simulations and thus seed the
    # growth of structures that would be observable. Ellipsoidal collapse can be determined as a
    # function of linear growth factor D, and are tabulated in interpollation tables prior to the
    # Peak Patch calculaiton so that the integral doesn't need to be solved for every halo.
    #
    # HEC interpollation tables are 4-dimensional and tabulate the necessary overdensity at each
    # redshift in a range for an ellipsoid of a given shape (quantified in two parameters
    # ellipticity and prolateness) to collapse by the time we observe it.
    # 
    # This function calculates approximately the largest overdensity required in these look-up
    # tables.
    # 
    # Arguments
    # 
    #     z : float
    #         The maximum redshift of the Peak Patch simulation.
    # 
    #     Omega_m : float
    #         Density fraction of energy that clusters under gravity (matter and dark matter).
    #
    #     Omega_k : float
    #         Effective density fractio of curvature.
    # 
    #     Omega_Lambda : float
    #         Density fraction of cosmological constant-like energy.

    # At z=0, this is just 1.686, so it is scaled by the inverse of the linear growth factor
    return 1.686 / D_of_z( z, Omega_m, Omega_k, Omega_Lambda )



# Peak Patch halo mass as a function of filter scale
def Mhalo_of_Rth( Rth, **kwargs ):
    # Convenience function for converting from the Peak Patch top hat lagrangian filter radius in
    # Mpc at which a halo is found and the resultant mass of the halo in solar masses. The halo mass
    # is taken to be equal to the matter density $\rho_{m,0}$ times the Lagrangian space volume
    # defined by a sphere of radius the filter scale $R_{th}$
    #     $ M_{halo} = \rho_{m,0} * 4/3 pi R_{th}^3 $
    # By default the calculation assumes Planck 2018 cosmology but this can be overridden using
    # keyword arguments outlined below.
    #
    # To calculate with a specified Hubble constant $H_0$ (km / s / Mpc) or reduced Hubble constnat
    # $h$, use keyword arguments `H_0' or `h'. E.g.:
    # >>> Mhalo_of_Rth( 1.5, H_0=70. )
    # >>> Mhalo_of_Rth( 1.5, h=0.700 )
    #
    # To calculate with specified matter density fraction $\Omega_{m,0}$, use keyword argument
    # `Omega_m'. E.g.:
    # >>> Mhalo_of_Rth( 1.5, Omega_m=0.25 )

    # Set $\rho_{m,0}$ if specified in keyword arguemnts
    if 'rho_m' in kwargs:
        rho_m = kwargs['rho_m']

    # If $\rho_{m,0}$ isn't specified, we need to check if other values upon which it depends (the
    # Hubble constant $H_0$ or reduced form $h$, or the matter density fraction $\Omega_{m,0}$) are
    # specified in keyword arguments
    else:
        # Set Hubble constant if in keyword arguments, or otherwise use Planck 2018 value
        if   'h'   in kwargs: H_0 = 100 * kwargs['h'] * u.km / u.s / u.Mpc
        elif 'H_0' in kwargs: H_0 =       kwargs['h'] * u.km / u.s / u.Mpc
        else                : H_0 = 67.36             * u.km / u.s / u.Mpc

        # Set matter density fraction if specified in keyword arguments, or use Planck 2018 value
        if 'Omega_m' in kwargs: Omega_m = kwargs['Omega_m']
        else                  : Omega_m = 0.3138

        # Gravitational constant G in units of seconds, solar masses and Mpc 
        G = const.G.to( u.Mpc**3 * u.M_sun**-1 * u.s**-2 )
        
        # Critical density in Msol / Mpc^3
        rho_crit = ( 3*H_0**2 / (8*np.pi*G) ).to( u.M_sun * u.Mpc**-3 )
            
        # mass density (without scipy units, in units of Msol / Mpc^3
        rho_m = rho_crit.value * ( Omega_m )
        
    return rho_m * 4/3 * np.pi * Rth**3



def Rth_of_Mhalo( Mhalo, **kwargs ):
    # Convenience function for converting to the Peak Patch top hat lagrangian filter radius in Mpc
    # at which a halo is found from the resultant mass of the halo in solar masses. The halo mass
    # is taken to be equal to the matter density $\rho_{m,0}$ times the Lagrangian space volume
    # defined by a sphere of radius the filter scale $R_{th}$
    #     $ M_{halo} = \rho_{m,0} * 4/3 pi R_{th}^3 $
    # By default the calculation assumes Planck 2018 cosmology but this can be overridden using
    # keyword arguments outlined below.
    #
    # To calculate with a specified Hubble constant $H_0$ (km / s / Mpc) or reduced Hubble constnat
    # $h$, use keyword arguments `H_0' or `h'. E.g.:
    # >>> Mhalo_of_Rth( 1.5, H_0=70. )
    # >>> Mhalo_of_Rth( 1.5, h=0.700 )
    #
    # To calculate with specified matter density fraction $\Omega_{m,0}$, use keyword argument
    # `Omega_m'. E.g.:
    # >>> Mhalo_of_Rth( 1.5, Omega_m=0.25 )

    # If $\rho_{m,0}$ isn't specified, we need to check if other values upon which it depends (the
    # Hubble constant $H_0$ or reduced form $h$, or the matter density fraction $\Omega_{m,0}$) are
    # specified in keyword arguments

    # Set $\rho_{m,0}$ if specified in keyword arguemnts
    if 'rho_m' in kwargs:
        rho_m = kwargs['rho_m']

    # If $\rho_{m,0}$ isn't specified, we need to check if other values upon which it depends (the
    # Hubble constant $H_0$ or reduced form $h$, or the matter density fraction $\Omega_{m,0}$) are
    # specified in keyword arguments
    else:
        # Set Hubble constant if in keyword arguments, or otherwise use Planck 2018 value
        if   'h'   in kwargs: H_0 = 100 * kwargs['h'] * u.km / u.s / u.Mpc
        elif 'H_0' in kwargs: H_0 =       kwargs['h'] * u.km / u.s / u.Mpc
        else                : H_0 = 67.36             * u.km / u.s / u.Mpc

        # Set matter density fraction if specified in keyword arguments, or use Planck 2018 value
        if 'Omega_m' in kwargs: Omega_m = kwargs['Omega_m']
        else                  : Omega_m = 0.3138

        # Gravitational constant G in units of seconds, solar masses and Mpc 
        G = const.G.to( u.Mpc**3 * u.M_sun**-1 * u.s**-2 )
        
        # Critical density in Msol / Mpc^3
        rho_crit = ( 3*H_0**2 / (8*np.pi*G) ).to( u.M_sun * u.Mpc**-3 )
            
        # mass density (without scipy units, in units of Msol / Mpc^3
        rho_m = rho_crit.value * ( Omega_m )
    
    return ( 3 / ( 4 * np.pi * rho_m ) * Mhalo )**(1./3.)



# Maximum halo radius as a function of redshift
def R_th_max( z , R_th_max_0=34. , zj=5. , R_th_max_j=2 ):
    # For a given redshift, there is a maximum size of halo that can have formed. For redshift 0,
    # using Planck 2018 Lambda CDM cosmology, this is typically around 34 Mpc. The maximum halo
    # radius as a function of redshift is roughly an exponential decay of the form
    #
    #     R_th_max(z) = R_th_max_0 e^( - k z )
    #
    # where
    #
    #     k = - 1 / z_j ln( R_th_max_j / R_th_max_0 )
    #
    # and R_th_max_j is another value of the maximum halo radius at a redshift z_j. For default
    # values I used R_th_max_j=2 and zj=5.
    return R_th_max_0 * np.exp( 1./zj * np.log( R_th_max_j/R_th_max_0 ) * z )



####################################################################################################
###   Non-Gaussianity scripts   ####################################################################
####################################################################################################

# Fit functions relating the log ratio of P_ping - P_gauss over P_gauss as a hyperbola
def lnD_of_lnk_n( lnk, lnk_pulse, lnD_0, lnn ):
    return lnD_0 + 3*lnn - 3*lnk_pulse + 3*lnk + 3*lnn**2/(lnk-lnn-lnk_pulse)
def D_of_k_n( k, k_pulse, D_0, n ):
    return (   D_0 
             * n**3 
             * (k/k_pulse)**3 
             * np.exp( 3*np.log(n)**2 / np.log(k/k_pulse/n) )
             * np.heaviside( n*k_pulse-k, 0 )
           )



####################################################################################################
###   Math scripts   ###############################################################################
####################################################################################################

def solid_angle_subtended_by_square( s, r, units=None ):
    # Well, this one is pretty self explanatory, it calculates the solid angle subtended by a square
    # of side length s lying on the z axis at a distance r from an observer.
    # 
    # You can solve this problem analytically, it's not pretty, but it comes out to the following:
    # $$ \Omega(s,r) = 2\pi - 8\arctan( ( 1 + s^2/(2r^2))^{-1/2} ) ) $$
    # Which of course goes to half the sky (2pi sr) when r goes to 0 and goes to 0 as r/s goes to 0.
    # 
    # Arguments:
    # 
    #     s : float or astropy quantity
    #         Side length of the square.
    # 
    #     r : float or astropy quantity
    #         The distance to the square. If r or s is an astropy quantity, both are converted to
    #         the same units as s and then just absolute values are taken. If one or neither have
    #         units, they are assumed to be of the same units.
    # 
    #     units : astropy quantity
    #         The units (of solid angle) in which to output the result (the default is steradians.

    # Set the default units of units is set to None
    if units == None:
        units = u.sr

    # If s and r have units, make sure they are units that make sense
    if isinstance(s,u.quantity.Quantity):
        if not s.unit.is_equivalent(u.meter):
            raise AttributeError('s must have units of length, not '+str(s.unit)+'.')
    if isinstance(r,u.quantity.Quantity):
        if not not r.unit.is_equivalent(u.meter):
            raise AttributeError('r must have units of length, not '+str(r.unit)+'.')

    # Make sure units has units of solid angle
    if not isinstance(units,u.core.Unit) and not isinstance(units,u.core.CompositeUnit):
        raise TypeError('units must be of type astropy.units.quantity.Quantity.')
    if not (1*units).unit.is_equivalent(u.sr):
        raise AttributeError('units must have units of solid angle, not '+str(units.unit)+'.')

    # Convert u and r to dimensionless ratios
    if isinstance(s,u.quantity.Quantity) and isinstance(r,u.quantity.Quantity):
        if s.unit == r.unit:
            s, r = s.value, r.value
        else:
            s, r = s.value, r.to(s.unit).value
    elif isinstance(s,u.quantity.Quantity) and isinstance(r,float):
        s, r = s.value, r
    elif isinstance(s,float) and isinstance(r,u.quantity.Quantity):
        s, r = s, r.value
    elif not isinstance(s,float) and not isinstance(r,float):
        raise TypeError('s and r must be of type float or astropy.units.quantity.Quantity.')

    # Calculate the solid angle and return it
    return ( ( 2*np.pi - 8*np.arctan( 1/np.sqrt( 1 + s**2/(2*r**2) ) ) ) * u.sr ).to( units )



def levi_civita(dim,*args):
    # The permutation pseudo tensor (or Levi-Civita symbol) of dimension dim. This takes `dim`
    # arguments args[0], args[1], ... args[dim-2], args[dim-1] which must be integers between 0 and
    # dim-1 (inclusive) and outputs +1 if the arguments are an even permutation of
    # 0, 1, ... dims-2, dims-1 or -1 if they are an odd permutation.
    # 
    # Arguments:
    # 
    #     dim : int
    #         The dimension of the Levi-Civita symbol.
    # 
    #     args[0], args[1], ... args[dim-2] , args[dim-1] : integers
    #         dim integers representing indices of the pseudo-tensor.
    # 
    # Examples:
    # >>> import peakpatchtools as pkp
    # >>> pkp.levi_civita(3,0,1,2)
    # 1
    # >>> pkp.levi_civita(4,0,1,3,4)
    # -1
    # >>> pkp.levi_civita(3,0,0,0)
    # 0

    # Import python tensor library
    import sympy

    # Make sure the correct number (dim) of arguments is passed
    if len(args) != dim:
        raise IndexError('the function must take a dimension dim and dim arguments.')

    # Make sure arguments are integers between 0 and dim-1 (inclusive)
    for i in range(dim):
        if args[i] != args[i]//1:
            raise TypeError('args should be integers.')
    if np.max(np.array(args)) >= dim or np.min(np.array(args)) < 0:
        raise ValueError('args must be between 0 and 1-dim (inclusive).')

    # Return +1 for even permutations, -1 for odd permutations or 0 otherwise
    return int( 
            sympy.prod(
                    sympy.prod( args[j]-args[i] for j in range(i+1,dim) ) 
                            / np.math.factorial(i) for i in range(dim)
                      )
              )



def fourier_top_hat_filter_2D(data, cutoff_frequency):
    # Apply a top-hat filter to a 2D array in Fourier space.
    # Parameters:
    # data (2D array): The input data array to be smoothed.
    # cutoff_frequency (float): The cutoff frequency for the top-hat filter (in cycles per unit distance).
    # 
    # Returns:
    # 2D array: The smoothed data.
    
    # Get the shape of the data
    nx, ny = data.shape

    # Compute the Fourier transform of the data
    data_ft = fft2(data)
    data_ft_shifted = fftshift(data_ft)

    # Create the frequency grid
    kx = np.fft.fftfreq(nx).reshape(-1, 1) * nx
    ky = np.fft.fftfreq(ny).reshape(1, -1) * ny
    k  = np.sqrt(kx**2 + ky**2)

    # Define the top-hat filter
    filter_mask = k <= cutoff_frequency

    # Apply the filter
    filtered_data_ft_shifted = data_ft_shifted * filter_mask

    # Inverse Fourier transform to go back to real space
    filtered_data_ft = ifftshift(filtered_data_ft_shifted)
    smoothed_data = np.real(ifft2(filtered_data_ft))

    return smoothed_data



# Log-log interpollator
def loglog_interp( x_interp, x, y, verbose=False ):

    # If x_interp is less than the minimum value of x, extrapolate as an exponential
    if x_interp < x[0]:

        # Throw a worning
        if verbose == True:
            print('Warning, x_interp out of range, extrapolating as log.')

        # Extrapolate as an exponential
        y_interp = y[0] * np.exp(   np.log( y[1]     / y[0] )
                                  / np.log( x[1]     / x[0] )
                                  * np.log( x_interp / x[0] ) )

    # If x_interp is greater than the maximum value of x, extrapolate as an exponential
    elif x_interp >= x[-1]:

        # Throw a warning
        if verbose == True:
            print('Warning, x_interp out of range, extrapolating as log.')

        # Extrapolate as an exponential
        y_interp = y[-2] * np.exp(   np.log( y[-1]    / y[-2] )
                                   / np.log( x[-1]    / x[-2] )
                                   * np.log( x_interp / x[-2] ) )

    else:

        for j in range(len(x)-1):
            if x[j] <= x_interp and x[j+1] > x_interp:
                j_interp = j
                break
        y_interp = y[j_interp] * np.exp(   np.log( y[j_interp+1]   / y[j_interp] )
                                         / np.log( x[j_interp+1]   / x[j_interp] )
                                         * np.log( x_interp        / x[j_interp] ) )
    return y_interp



####################################################################################################
###   Image manipulation: scripts for image handling   #############################################
####################################################################################################

def animator( images, out_fig, figsize=[], dpi=300, fps=2, **kwargs ):
    # Function for stitching together individual images into a gif.
    #
    # Arguments:
    # 
    #     images : list of str
    #         List of strings indicating the path to an image file.
    # 
    #     out_fig : str
    #         String indicating the directory and file name at which to save the animation.
    # 
    # Optional arguments:
    # 
    #     figsize : list of 2 integers
    #         The size of the images in inches [width,height]. If an empty list [] is passed, then
    #         the native image size (in pixels) is determined and figsize is set accordingly taking
    #         into account dpi. Default is [].
    # 
    #     dpi : int
    #         Resolution of the images (dots per inch) so the resultant image will be figsize[0]*dpi
    #         by figsize[1]*dpi (length*width). Default is 300.
    # 
    #     fps : int
    #         The number of frames per second for the animation. Defualt is 2.
    # 
    # Key Word Arguments:
    # 
    # There are three options for text: either a single string of text for all frames, a unique
    # string of text for each frame, or a list of unique strings of text for each frame. Accordingly
    # these can be handled using kwargs that must be lists of length 1, or the number of images. If
    # there are to be multiple strings per frame, then each of thees key word arguments must be a
    # list of length the number of images with each element a list of length the number of
    # individual strings to display in each frame of the animation.
    # 
    #     text : str, or list (or list of lists) of strings, length 1 or equal to len(images)
    #         Strings to print on each frame of the animation.
    # 
    #     text_pos : list of 2 numbers, or list (or list of lists) of lists of 2 numbers
    #         The x,y coordinates of the text strings (in pixels with (0,0) being the top left,
    #         (figsize[0]*dpi,figsize[1]*dpi) being the bottom right). Default is [0,0].
    #
    #     text_color : list (or list of lists) of strings
    #     text_bg_color : list (or list of lists) of strings
    # 
    # Examples:
    # 
    # 1) If a single string "string" is to be displayed on all frames:
    # >>> animator( ['1.png','2.png'], '3.gif', text='string', text_pos=[0,0] )
    # 
    # 2) If a the strings 'one' and 'two' are to be displayed in the middle of a two-frame gif with
    # dimensions 1200x1200 pixels:
    # >>> animator( ['1.png','2.png'], '3.gif', figsize=[4,4], dpi=300, text=['one','two'],
    # ... text_pos=[600,600] )
    # 
    # 3) combining the two cases above where the first should be font size 24 and the second 12 with
    # font colors black and white respectively:
    # >>> animator( ['1.png','2.png'], '3.gif', figsize=[4,4], dpi=300, 
    # ... text=[['string','one'],['string','two']], text_pos=[[[0,0],[600,600]],[[0,0],[600,600]]],
    # ... text_size=[[24,12],[24,12]], text_color=[['k','w'],['k','w']] )
    #
    # There are also a few key word arguments to edit the animation settings
    # 
    #     anim_writer : str
    #         The only supported writer right now is 'PillowWriter', so obviously that's the
    #         default.
    # 
    #     anim_repeat : bool
    #         Whether the gif should repeat. Default is True
    # 
    #     anim_repeat_delay : int
    #         The delay (in miliseconds) after the end of the animation before it repeats (if
    #         anim_repeat==True). The default is to set it to the same as the delay between frames,
    #         1000/fps.
    # 
    #     anim_blit : bool
    #         Whether to use bit blitting (wikipedia.org/wiki/Bit_blit), which optimizes drawing.
    #         Default is True.

    import matplotlib.image as mgimg
    from matplotlib import animation

    # Set size of images read in using imread
    if figsize == []:
        img = mgimg.imread(images[0])
        height,width,*channels = img.shape
        figsize = [width/dpi,height/dpi]
    fig,ax = plt.subplots( figsize=figsize, dpi=dpi )

    # Number of frames
    N_frames = len(images)

    # If the text key word argument is passed, ensure all text-related kwargs properly formatted
    if 'text' in kwargs:

        # If just a string is passed, set it to a list
        if type(kwargs['text']) == str: kwargs['text'] = [ kwargs['text'] ]

        # If other keyword arguemnts not present, use default values
        if 'text_pos'      not in kwargs: kwargs['text_pos'     ] = [[0,0]]
        if 'text_size'     not in kwargs: kwargs['text_size'    ] = [12]
        if 'text_color'    not in kwargs: kwargs['text_color'   ] = ['k']
        if 'text_bg_color' not in kwargs: kwargs['text_bg_color'] = ['none']

        # All text-related key word arguments
        text_keys = ('text','text_pos','text_size','text_color','text_bg_color')

        # Each text-related key word argument in the dictionary kwargs must end up formatted as a
        # list of the form kwargs[key][ frame# ][ string# ] where frame# is the number of frames in
        # the animation given by len(images) and string# is the number of text strings printed on
        # each frame.
        l = ( 1, N_frames )
        for key in text_keys:
                
            # If just one value is passed not as a list or tuple
            if key == 'text_pos' and type(kwargs[key][0]) not in (list,tuple):
                kwargs[key][0] = [ kwargs[key][0] ]
            elif type(kwargs[key]) not in (list,tuple):
                kwargs[key] = [ kwargs[key] ]

            # Check that argument lists have the allowed lengths
            if len(kwargs[key]) not in l:
                raise IndexError(f'{key} must be list of length 1 or len(images).')

        # If one of the text key word arguments is a list of length 1 (meaning this 1 value is to be
        # used for all frames), turn it into a list of length len(images) where all elements have
        # the value of the initial lone element.
        if l[1] > 1:
            for key in text_keys:
                if len(kwargs[key]) == 1: 
                    kwargs[key] = [ kwargs[key][0] for j in range(l[1]) ]

        # Determine the number of strings per frame
        if type(kwargs['text'][0]) not in (list,tuple):
            N_strings = 1
        else:
            N_strings = len( kwargs['text'][0] )

        # Even if there is only one string to be added to the figure, to make it so that we only
        # need to make one imgplot.axes.text call, even if len(images)==1, the arguments should
        # still be formatted as a list of lists.
        for key in text_keys:
            if key == 'text_pos':
                if type(kwargs[key][0][0]) not in (list,tuple):
                    for j in range(l[1]):
                        kwargs[key][j] = [ kwargs[key][j] for k in range(N_strings) ]
            else:
                if type(kwargs[key][0]) not in (list,tuple):
                    for j in range(l[1]):
                        kwargs[key][j] = [ kwargs[key][j] for k in range(N_strings) ]

    # Make each frame of the animation
    animation_frames = []
    for p in range( N_frames ):

        # Read image file and plot it as a matplotlib object
        img     = mgimg.imread( images[p] )
        imgplot = ax.imshow( img, aspect='auto' )

        # Remove axes and white space around image
        plt.tight_layout()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.spines['top'   ].set_visible(False)
        ax.spines['right' ].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.axes.spines['left'  ].set_visible(False)

        # If the text key word argument is passed, add text to image
        if 'text' in kwargs:
            text=[ 0 for j in range(N_strings) ]
            for j in range(N_strings):

                text[j] = ax.text( kwargs['text_pos'][p][j][0], kwargs['text_pos'][p][j][1],
                                   kwargs['text'][p][j], fontsize=kwargs['text_size'][p][j],
                                   color=kwargs['text_color'][p][j],
                                   backgroundcolor=kwargs['text_bg_color'][p][j] )

            # Add frame to list of animation frames (both the imshow plot, imgplot, and the text
            # must be passed so that both are updated each frame
            animation_frames.append( [imgplot]+text )

        # If there is no text, only the imshow plot imgplot needs to be updated with each frame
        else:
            animation_frames.append([imgplot])

    # Check for other kwargs related to animation
    if 'anim_writer' in kwargs: anim_writer = kwargs['anim_writer']
    else                      : anim_writer = 'PillowWriter'
    if 'anim_blit' in kwargs: blit = kwargs['anim_blit']
    else                    : blit = True
    if 'anim_repeat' in kwargs: repeat = kwargs['anim_repeat']
    else                      : repeat = True
    if 'anim_repeat_delay' in kwargs: repeat_delay = kwargs['anim_repeat_delay']
    else                            : repeat_delay = 1000/fps

    # Make animation
    anim = animation.ArtistAnimation( fig=fig, artists=animation_frames, interval=1000/fps,
            blit=blit, repeat=repeat, repeat_delay=repeat_delay )

    # Save the animation
    if anim_writer == 'PillowWriter':
        from matplotlib.animation import PillowWriter
        writer = PillowWriter( fps=fps )
        anim.save( out_fig, writer=writer, dpi=dpi )
    else:
        raise ValueError('this value of anim_writer is not supported.')



####################################################################################################
###   Time Keeping: scripts for timing functions   #################################################
####################################################################################################

# Time objects
class my_time:
    def __init__(self,t):
        # t is a time in seconds
        self.t =           t
        self.h = int(  abs(t) // 3600        )
        self.m = int( (abs(t) %  3600) // 60 )
        self.s =       abs(t) %  60

        # sign
        if t >= 0: self.sign=''
        else:      self.sign='-'

    # Types
    def __repr__(self):  return f'{self.sign}{self.h:02}:{self.m:02}:{self.s:06.3f}'
    def __str__(self):   return f'{self.sign}{self.h:02}:{self.m:02}:{self.s:06.3f}'
    def __float__(self): return float(self.t)
    def __int__(self):   return int(self.t)
   
    # Math operators
    def __add__     (self,arg): return my_time( self.t +  arg.t )
    def __sub__     (self,arg): return my_time( self.t -  arg.t )
    def __mul__     (self,arg): return my_time( self.t *  arg   )
    def __truediv__ (self,arg): return my_time( self.t /  arg   )
    def __floordiv__(self,arg): return my_time( self.t // arg   )
    def __mod__     (self,arg): return my_time( self.t %  arg   )
    def __pow__     (self,arg): return my_time( self.t ** arg   )

    # Unitary math operators
    def __neg__(self): return my_time(    -self.t  )
    def __pos__(self): return my_time(    +self.t  )
    def __abs__(self): return my_time( abs(self.t) )
   
    # Outputting just components of time
    def hr(self):        return self.t/3600
    def m(self):       return self.t/60
    def s(self):       return self.t



# Timer
class timer:

    def __init__(self):
        self.times = []

    def start(self):
        now = datetime.datetime.now()
        self.times = [ my_time( time.mktime(now.timetuple())+now.microsecond*1e-6 ) ]
        self.dates = [ now ]
        return self.times[0]

    def lap(self):
        if len(self.times) == 0:
            raise IOError('timer.start must be called before timer.lap.')
        now = datetime.datetime.now()
        self.times += [ my_time( time.mktime(now.timetuple())+now.microsecond*1e-6 ) ]
        self.dates += [ now ]
        return self.times[-1]-self.times[-2]

    def stop(self):
        now = datetime.datetime.now()
        self.times += [ my_time( time.mktime(now.timetuple())+now.microsecond*1e-6 ) ]
        self.dates += [ now ]
        return self.times[-1]-self.times[0]

    def timestamp(self):
        return self.dates[-1].date()



####################################################################################################
###   Wordplay: scripts for text handling   ########################################################
####################################################################################################

def justify_text( string, width=100, first_line_indent=0, first_line_buffer=0 ):
    # This is a script for justifying text.
    # 
    # Parameters
    # 
    #     string : str
    #         The figure in which to plot the field slice.
    # 

    # Split any hyphenated words and / words
    def split_hyphenated(str_lst,hyphen='-'):
        # Accepts str_lst, a list of strings, typically a sentence that has been split into a list 
        # using the .slit() intrinsic function. Finds any words that have a specific character in
        # them `hyphen` (with the default being a hyphen '-') and separates them into two strings,
        # where the first retains the `hyphen' character.

        # Find all instances of character `hyphen'
        j      = 0
        length = len(str_lst)
        while j < length:

            # If a word has a hyphen
            if hyphen in str_lst[j]:

                # Split the word into the two pieces on either side of each `hyphen', then add a
                # `hyphen' back onto each word preceding it
                str_sublst = str_lst[j].split(hyphen)
                offset     = len(str_sublst)-1
                for k in range(offset):
                    str_sublst[k] = str_sublst[k]+hyphen

                # Update number of distinct words in the list
                length += offset

                # Reform the full list
                if j+1 < len(str_lst):
                    str_lst = str_lst[:j] + str_sublst + str_lst[j+1:]
                else:
                    str_lst = str_lst[:j+offset] + str_sublst

                # Increment j by the number of new words added
                j += offset

            # Increment j by one to look at the next new word
            j+=1

        return str_lst

    # Set the indentation of the first line
    if type(first_line_indent) == int:
        first_line_indent = ' '*first_line_indent
    elif type(first_line_indent) != str:
        raise TypeError('first_line_indent must be of type str or int.')
    if len(first_line_indent) > width:
        raise ValueError('first_line_indent must be less than or equal to width, obviously.')

    # Add a buffer to the first line (this is useful for error messages which print the error type
    # before the message you give
    if type(first_line_buffer)==str:
        if ( ('Error' in first_line_buffer or 'Warning' in first_line_buffer) and
                ':' not in first_line_buffer and ' ' not in first_line_buffer ):
            first_line_buffer = len(first_line_buffer)+2
        else:
            first_line_buffer = len(first_line_buffer)
    elif type(first_line_buffer)!=int:
        raise TypeError('first_line_buffer must be of type str or int.')

    # Split paragraphs (separated by '\n')
    paragraph_list = string.split('\n')

    # Output
    string = ''

    # Loop over paragraphs
    for i in range(len(paragraph_list)):

        # Split string into words
        string_list = paragraph_list[i].split()

        # Split the list of strings wherever there is a hyphen or forward slash
        string_list = split_hyphenated(string_list)
        string_list = split_hyphenated(string_list,hyphen='/')

        if i == 0:
            # Make sure buffer isn't larger than the line width
            if first_line_buffer + len(first_line_indent) >= width:
                string+='\n'
                first_line_buffer=0

            # Add buffers and indentation to first line
            line_length = first_line_buffer + len(first_line_indent)
            string      = first_line_indent

        # Add to output string
        for j in range(len(string_list)):

            # For the first word in the line
            if j==0 or line_length==0:
                if len(string_list[j]) + line_length <= width:
                    string      += string_list[j]
                    line_length += len(string_list[j])

            # If the previous word had a hyphen or slash in it
            elif string[-1] in '-/' and len(string_list[j]) + line_length <= width:
                    string      += string_list[j]
                    line_length += len(string_list[j])

            # Any other words that aren't going to overflow the width limit
            elif len(string_list[j]) + line_length + 1 <= width:
                if string[-1] in '-/':
                    string      += string_list[j]
                    line_length += len(string_list[j])
                else:
                    string      += ' '+string_list[j]
                    line_length += len(string_list[j])+1

            # If a single word is longer than the width limit
            elif len(string_list[j]) > width:
                if j>0: string+='\n'
                while len(string_list[j]) > width:
                    string         += string_list[j][:width] + '\n'
                    string_list[j]  = string_list[j][width:]
                    line_length     = len(string_list[j])
    
            # If word overflows the width, start a new line
            else:
                string      += '\n' + string_list[j]
                line_length  = len(string_list[j])

        # Add newline character if the input string has a '\n' in it
        if i+1 < len(paragraph_list):
            string += '\n'

    return string



# Find all instances of a character in a string
def findall(string,character):
    return [ i for i, character_i in enumerate(string) if character_i == character ]



####################################################################################################
###   Other related scripts   ######################################################################
####################################################################################################

# Planck colour scale for matplotlib
def planck_cmap():
    # Use matplotlib colour map that approximates the Planck colour map
    # 
    # Adapted from script by Andrea Zonca, see their blog post for more information:
    # http://zonca.github.io/2013/09/Planck-CMB-map-at-high-resolution.html

    from matplotlib.colors import ListedColormap

    # Load the RGB file approximating the Planck colour map, retrieved from:
    # https://github.com/zonca/paperplots/raw/master/data/Planck_Parchment_RGB.txt
    planck_cmap_file = peak_patch_dir+'/tables/Planck_Parchment_RGB.txt'
    colombi1_cmap    = ListedColormap( np.loadtxt(planck_cmap_file) / 255.0 )

    # Set colour of missing pixels to grey
    colombi1_cmap.set_bad((0,0,0,0))

    return colombi1_cmap



# Planck colour scale for matplotlib
def planck_cmap_mollview():
    # Use matplotlib colour map that approximates the Planck colour map
    # 
    # Adapted from script by Andrea Zonca, see their blog post for more information:
    # http://zonca.github.io/2013/09/Planck-CMB-map-at-high-resolution.html

    from matplotlib.colors import ListedColormap

    # Load the RGB file approximating the Planck colour map, retrieved from:
    # https://github.com/zonca/paperplots/raw/master/data/Planck_Parchment_RGB.txt
    planck_cmap_file = peak_patch_dir+'/tables/Planck_Parchment_RGB.txt'
    colombi1_cmap    = ListedColormap( np.loadtxt(planck_cmap_file) / 255.0 )

    # Set colour of missing pixels to grey
    colombi1_cmap.set_bad('grey')

    # Set background colour to white (necessary if you want to use the colour map directly with
    # hp.mollview(m, cmap=planck_cmap())
    colombi1_cmap.set_under((0,0,0,0))
    # Note that mollview stupidly puts the background (outside the ellipse) as some value below
    # vmin, so you need this to make it not just appear the same colour as the minimum of your
    # colourbar. But then if you change vmin, you can end up with any clipped values showing up as
    # white, so don't use this map if you are adjusting vmin.

    return colombi1_cmap



# Planck colour scale for matplotlib
def earth_tones_cmap():
    # Use matplotlib colour map that approximates the colours of the Earth

    from matplotlib.colors import ListedColormap

    # Load the RGB file approximating Earth colours
    earth_tones_cmap_file = peak_patch_dir+'/tables/earth_tones.txt'
    colombi1_cmap         = ListedColormap( np.loadtxt(earth_tones_cmap_file) / 255.0 )

    # Set colour of missing pixels to grey
    colombi1_cmap.set_bad('gray')

    # Set background colour to white (necessary if you want to use the colour map directly with
    # hp.mollview(m, cmap=planck_cmap())
    colombi1_cmap.set_under('white')

    return colombi1_cmap
