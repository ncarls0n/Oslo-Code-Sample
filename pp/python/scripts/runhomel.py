#!/usr/bin/env python
import os

# Compilers
fortrancompiler = 'gfortran'
python2         = 'python2'
python3         = 'python3'

# Copiler options (note that this string must begin and end with a space)
fortranoptions = ' -O4 -w '
# -Oi is optimization with -O1 being minimum and -O4 being maximum
# -w flag inhibits all warnings

# Directories that HomogeneousEllipsoid.f90 makes calls to
dir_python          = '.'
dir_peakpatch       = dir_python+'/..'
dir_src             = dir_peakpatch+'/src'
dir_homel           = dir_src+'/modules/HomogeneousEllipsoid'
dir_Solvers         = dir_src+'/modules/Solvers'
dir_External        = dir_src+'/modules/External'
dir_GlobalVariables = dir_src+'/modules/GlobalVariables'
dir_cosmology       = dir_src+'/cosmology'

# HomogeneousEllipsoid dependencies
f90_Solvers          = dir_Solvers+'/Solvers.f90'
f90_intreal_types    = dir_External+'/intreal_types.f90'
f90_input_parameters = dir_GlobalVariables+'/input_parameters.f90'
f90_params           = dir_GlobalVariables+'/params.f90'
f90_cosmoparams      = dir_GlobalVariables+'/cosmoparams.f90'
f90_psubs_Dlinear    = dir_cosmology+'/psubs_Dlinear.f90'
f90_homel            = dir_homel+'/HomogeneousEllipsoid.f90'

# All scripts
scripts_f90 = (  f90_Solvers         +' '+f90_intreal_types+' '
                +f90_input_parameters+' '+f90_params       +' '
                +f90_cosmoparams     +' '+f90_psubs_Dlinear+' '
                +f90_homel
              )

# Executable
exe_homel = 'run_homel'

# Compile and Link .f90 files
#os.system
print(  'cd '+dir_homel+';'
           +fortrancompiler+fortranoptions+scripts_f90+' -o '+exe_homel
         )


"""
That will compile and link all the necessary .f90 scripts.
Next you need to write a fortran script to run a homel, probably start with Dick's.

"""
