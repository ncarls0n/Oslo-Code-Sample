import os,sys

# Read in command line arguments
field_file = sys.argv[1]
nbuff      = int(sys.argv[2])
s_box      = float(sys.argv[3])

# Modules to load on CITA computers
intel_module = 'intel'
fftw_module  = 'fftw'

fortrancompiler = 'ifort'
fortranoptions  = ('-w'                  # suppress all warnings
                  +' -fpp'               # use the Fortran preprocessor
                  +' -qopenmp -parallel')# use OpenMP compiler directives for multi-threading
                  #+' -traceback'        # make debugging easier with traceback
fftw_include    = '-I/cita/modules/fftw/3.3.10/include'
fftw_lib_dir    = '-L/cita/modules/fftw/3.3.10/lib'
fftw_libs       = '-lfftw3f_omp -lfftw3f_threads -mkl'# -lfftw3f_mpi -lfftw3f -lm'
# Paths to libraries and inclusions can be found on CITA servers by running e.g. `module show fftw`

# Fortran script and executable
f90_gps = 'get_power_omp.f90'
exe_gps = 'get_power_omp'
os.system('chmod -x {0}'.format( f90_gps ))

# Run scripts
os.system(('module purge;module load {0} {1};'+
           '{2} {3} {4} {5} {6} {7} -o {8};'+
           './{8} << EOF\n"{9}"\n{10}\n{11}\nEOF')
    #        {0}            {1}           {2}               {3}              {4}
    .format( intel_module , fftw_module , fortrancompiler , fortranoptions , fftw_include ,
    #        {5}            {6}           {7}               {8}              {9}
             fftw_lib_dir , fftw_libs   , f90_gps         , exe_gps        , field_file   ,
    #        {10}           {11}
             nbuff        , s_box
        )
    )
