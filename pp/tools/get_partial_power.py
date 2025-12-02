import os,sys

# Read in command line arguments
field_file = sys.argv[1]
nbuff      = int(sys.argv[2])
s_box      = float(sys.argv[3])

# Modules to load on CITA computers
intel_module = 'intel/19.1.3'
fftw_module  = 'fftw/3.3.10-intelmpi'

# Paths etc. can be found on CITA servers by running `module show fftw` etc.
fortrancompiler = 'ifort'
fortranoptions  = '-w -traceback' # -fpp'
fftw_include    = '-I/cita/modules/fftw/3.3.10-intelmpi/include'
mpi_include     = '-I/cita/modules/intel/19.1.3/compilers_and_libraries/linux/mkl/include'
fftw_lib_dir    = '-L/cita/modules/fftw/3.3.10-intelmpi/lib'
mpi_lib_dir     = '-L/cita/modules/intel/19.1.3/impi/2019.9.304/intel64/lib'
fftw_libs       = '-lfftw3f'# -lfftw3f_mpi -lfftw3f_omp -lfftw3f_threads -lm'

# Fortran script and executable
f90_gps = 'get_partial_power.f90'
exe_gps = 'get_partial_power'
os.system('chmod -x {0}'.format( f90_gps ))

# Run scripts
os.system(('module purge;module load {0} {1};'+
           '{2} {3} {4} {5} {6} {7} {8} {9} -o {10};'+
           './{10} << EOF\n"{11}"\n{12}\n{13}\nEOF')
    .format( intel_module,fftw_module,fortrancompiler,fortranoptions,
             mpi_include,mpi_lib_dir,fftw_include,fftw_lib_dir,fftw_libs,
             f90_gps,exe_gps,field_file,nbuff,s_box
        )
    )
