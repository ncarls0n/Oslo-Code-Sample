import os,sys

# Read in command line arguments
field_file = sys.argv[1]
nbuff      = int(sys.argv[2])
s_box      = float(sys.argv[3])

# Modules to load on CITA computers
NiaEnv_module = 'NiaEnv/2019b'
intel_module  = 'intel/2019u4'
mpi_module    = 'intelmpi/2019u4'
fftw_module   = 'fftw/3.3.8'
#gsl_module    = 'gsl/2.5'
modules = '{0} {1} {2} {3}'.format(NiaEnv_module,intel_module,mpi_module,fftw_module)

# Paths etc. can be found on CITA servers by running `module show fftw` etc.
fortrancompiler = 'ifort'
fortranoptions  = '-w -traceback' # -fpp'
fftw_include    = '-I/scinet/niagara/software/2019b/opt/intel-2019u4-intelmpi-2019u4/fftw/3.3.8/include'
mpi_include     = ''
fftw_lib_dir    = '-L/scinet/niagara/software/2019b/opt/intel-2019u4-intelmpi-2019u4/fftw/3.3.8/lib'
mpi_lib_dir     = ''#'-L/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/mpi/intel64/lib'
fftw_libs       = '-lfftw3f_mpi -lfftw3f'

# Fortran script and executable
f90_gps = 'get_powerspectrum.f90'
exe_gps = 'get_powerspectrum'
os.system('chmod -x {0}'.format( f90_gps ))

# Run scripts
os.system(('module purge;module load {0};'+
           '{1} {2} {3} {4} {5} {6} {7} {8} -o {9};'+
           './{9} << EOF\n"{10}"\n{11}\n{12}\nEOF')
    .format( modules,fortrancompiler,fortranoptions,
             mpi_include,mpi_lib_dir,fftw_include,fftw_lib_dir,fftw_libs,
             f90_gps,exe_gps,field_file,nbuff,s_box
        )
    )
