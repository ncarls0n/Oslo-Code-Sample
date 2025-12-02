import os
import sys

field_file = sys.argv[1]

fortrancompiler = 'gfortran'
fortranoptions  = '-w'#'-04 -w'
fftw_include    = '-I/usr/local/Cellar/fftw/3.3.10/include'
mpi_include     = ''#'-I/usr/local/Cellar/open-mpi/4.1.2/include'
fftw_lib_dir    = '-L/usr/local/Cellar/fftw/3.3.10/lib'
mpi_lib_dir     = ''#'-L/usr/local/Cellar/open-mpi/4.1.2/lib'
fftw_libs       = '-lm'#'-lfftw3f_mpi -lfftw3f_omp -lfftw3f_threads -lm'

f90_gps = 'get_powerspectrum.f90'
exe_gps = 'get_powerspectrum'

os.system('chmod -x {0}'.format( f90_gps ))
os.system('{0} {1} {2} {3} {4} {5} {6} {7} -o {8}'.format(
    fortrancompiler,fortranoptions,mpi_include,mpi_lib_dir,fftw_include,
    fftw_lib_dir,fftw_libs,f90_gps,exe_gps ))

os.system('./{0} << EOF\n"{1}"\nEOF'.format( exe_gps,field_file ))
