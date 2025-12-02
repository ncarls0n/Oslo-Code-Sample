import sys,os

# Fortran compiler options
compiler  = 'ifort'
options   = '-04 -w'

# Paths
pks2map   = '.'
src       = pks2map   + '/..'
peakpatch = src       + '/..'
tables    = peakpatch + '/tables'

# Fortran scripts
cosmology = pks2map + '/cosmology.f90'
bbps_prof = pks2map + '/bbps_profile.f90'
inte_prof = pks2map + '/integrate_profiles.f90'
profiles  = pks2map + '/profiles.f90'
maptable  = pks2map + '/maptable.f90'
textlib   = src     + '/modules/External/textlib.f90'
program   = pks2map + '/make_maptable.f90'
scripts   = '{0} {1} {2} {3} {4} {5} {6}'.format(cosmology,bbps_prof,inte_prof,profiles,maptable,textlib,program)
exe       = pks2map + '/make_maptable'

# Read output table file
if len(sys.argv)==1:
    bbps_out = tables+'/bbps_model_1.tab'
else:
    bbps_out = sys.argv[1]

# If bbps_out exits, erase it
os.system( 'if [ -f {0} ] ; then rm -f {0} ; fi'.format(bbps_out) )

# Compile make_maptable.f90
os.system( 'cd {0};rm -rf *.o *.mod'.format(pks2map) )
os.system( 'cd {0};{1} {2} {3} -o {4}'.format( pks2map, compiler, options, scripts, exe ) )

# Run executable make_maptable
os.system( 'cd {0};./{1} {2} 1'.format(pks2map,exe,bbps_out) )

# Cleanup
os.system( 'cd {0};rm -rf *.o *.mod'.format(pks2map) )
