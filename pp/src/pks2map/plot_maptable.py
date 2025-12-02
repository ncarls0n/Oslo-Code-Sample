import numpy as np
import matplotlib.pyplot as plt 
import sys,os
plt.rcParams['text.usetex'] = True

# Paths
pks2map   = '.' 
src       = pks2map   + '/..'
peakpatch = src       + '/..'
tables    = peakpatch + '/tables'

# Read output table file
if len(sys.argv)==1:
    bbps_file = tables+'/bbps_model_1.tab'
else:
    bbps_file = sys.argv[1]

# Open the file
bbps_in = open( bbps_file, 'rb' )

# Read in header
[ N_chi, N_M, N_R ] = [ np.fromfile(bbps_in,dtype=np.int32,count=1)[0] for i in range(3) ]
[ chi_range, M_range, R_range ] = [ np.fromfile(bbps_in,dtype=np.float32,count=2) for i in range(3) ]

chi = np.logspace( np.log10(chi_range[0]) , np.log10(chi_range[1]) , N_chi )
M   = np.logspace( np.log10(  M_range[0]) , np.log10(  M_range[1]) , N_M   )   
R   = np.logspace( np.log10(  R_range[0]) , np.log10(  R_range[1]) , N_R   )   

# Read in table as 3 × N_chi × N_M × N_R array
table = np.reshape( np.fromfile(bbps_in,dtype=np.float32,count=3*N_chi*N_M*N_R) , (N_R,N_M,N_chi,3) ).T 
dchi  = np.concatenate( ( chi[0:1] , chi[1:]-chi[:-1] ) )
dy    = table[0,:,-1,-1] 
del(table)

fig,ax = plt.subplots()
ax.plot( chi, dy/dchi**3 )

bbps_in = open( tables+'/bbps_1.tab', 'rb' )
[ N_chi, N_M, N_R ] = [ np.fromfile(bbps_in,dtype=np.int32,count=1)[0] for i in range(3) ]
[ chi_range, M_range, R_range ] = [ np.fromfile(bbps_in,dtype=np.float32,count=2) for i in range(3) ]
chi = np.logspace( np.log10(chi_range[0]) , np.log10(chi_range[1]) , N_chi )
M   = np.logspace( np.log10(  M_range[0]) , np.log10(  M_range[1]) , N_M   )
R   = np.logspace( np.log10(  R_range[0]) , np.log10(  R_range[1]) , N_R   )
table = np.reshape( np.fromfile(bbps_in,dtype=np.float32,count=3*N_chi*N_M*N_R) , (N_R,N_M,N_chi,3) ).T
dchi  = np.concatenate( ( chi[0:1] , chi[1:]-chi[:-1] ) )
dy    = table[0,:,-1,-1] 
ax.plot( chi, dy/dchi**3 )

ax.set_xlabel(r'comoving distance $\chi$ [Mpc]')
ax.set_ylabel(r'Differential Compton-y $d^3y/d\chi^3$')
ax.set_xscale('log')
ax.set_yscale('log')

plt.show()
