import os,sys
import numpy as np

# Usage:
# 
#     python3 make_powerspectrum.py s_box n n_buff

s_box  = float( sys.argv[1] )
n      =   int( sys.argv[2] )
n_buff =   int( sys.argv[3] )

# Load the relevant columns of the correlation matrix file
# column -17 is k, colun -6 is P_chichi(k)
p = np.loadtxt('corr.out',usecols=(-17,-6))
p = p[1:-1]
#os.system('rm -f corr.out')

# Scale DFT frequencies such that they sample the wavenumber of the
# corresponding continuous field
n_k    = len(p)
dk     = 2*np.pi/s_box
k_max  = (n/2+1)*dk
p[:,0] = (p[:,0]-p[0,0])/(p[-1,0]-p[0,0]) * k_max + dk/2

# Scale Fourier modes such that they sample the power spectrum of the
# corresponding continuous field
p[:,1] = ( 1/(n*dk*(n-2*n_buff)) )**3 * p[:,1]

# Save field as 32-bit floating point numbers
np.savetxt('corr.dat',p,fmt='%.5e')

