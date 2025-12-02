import numpy as np
import pyfftw
import sys,gc
from scipy.special import spherical_jn
import scipy.integrate as integrate
###########################################################################
#                                                                         #
#                  READ SIGMA-8 from SPECTRUM of 3D FIELD                 #
#                                                                         #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Cosmological parameters
# see Plank 2018 results. VI cosmological Parameters, Table 2, column
# TT,TE,EE+lowE+lensing
h = 0.6735 # dimensionless Hubble parameter H_0/(100 km/s/Mpc)

filein  = sys.argv[1]
if len(sys.argv)>2: kcol,pcol=int(sys.argv[1]),int(sys.argv[2])
else              : kcol,pcol=0,1

#neff    = int(sys.argv[2])
#sbox    = float(sys.argv[3])
#dk      = 2*np.pi/s
#kmax    = neff*dk
#nyquist = int(neff/2)+1
#klen  = int(3**.5*(nyquist-1))

data  = np.loadtxt( filein, delimiter=' ' )
k_in  = data[:,kcol]
Pk_in = data[:,pcol]
del(data)

# Real-space top-hat kernel
def W(u):return 3*u**-1*spherical_jn(1,u)
# Should really have the form:
# def W(u):return 3/(2*np.pi)**3*u**-1*spherical_jn(1,u)
# buf for some reason when I run it like that the result is off by more than just the factor of (2*np.pi)**3, which I can't make sense of...

# Power spectrum as a function using NumPy's interpolator
index = np.arange(len(k_in))
def P(k):return np.interp( np.interp(k,k_in,index), index, Pk_in )
    #k_in[i] <= k <= k_in[i+1]
    #i = np.searchsorted(k_in,k,'right')-1
    #return Pk_in[i]+(Pk_in[i+1]-Pk_in[i])/(k_in[i+1]-k_in[i])*(k-k_in[i])

# Read d(sigma_8^2)/dk
#def dvar8_dk(k):return (2*np.pi**2)**-1*(2*np.pi)**-6*k**2*W(8/h*k)**2*P(k)
def dvar8_dk(k):return (2*np.pi**2)**-1*k**2*W(8/h*k)**2*P(k)

import matplotlib.pyplot as plt
k_out = np.logspace(np.log10(k_in[1]), np.log10(k_in[-2]), 1000 )

fig,ax = plt.subplots(1)
ax.plot(k_out, P(k_out), label=r'$P(k)$ interpolated') 
ax.plot(k_out, W(8/h*k_out)**2, label=r'$\tilde{W}^2(k)$' )
ax.plot(k_out, dvar8_dk(k_out), label=r'$\frac{d\sigma^2}{dk}$')
#ax.plot(k_out, P(k_out)*W(8*k_out)**-2, label=r'$\tilde{W}^{-2}(k)\frac{d\sigma^2}{dk}$')
ax.plot(k_in, Pk_in, label=r'$P(k)$', ls='--')
ax.set_yscale('log')
ax.set_xscale('log')

sigma8 = (2*np.pi)**1.5*integrate.quad( dvar8_dk, k_in[0], k_in[-1], limit=100 )[0]**.5
print(sigma8)
plt.legend()
plt.show()
