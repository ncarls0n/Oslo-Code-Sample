import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex']=True
import sys,os

###########################################################################
#                                                                         #
#                 HARMONIC POWER SPECTRA FROM WebSky MAPS                 #
#                                                                         #
# This script reads in a WebSky map and runs the python wrapper for the   #
# HEALPix module anafast.f90, which computes the harmonic power spectra,  #
# $C_\ell$ for that map.                                                  #
#                                                                         #
# USAGE:                                                                  #
#     python3 <...>/peak-patch/tools/get_websky_Cl.py <...>/<map>.fits    #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Read in HEALPix skymap file

g_dir = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.02.05_SBsuite/ng0/cenz6500Mpc/maps/'
ng_dir= '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.02.05_SBsuite/ng7_m40_test/cenz6500Mpc/maps/'
files = [ 'cib100.fits','cib143.fits','cib217.fits','cib353.fits','cib545.fits','cib857.fits' ]
labels= [ '100 GHz',    '143GHz'     ,'217 GHz',    '353 GHz'    ,'545 GHz'    ,'857 GHz'     ]
# Read skymap file

h   = 6.62607015e-34 # J s
k_B = 1.380649e-23 # J/K
T_CMB = 2.7255 # K
def x(nu):    return h*nu/(k_B*T_CMB)
def dT_y(nu): return T_CMB*(x(nu)*(np.exp(x(nu))+1)/(np.exp(x(nu))-1)-4)

fig,ax = plt.subplots(nrows=1,ncols=1)
ax.set_xscale('log')
ax.set_yscale('log')

for j in range(1):#len(labels)):
    skymap = hp.read_map( g_dir+files[j] )

    # Compute harmonic power spectrum
    C_l = hp.sphtfunc.anafast(skymap)
    l   = np.arange(len(C_l))
    D_l = 1e12*l*(l+1)*C_l/(2*np.pi)

    ax.plot(l,D_l,color='k')# , label = labels[j] )

#for j in range(len(labels)):
#    skymap = hp.read_map( ng_dir+files[j] )
#    C_l = hp.sphtfunc.anafast(skymap)
#    l   = np.arange(len(C_l))
#    D_l = 1e12*l*(l+1)*C_l/(2*np.pi)
#
#    ax.plot(l,D_l,ls=':',color='tab:red', label = labels[j])


plt.legend()
plt.ylabel("Power $\sim C_{\ell}^{TT}$")
plt.xlabel("Multipole moment $\sim \ell$")
plt.savefig('C_ell_TT_CIB.pdf')





