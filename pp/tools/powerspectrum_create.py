# powerslectrum_create.py

###########################################################################
##  This script uses CAMB (Code for Anisotropies in the Microwave        ##
##  Background) to make tables of power spectra and transfer functions   ##
##  used to construct initial condition fields for Peak Patch.           ##
###########################################################################

import numpy as np
import pylab as py
import matplotlib.pyplot as plt
#from hmf_parameters import *
#from hmf_pp import * missing modules used for calculating sigma8, no necessary to make power spectra
from scipy.interpolate import *

# When I downloaded CAMB using `pip install camb`, it installed for
# python3.9, but I am still mostly using python3.8.2 and haven't gotten
# arround to updating so I had to append it to the system path (as just
# ading it to PATH in terminal was causing other problems).
import subprocess
machine = subprocess.check_output('hostname')
if machine in [b'homes-MacBook-Air.local\n',b'Nates-MacBook-Pro.local\n']:
    import sys
    sys.path.append('/usr/local/lib/python3.9/site-packages/')

# Import CAMB
import camb
from camb import model, initialpower
# modules for pycamb:
# module load gcc/5.2.0 intel/15.0.2 openmpi/gcc/1.10.3 python/2.7.8


###########################################################################
##  Set cosmology to Plack's latest results, see "Planck 2018 results.   ##
##  VI Cosmological Parameters", Table 2, column TT,TE,EE+lowE+lensing   ##
###########################################################################

# Cosmological energy density fractions
omegab = 0.0493        # of baryonic matter
omegac = 0.2645        # of cold dark matter
omegam = omegab+omegac # of all energy that clusters under gravity
omegal = 1-omegam      # of cosmological constant energy
omegak = 0.0

# Additional parameters
h      = 0.6735 # dimensionless Hubble const. "little h", H_0/100 km/s/Mpc
ns     = 0.9649 # spectral index
sigma8 = 0.8111 # standard deviation of field smoothed at 8 Mpc/h
rho_c  = 2.7754e11 # critical energy density in units of h^2 M_sol Mpc^-3
rho    = rho_c*omegam*h**2 # average matter-like energy density at present
z      = 0.0 # redshift at the present
As     = 2.100e-9 # primordial super-horizon comoving curvature power spectrum amplitude
tau    = 0.0544 # optical depth of reionization
m_nu   = 0.0 # massive neutrino mass

###########################################################################
##  Initialize CAMB calculation of power spectra and transfer functions  ##
##  by creating an object of class CAMBparams, pars, and setting it to   ##
##  the values corresponding to our cosmology.                           ##
###########################################################################

# Range of wavenumbers to be used
kmax_h   = 5.0e3/h  # max wavenumber / h
kmin_h   = 5.0e-6/h # min wavenumber / h
nkpoints = 1000     # number of points in power spectra tables

# Set cosmological parameters for CAMB power spectra in above cosmology
pars = camb.CAMBparams() # Make parameter object
pars.set_cosmology(H0=h*100, ombh2=omegab*h**2, omch2=omegac*h**2,
    mnu=m_nu, omk=omegak, tau=tau) # Set with above cosmology
pars.set_dark_energy() # re-set default DE equation of state
pars.InitPower.set_params(ns=ns, As=As) # sets primordial power spectrum

# Set parameters for calculating linear matter power spectrum at redshift z
# and maximum wavenumber kmax_h*h
pars.set_matter_power(redshifts=[z], kmax=kmax_h*h)
pars.NonLinear = model.NonLinear_none # linear matter power spectrum only
pars.PK_WantTransfer = 1 # Tell CAMB to calculate the matter
pars.WantTransfer    = 1 #     power transfer function
pars.Transfer.kmax   = kmax_h*h # Maximum wavenumber for transfer function


###########################################################################
##  Compute matter power spectrum P(k)                                   ##
###########################################################################

# Calculate results for specified parameters
results = camb.get_results(pars) # object of class CAMBdata

# Calcualte matter power spectrum for CAMBdata object results
kh, zz, pk = results.get_matter_power_spectrum(minkh=kmin_h,
    maxkh=kmax_h, npoints=nkpoints)
# kh[j] is an array of wavenumber k/h
# zz[i] is an array of redshift z
# pk[i,j] is a matrix of powerspectrum as a function of k/h and z

# Get sigma_8 (standard deviation for overdensity smoothed at 8 Mpc/h)
# value "today", so z=0. Must be called after get_matter_power_spectrum()
s8 = np.array(results.get_sigma8())


###########################################################################
##  Compute Transfer function T()                                   ##
###########################################################################

# Calculate matter transfer function
transfer = results.get_matter_transfer_data()

# Transfer functions as NumPy array Tf[n,i,j]
Tf = transfer.transfer_data # transfer funcs as function of z and q-modes
# divided by k^2 so that they are roughly constant at low k on super-
# horizon scales
# where n=6 is the total matter power, index i corresponds to the index of
# q modes calculated, and j corresponds to the index of redshift z[j]
kk = transfer.q/h # q mode over h
Tf_m = Tf[6,:,0]  # total matter power at z=0

# Normalize matter power spectrum by setting sigma_8 equal to that of CMB
print("sigma_8 pre normalization = ", s8)
norm = (sigma8/s8)**2 # normalization constant
k    = kh * h         # exact wavenumber k, not scaled by h
pk   = norm * pk[0,:] / ( 2.*np.pi * h)**3 # normalized P_m(z=0,k)
# Note here that we divide by h^3 to get pk in units of Mpc^3 as opposed to (Mpc/h)^3, the factor of (2pi)^3 is just a weird factor that Peak Patch expects to be there for no particular reason...

# #plt.loglog(k,pk)
# #plt.show()
# # Calculate mass and sigma
# Ma, sigma = calculate_sigma(k, pk, omegam, h)
# sigma_M = interp1d(10**Ma,sigma)
# print "sigma_8 post normalization is ", sigma_M(4./3*np.pi*(8.0/h)**3*rho)

# Primordial zeta power spectrum
ko     = 0.05
pkzeta = 2*np.pi**2*As/k**3 * (k/ko)**(ns-1)

# Light field power spectra
Achi = (5.e-7)**2
pkchi = 2*np.pi**2*Achi/k**3 #in units of sigmas
pkchi = pkchi/(2*np.pi)**3 #for pp power spectra
#Get transfer function
Trans = np.sqrt(pk/pkzeta)

# New model (Nate, 22 April 2022) spatially localized intermittent non
# -Gaussianity
#A2     = 1.6e-19 
#R2     = 6.4e-1  # Mpc/h
#pkchi2 = A2*(k*R2)**2 * np.exp( -k**2*R2**2 )
#pkchi2 = 2*np.pi**2*k**-3 * pkchi2

#R   = (10**Ma*3/4/np.pi/3.4e10)**(1./3)

np.savetxt("power.dat",np.transpose([k,pk,Trans,pkchi]),fmt='%1.4e')
#np.savetxt("power.dat",np.transpose([k,pk,Trans,pkchi,pkchi2]),fmt='%1.4e')

#np.savetxt("sigma.dat",np.transpose([Ma,sigma]),fmt='%1.4e')
