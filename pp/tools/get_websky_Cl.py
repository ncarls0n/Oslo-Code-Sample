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
if len(sys.argv)!=2:
    raise SyntaxError('You\'ve passed the wrong number of arguments at com'
        +'mand line.\nCorrect usage:\n    python3 <...>/peak-patch//tools/'
        +'get_websky_Cl.py <...>/<map>.fits')
filein = sys.argv[1]

# Read model from skymap file
if   'tsz' in filein.lower(): model='tSZ'
elif 'ksz' in filein.lower(): model='kSZ'
elif 'tau' in filein.lower(): model='tau'
elif 'kap' in filein.lower(): model='kappa'
else:
    raise Warning('type of map passed at command line not recognized, assu'
        +'ming map\ntype is tSZ Compton-y map.')

# Read skymap file
skymap = hp.read_map(filein)

###########################################################################
#                                                                         #
# CONSTANTS & CODE PARAMETERS                                             #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# 2018 CODATA - NIST SP 959 (June 2019) 
h   = 6.62607015e-34 # J s
k_B = 1.380649e-23 # J/K

# CMB temperature from Planck
T_CMB = 2.7255 # K

# Peak Patch directory
pkpdir = os.path.dirname(os.path.realpath(__file__))+'/..'

# Convert to dimensionless C_l
def x(nu):    return h*nu/(k_B*T_CMB)
def dT_y(nu): return T_CMB*(x(nu)*(np.exp(x(nu))+1)/(np.exp(x(nu))-1)-4)

###########################################################################
#                                                                         #
# COMPUTE HARMONIC POWER SPECTRA C_l                                      #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#if   model=='tSZ':   skymap=np.log10(map+1e-9)
#elif model=='kSZ':   skymap=skymap*T_CMB
#elif model=='kappa': skymap=hp.smoothing(skymap,fwhm=np.radians(1./60*15))

# Compute harmonic power spectrum
C_l = hp.sphtfunc.anafast(skymap)
l   = np.arange(len(C_l))
D_l = 1e12*l*(l+1)*C_l/(2*np.pi)

fig,ax = plt.subplots(nrows=1,ncols=1)
ax.set_xscale('log')
ax.set_yscale('log')

if model=='tSZ':

    ax.plot( l, D_l, label=r'WebSky simulation' )

    # Plot labels
    xlabel = r'$\ell$'
    ylabel = r'$10^{12}\ell(\ell+1)C_\ell^{yy}/(2\pi)$'

    # Read in C_l data from Planck 2015 Results XXII Table 2
    planck_file = pkpdir+'/tables/planck2015_tsz_cl.txt'
    planck_data = np.loadtxt(planck_file,delimiter='\t')
    planck_ebar = ax.errorbar( planck_data[:,2], planck_data[:,3],
        yerr=(planck_data[:,4]**2+planck_data[:,5]**2)**.5, ls='none',
        marker='o', ms=8, mew=.5, mec='w', mfc='k', elinewidth=2,
        ecolor='k', label=r'$Planck$ 2015' )
    planck_ebar[-1][0].set_linestyle('-')

    # Read in C_l data from Boillet et al 2018 (arXiv:1712.00788) Table 4
    boillet_file = pkpdir+'/tables/boillet2018_tsz_cl.txt'
    boillet_data = np.loadtxt(boillet_file,delimiter='\t')
    boillet_ebar = ax.errorbar( boillet_data[:,0], boillet_data[:,1],
        yerr=boillet_data[:,2], ls='none', marker='p', ms=8,
        mew=.5, mec='w', mfc='tab:gray', elinewidth=1, ecolor='tab:gray',
        label=r'Boillet $et$ $al.$ 2018' )
    planck_ebar[-1][0].set_linestyle('-')

    # ACT
    nu_act   = 148e9 # Hz
    D_act    = 3.4e-12 * dT_y(nu_act)**-2
    dD_act   = 1.4e-12 * dT_y(nu_act)**-2
    act_ebar = ax.errorbar( 3000, 1e12*D_act, yerr=1e12*dD_act, ls='none',
        marker='s', ms=8, mew=.5, mec='w', mfc='tab:red',
        elinewidth=1, ecolor='tab:red', label=r'ACT 2013' )

    # SPT
    nu_spt   = 143e9 # Hz
    D_spt    = 4.08e-12 * dT_y(nu_spt)**-2
    lD_spt   = .67e-12  * dT_y(nu_spt)**-2
    uD_spt   = .58e-12  * dT_y(nu_spt)**-2
    spt_ebar = ax.errorbar( 3000, 1e12*D_spt,
        yerr=1e12*np.array([(lD_spt,uD_spt)]).T,
        ls='none', marker='X', ms=8, mew=.5, mec='w', mfc='tab:green',
        elinewidth=1, ecolor='tab:green', label=r'SPT 2015' )

elif model=='kSZ':
    
    ax.plot( l, D_l, label=r'WebSky $f_{\text{NL}}=1000$' )
    # Plot labels
    xlabel=r'$\ell$'
    ylabel=r'$\ell(\ell+1)C_\ell^{TT}/(2\pi)$ $[\mu$K$]$'

    skymap = hp.read_map('/Users/Nate/Desktop/ng2_fnl1e2/maps/4000Mpc_n236_nb20_nt10_ksz_13579_hp.fits')
    C_l = hp.sphtfunc.anafast(skymap)
    l   = np.arange(len(C_l))
    D_l = 1e12*l*(l+1)*C_l/(2*np.pi)
    ax.plot( l, D_l, label=r'WebSky $f_{\text{NL}}=100$' )

    skymap = hp.read_map('/Users/Nate/Desktop/ng2_fnl0/maps_nside2048/4000Mpc_n236_nb20_nt10_ksz_13579_hp.fits')
    C_l = hp.sphtfunc.anafast(skymap)
    l   = np.arange(len(C_l))
    D_l = 1e12*l*(l+1)*C_l/(2*np.pi)
    ax.plot( l, D_l, label=r'WebSky gaussian' )
elif model=='ksz':
    # Plot labels
    xlabel=r'$\ell$'
    ylabel=r'$\ell(\ell+1)C_\ell^{TT}/(2\pi)$'



ax.set_xlim(5,5000)
#ax.set_ylim(5e-4,2)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.legend()
plt.show()













