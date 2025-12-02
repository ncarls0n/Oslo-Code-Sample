from pixell import enmap, curvedsky, utils, powspec#, enplot
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex']=True
import sys,os
import healpy as hp

###########################################################################
#                                                                         #
#                 HARMONIC POWER SPECTRA FROM FlatSky MAPS                #
#                                                                         #
# This script reads in a WebSky flat-sky map covering a roughly square    #
# field of view of `fov` by `fov` square degrees and computes the angular #
# power spectrum coefficients $C_\ell$ for the given response function.   #
#                                                                         #
# USAGE:                                                                  #
#     python3 <...>/peak-patch/tools/get_websky_Cl.py <...>/<map>.        #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Read in HEALPix or flat-sky skymap file
filein = sys.argv[1]

# Read model from skymap file
if   'tsz' in filein.lower(): model='tSZ'
elif 'ksz' in filein.lower(): model='kSZ'
elif 'tau' in filein.lower(): model='tau'
elif 'kap' in filein.lower(): model='kappa'
elif 'cib' in filein.lower(): model='CIB'
else:
    raise Warning('type of map passed at command line not recognized, assu'
        +'ming map\ntype is tSZ Compton-y map.')

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
pkpdir = os.path.dirname(os.path.realpath(__file__))+'/../../../'

# Convert to dimensionless C_l
def x(nu):    return h*nu/(k_B*T_CMB)
def dT_y(nu): return T_CMB*(x(nu)*(np.exp(x(nu))+1)/(np.exp(x(nu))-1)-4)

# Bin by 
def bin_l(data,modlmap,bin_edges):
    digitized = np.digitize(np.ndarray.flatten(modlmap), bin_edges,right=True)
    return np.bincount(digitized,(data).reshape(-1))[1:-1]/np.bincount(digitized)[1:-1]

if model=='CIB':
    # Read in HEALPix formatted 1D array
    cibmap = hp.read_map(filein)
    nmap   = len(cibmap)
    nside  = int( (n/12)**.5 )
    fov = sys.argv[2]

    # Angular coordinates for each HEALPix pixel
    theta,phi = hp.pixelfunc.pix2map(nside,np.array(cibmap)).T

    # Angle from z axis sweeping toward the x and y axes
    theta_x   = np.arccos(np.sin(theta)*np.cos(phi))
    theta_y   = np.arccos(np.sin(theta)*np.sin(phi))

    # Mask out region of HELPix map outside the square field of view
    cibmap[ np.argwhere(theta_x>fov/2) ] = np.UNSEEN
    cibmap[ np.argwhere(theta_y>fov/2) ] = np.UNSEEN
    hp.pixelfunc.mask_bad(cibmap)

    # Compute harmonic power spectrum
    C_l = hp.sphtfunc.anafast(cibmap)
    l   = np.arange(len(C_l))
    D_l = 1e12*l*(l+1)*C_l/(2*np.pi)

else:
    # Read in data
    read_file=open(filein,"rb")
    # Read in 16 byte header (2 ints, 2 floats)
    npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
    npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
    fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
    fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
    fov  = float(fov/2/np.pi*360) #convert from rad to deg

    # Read map in peak patch format
    ppmap = np.fromfile(read_file,dtype=np.float32,count=npix**2).reshape(npix,npix)
    
    # Field of view 
    fov_box   = np.array([ [-fov/2,-fov/2] , [fov/2,fov/2] ]) * utils.degree
    shape,wcs = enmap.geometry(pos=fov_box,res=fov/npix*utils.degree,proj='car')
    map_IQU   = enmap.zeros(shape[-2:], wcs=wcs, dtype=np.float32) + ppmap

    #k_map_TEB = enmap.map2harm(map_IQU, nthread=0, normalize='phys')
    #k_map_T   = k_map_TEB[0]
    #k_map_E   = k_map_TEB[1]
    #k_map_B   = k_map_TEB[2]
    k_map_T = enmap.map2harm(map_IQU, nthread=0, normalize='phys')

    # 2D TT harmonic power spectrum
    TT_spec_2D = enmap.calc_ps2d(k_map_T)
    #TE_spec_2D = enmap.calc_ps2d(k_map_T, k_map_E)

    # Binning to make 1D C_\ell power spectrum
    modlmap       = TT_spec_2D.modlmap()
    # minell,maxell = modlmap[np.nonzero(modlmap)].min() , modlmap.max()
    minell,maxell = 5.e4 , 1.e5
    bin_edges = np.arange(minell,maxell,40)  # 40 bins for ℓ in (0,1000)
    ell_bin_centers = (bin_edges[1:] + bin_edges[:-1])/2.
    binned_power = bin_l(TT_spec_2D,modlmap,bin_edges)

    l   = ell_bin_centers
    C_l = binned_power
    D_l = 1e12*l*(l+1)*C_l/(2*np.pi)

# Plot the C_ls
fig,ax = plt.subplots()

if model=='tSZ':
    ax.plot(l, D_l, label='WebSky')

    # Plot labels
    xlabel = r'$\ell$'
    ylabel = r'$10^{12}\ell(\ell+1)C_\ell^{yy}/(2\pi)$'

    compare2data=False
    if compare2data==True:
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

plt.xlabel("Multipole moment $\sim \ell$")
elif model=='CIB':
    ax.plot(ell_bin_centers, binned_power, label='WebSky')

ax.set_xscale('log')
ax.set_yscale('log')
plt.legend()
ax.set_ylabel("C_l or whateber")
plt.savefig('C_ell_TT'+model+'.pdf')
###########################################################################
#                                                                         #
# COMPUTE HARMONIC POWER SPECTRA C_l                                      #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#if   model=='tSZ':   skymap=np.log10(map+1e-9)
#elif model=='kSZ':   skymap=skymap*T_CMB
#elif model=='kappa': skymap=hp.smoothing(skymap,fwhm=np.radians(1./60*15))

# Compute harmonic power spectrum
# C_l = hp.sphtfunc.anafast(skymap)
# l   = np.arange(len(C_l))
# D_l = 1e12*l*(l+1)*C_l/(2*np.pi)
# 
# fig,ax = plt.subplots(nrows=1,ncols=1)
# ax.set_xscale('log')
# ax.set_yscale('log')
# 
# if model=='tSZ':
# 
#     ax.plot( l, D_l, label=r'WebSky simulation' )
# 
#     # Plot labels
#     xlabel = r'$\ell$'
#     ylabel = r'$10^{12}\ell(\ell+1)C_\ell^{yy}/(2\pi)$'
# 
#     # Read in C_l data from Planck 2015 Results XXII Table 2
#     planck_file = pkpdir+'/tables/planck2015_tsz_cl.txt'
#     planck_data = np.loadtxt(planck_file,delimiter='\t')
#     planck_ebar = ax.errorbar( planck_data[:,2], planck_data[:,3],
#         yerr=(planck_data[:,4]**2+planck_data[:,5]**2)**.5, ls='none',
#         marker='o', ms=8, mew=.5, mec='w', mfc='k', elinewidth=2,
#         ecolor='k', label=r'$Planck$ 2015' )
#     planck_ebar[-1][0].set_linestyle('-')
# 
#     # Read in C_l data from Boillet et al 2018 (arXiv:1712.00788) Table 4
#     boillet_file = pkpdir+'/tables/boillet2018_tsz_cl.txt'
#     boillet_data = np.loadtxt(boillet_file,delimiter='\t')
#     boillet_ebar = ax.errorbar( boillet_data[:,0], boillet_data[:,1],
#         yerr=boillet_data[:,2], ls='none', marker='p', ms=8,
#         mew=.5, mec='w', mfc='tab:gray', elinewidth=1, ecolor='tab:gray',
#         label=r'Boillet $et$ $al.$ 2018' )
#     planck_ebar[-1][0].set_linestyle('-')
# 
#     # ACT
#     nu_act   = 148e9 # Hz
#     D_act    = 3.4e-12 * dT_y(nu_act)**-2
#     dD_act   = 1.4e-12 * dT_y(nu_act)**-2
#     act_ebar = ax.errorbar( 3000, 1e12*D_act, yerr=1e12*dD_act, ls='none',
#         marker='s', ms=8, mew=.5, mec='w', mfc='tab:red',
#         elinewidth=1, ecolor='tab:red', label=r'ACT 2013' )
# 
#     # SPT
#     nu_spt   = 143e9 # Hz
#     D_spt    = 4.08e-12 * dT_y(nu_spt)**-2
#     lD_spt   = .67e-12  * dT_y(nu_spt)**-2
#     uD_spt   = .58e-12  * dT_y(nu_spt)**-2
#     spt_ebar = ax.errorbar( 3000, 1e12*D_spt,
#         yerr=1e12*np.array([(lD_spt,uD_spt)]).T,
#         ls='none', marker='X', ms=8, mew=.5, mec='w', mfc='tab:green',
#         elinewidth=1, ecolor='tab:green', label=r'SPT 2015' )
# 
# elif model=='kSZ':
#     
#     ax.plot( l, D_l, label=r'WebSky $f_{\text{NL}}=1000$' )
#     # Plot labels
#     xlabel=r'$\ell$'
#     ylabel=r'$\ell(\ell+1)C_\ell^{TT}/(2\pi)$ $[\mu$K$]$'
# 
#     skymap = hp.read_map('/Users/Nate/Desktop/ng2_fnl1e2/maps/4000Mpc_n236_nb20_nt10_ksz_13579_hp.fits')
#     C_l = hp.sphtfunc.anafast(skymap)
#     l   = np.arange(len(C_l))
#     D_l = 1e12*l*(l+1)*C_l/(2*np.pi)
#     ax.plot( l, D_l, label=r'WebSky $f_{\text{NL}}=100$' )
# 
#     skymap = hp.read_map('/Users/Nate/Desktop/ng2_fnl0/maps_nside2048/4000Mpc_n236_nb20_nt10_ksz_13579_hp.fits')
#     C_l = hp.sphtfunc.anafast(skymap)
#     l   = np.arange(len(C_l))
#     D_l = 1e12*l*(l+1)*C_l/(2*np.pi)
#     ax.plot( l, D_l, label=r'WebSky gaussian' )
# 
# elif model=='ksz':
#     # Plot labels
#     xlabel=r'$\ell$'
#     ylabel=r'$\ell(\ell+1)C_\ell^{TT}/(2\pi)$'
# 
# ax.set_xlim(5,5000)
# # ax.set_ylim(5e-4,2)
# ax.set_xlabel(xlabel)
# ax.set_ylabel(ylabel)
# plt.legend()
# plt.show()













