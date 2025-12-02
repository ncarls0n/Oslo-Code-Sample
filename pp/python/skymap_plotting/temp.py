from pixell import enmap, curvedsky, utils, powspec#, enplot
import numpy as np
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

gauss_file = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.02.05_SBsuite/ng0/cenz6500Mpc/maps/1000Mpc_n580_nb40_nt2_tsz_76543_fs.map'
filein = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.02.05_SBsuite/ng7_m40_test/cenz6500Mpc/maps/1000Mpc_n580_nb40_nt2_tsz_76543_fs.map'

# Read model from skymap fil
model='tSZ'
#model='kSZ'
#model='tau'
#model='kappa'
#model='CIB'


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

def bin_l(data,modlmap,bin_edges):
    digitized = np.digitize(np.ndarray.flatten(modlmap), bin_edges,right=True)
    return np.bincount(digitized,(data).reshape(-1))[1:-1]/np.bincount(digitized)[1:-1]

minell,maxell = 5.e4 , 1.e5
read_file=open(filein,"rb")
npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
fov  = float(fov/2/np.pi*360) #convert from rad to deg
ppmap = np.fromfile(read_file,dtype=np.float32,count=npix**2).reshape(npix,npix)
fov_box   = np.array([ [-fov/2,-fov/2] , [fov/2,fov/2] ]) * utils.degree
shape,wcs = enmap.geometry(pos=fov_box,res=fov/npix*utils.degree,proj='car')
map_IQU   = enmap.zeros(shape[-2:], wcs=wcs, dtype=np.float32) + ppmap
k_map_T = enmap.map2harm(map_IQU, nthread=0, normalize='phys')
TT_spec_2D = enmap.calc_ps2d(k_map_T)
modlmap       = TT_spec_2D.modlmap()
#minell,maxell = modlmap[np.nonzero(modlmap)].min() , modlmap.max()
bin_edges = np.arange(minell,maxell,40)  # 40 bins for ℓ in (0,1000)
ell_bin_centers = (bin_edges[1:] + bin_edges[:-1])/2.
binned_power = bin_l(TT_spec_2D,modlmap,bin_edges)
l   = ell_bin_centers
C_l = binned_power
D_l = 1e12*l*(l+1)*C_l/(2*np.pi)


read_file=open(gauss_file,"rb")
npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
fov  = float(fov/2/np.pi*360) #convert from rad to deg
ppmap = np.fromfile(read_file,dtype=np.float32,count=npix**2).reshape(npix,npix)
fov_box   = np.array([ [-fov/2,-fov/2] , [fov/2,fov/2] ]) * utils.degree
shape,wcs = enmap.geometry(pos=fov_box,res=fov/npix*utils.degree,proj='car')
map_IQU   = enmap.zeros(shape[-2:], wcs=wcs, dtype=np.float32) + ppmap
k_map_T = enmap.map2harm(map_IQU, nthread=0, normalize='phys')
TT_spec_2D = enmap.calc_ps2d(k_map_T)
modlmap       = TT_spec_2D.modlmap()
bin_edges = np.arange(minell,maxell,40)  # 40 bins for ℓ in (0,1000)
ell_bin_centers = (bin_edges[1:] + bin_edges[:-1])/2.
binned_power = bin_l(TT_spec_2D,modlmap,bin_edges)
lgauss   = ell_bin_centers
C_lgauss = binned_power
D_lgauss = 1e12*l*(l+1)*C_l/(2*np.pi)




# # Read some CAMB power spectrum
# ps,_ = powspec.read_camb_scalar("test_scalCls.dat")

# Plot the C_ls
fig,ax = plt.subplots()

if model=='tSZ':

    ax.plot(lgauss,D_lgauss, label='Gaussian')
    
    ax.plot(l, D_l, ls=':', label='non-Gaussian')

    # Plot labels
    xlabel = r'$\ell$'
    ylabel = r'$10^{12}\ell(\ell+1)C_\ell^{yy}/(2\pi)$'
else:
    ax.plot(ell_bin_centers, binned_power, label='WebSky')
# plt.plot(np.arange(12,1000), ps[0,1,12:1000], "-", label="theory")


ax.set_xscale('log')
ax.set_yscale('log')
plt.legend()
# plt.yscale("log")
plt.ylabel("Power $\sim C_{\ell}^{TT}$")
plt.xlabel("Multipole moment $\sim \ell$")
#plt.show()
plt.savefig('C_ell_TT'+model+'.pdf')
