import numpy as np
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
import sys,os
"""
USAGE:

cd peak-patch
python3 python/run_checks/halo_mass_function.py

"""
# Additional parameters
h       = .6735
Omega_m = .3138
rhocrit = 2.775e11*h**2
rho_m   = Omega_m * rhocrit
s_box   = 500. # Mpc

def halo_mass(R): # R is comoving radius in Mpc
    return 4/3*np.pi*R**3*rho_m # halo mass in Msol

#files = (r'/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.04.29_LIM_prod'
#    +r'uction_runs/ng{0}_cenz7150Mpc_C{1}/output/500Mpc_n432_nb'
#    +r'74_nt3_merge.pksc.13579_M_gt_{2}Msol')
files = (r'/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.04.29_LIM_prod'
    +r'uction_runs/ng{0}_cenz7150Mpc_C{1}/output/500Mpc_n432_nb'
    +r'74_nt3_merge.pksc.13579_M_gt_{2}Msol')
s   = '_mlambda'
a_e = np.array([ '1e-51', '2e-51', '2.23e-51', '2.5e-51', '2.75e-51', '3e-51', '4e-51' ])
halo_cats = [ files.format( 0, ''              , '1.129e+11' ) ,
              files.format( 7, s+'25_a1e-51'   , '1.129e+11' ) ,
              files.format( 7, s+'25_a2e-51'   , '1.129e+11' ) ,
              files.format( 7, s+'25'          , '1.129e+11' ) ,
              files.format( 7, s+'25_a2.5e-51' , '1.129e+11' ) ,
              files.format( 7, s+'25_a2.75e-51', '1.129e+11' ) ,
              files.format( 7, s+'25_a3e-51'   , '1.129e+11' ) ,
              files.format( 7, s+'25_a4e-51'   , '1.129e+11' ) ]

fig,ax = plt.subplots(2,1, gridspec_kw={'height_ratios':[3,1]})

for i in range(len(halo_cats)):

    in_cat = open( halo_cats[i] , 'rb' )

    # Read catalogue header
    N_halos    = np.fromfile( in_cat, dtype=np.int32,count=1 )[0]
    RTHmax     = np.fromfile( in_cat, dtype=np.float32,count=1 )[0]
    zin        = np.fromfile( in_cat, dtype=np.float32,count=1 )[0]
    m_halo_max = halo_mass(RTHmax)
    
    # All catalogue parameters are 32-bit (or 4-byte) integers or floats, thus
    size_bytes = os.path.getsize(halo_cats[i])
    N_columns = int( (size_bytes-12)/4/N_halos ) # is the number of columns
    catalogue = np.reshape( np.fromfile(in_cat,dtype=np.float32,count=N_columns*N_halos) , (N_halos,N_columns) )

    x,y,z = catalogue[:,:3].T
    m_halo = halo_mass(catalogue[:,6]) # Msol
    del(catalogue) # Free up memory

    n_bins      = 20
    #m_min,m_max = np.min(m_halo) , np.max(m_halo)
    m_min,m_max = 1.129e11 , 1e14
    bin_edges   = np.array([ m_min*(m_max/m_min)**(j/n_bins) for j in range(n_bins+1) ])
    bin_centres = (bin_edges[:-1]+bin_edges[1:])/2

    hist = np.histogram( m_halo, bins=bin_edges )[0]

    # Rather than dN=n(M), I want dN/dlnM = M dN/dM
    hist = bin_centres**2 * hist/(bin_edges[1:]-bin_edges[:-1]) / s_box**3

    if i==0:
        label=r'Gaussian'
        histG=hist
    else:
        label=r'$a_e = {0}\times10^{{-51}}$'.format(a_e[i-1][:-4])

    ax[0].plot( bin_centres[1:], hist[1:], label=label )
    ax[1].plot( bin_centres[1:], np.log(hist[1:]/histG[1:])        )

#ax.plot(m, n_m, label=r'Peak-Patch')
ax[0].set_xscale('log'), ax[0].set_yscale('log')
ax[0].set_xlim([1.5e11,1e14])
ax[1].set_xscale('log')# ax[1].set_yscale('log')
ax[1].set_xlabel(r'$M ~ [M_{{\odot}}]$')
ax[1].set_xlim([1.5e11,1e14])
ax[1].set_ylim([-.9,0.9])
#ax.set_ylabel(r'$n$')# $[$Mpc$^{-3}]$')
ax[0].set_ylabel(r'$M dN / d\ln M ~ [M_\odot \textrm{Mpc}^{-3}]$')# $[$Mpc$^{-3}]$')
ax[1].set_ylabel(r'$\ln \left[ \frac{ dN/d\ln M }{ dN_{gaussian}/d\ln M} \right]$')
ax[0].legend()
#fig.set_size_inches(14,6)
#fig.savefig(run_dir+'/halo_mass_function.png',dpi=100)
outfile = os.path.dirname(halo_cats[0])+'/n_of_m_vs_mlambda.pdf'
outfile = '/cita/d/www/home/njcarlson/cifar_meeting/n_of_m_vs_ae_sigma8planck.pdf'
fig.savefig(outfile,bbox_inches='tight')
