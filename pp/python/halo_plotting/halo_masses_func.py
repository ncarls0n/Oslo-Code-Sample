import numpy as np
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
import sys
import os

Omx=.2645
OmB=.0493
h=.6735
boxsize=4000.

# Additional parameters
Omm      = Omx+OmB        # Omega_m, fraction of energy that can cluster

m_per_Mpc = 3.0856775814913673e22 # 3.08...e22 m = 1 Mpc
H_0       = h*(100*1e3/m_per_Mpc) # s^-1
G         = 4.51724e-48           # Mpc^3 Msol^-1 s^-2
rhocrit   = 3*H_0**2/(8*np.pi*G)  # Msol Mpc^-3
rho_m     = rhocrit*Omm           # Msol Mpc^-3

def halo_mass(R): # R is comoving radius in Mpc
    return 4/3*np.pi*R**3*rho_m # halo mass in Msol

# Plot n(m)
fig,ax = plt.subplots(nrows=1,ncols=1)

merge0= '/mnt/scratch-lustre/njcarlson/peak-patch-runs/merges/4000Mpc_n236_nb20_nt10_merge.pksc.13579_fnl'
merge = [
    merge0+'1e4',
    merge0+'5e3',
    merge0+'1e3',
    merge0+'5e2',
    merge0+'0'
    ]

names = ['fNL=1e4','fNL=5e3','fNL=1e3','fNL=5e2','Gaussian']

filters    = np.fromfile( '{0}/filter.dat'.format(os.path.dirname(merge0)), sep=' ')
filters    = np.reshape( filters[1:], (int(filters[0]),4) )
filters    = filters[:,2]
m_halo_min = halo_mass( filters[0] )

for j in range(len(merge)):

    # Read size of halo catalogue in bytes
    size_bytes = os.path.getsize(merge[j])

    # Read catalogue header
    N_halos    = np.fromfile( merge[j], dtype=np.int32,count=1 )[0]
    RTHmax     = np.fromfile( merge[j], dtype=np.float32,count=1 )[0]
    zin        = np.fromfile( merge[j], dtype=np.float32,count=1 )[0]
    m_halo_max = halo_mass(RTHmax)

    # All catalogue parameters are 32-bit (or 4-byte) integers or floats, thus
    N_columns = int( (size_bytes-12)/4/N_halos ) # is the number of columns
    catalogue = np.fromfile(merge[j],dtype=np.float32,count=N_columns*N_halos)
    catalogue = np.reshape(catalogue,(N_halos,N_columns))
    m_halo    = halo_mass(catalogue[:,6]) # Msol
    del(catalogue) # Free up memory

    # Discretization
    m = np.logspace( np.log10(m_halo_min*.5), np.log10(m_halo_max*1.05), 100 )
    n_m = np.zeros(len(m))
    for i in range(N_halos):
        i_upper  = np.searchsorted(m,m_halo[i],side='left')
        if i_upper==0:
            n_m[i_upper]+=1
        elif i_upper>0 and i_upper<len(m):
            i_interp = i_upper - (m[i_upper]-m_halo[i])/(m[i_upper]-m[i_upper-1])
            n_m[i_upper]   += 1-(i_upper-i_interp)
            n_m[i_upper-1] +=    i_upper-i_interp
    n_m = n_m/boxsize**3

    ax.plot(m, n_m, label=names[j] )

## Press-Schechter formalism for halo mass function
#def M_star(M,ns,rho_m,sigma):
#    return ( (rho_m**((4-ns)/3)/(2*sigma**2))**(3/(2+ns))
#
#def dn_dM(M,ns,rho_m,sigma):
#    m_star = M_star(M,ns,rho_m,sigma)
#    return ( np.pi**-.5 * (2+ns)/3 * rho_m*M**-2 * (M/m_star)**((2+ns)/6)
#        * np.exp( -(M/m_star)**((2+ns)/3) ) )
#
#ax.plot(m, dn_dM(m,ns,rho_m,sigma), label=r'Press-Schechter')
#
# BCEK formalism for halo mass function
# def BCEK

print(m_halo_min,m_halo_min*.5)
print(m_halo_max,m_halo_max*1.05)

ax.set_xlim([m_halo_min*.5,m_halo_max*1.05])
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'$M/M_{{\odot}}$')
ax.set_ylabel(r'$dn/dM$ $[$Mpc$^{-3}]$')
ax.legend()
fig.set_size_inches(14,6)
fig.savefig('halo_masses_func.png',dpi=100)



