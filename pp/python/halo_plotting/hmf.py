import numpy as np
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
import sys,os
"""
USAGE:

cd peak-patch
python3 python/run_checks/halo_mass_function.py <...>/<run_name>

"""
# Additional parameters
h       = .6735
Omega_m = .3138
rhocrit = 2.775e11*h**2
rho_m   = Omega_m * rhocrit

def halo_mass(R): # R is comoving radius in Mpc
    return 4/3*np.pi*R**3*rho_m # halo mass in Msol

files = (r'/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.04.29_LIM_prod'
    +r'uction_runs/ng7_cenz7150Mpc_C_mlambda{0}_a{1}/output/500Mpc_n432_nb'
    +r'74_nt3_merge.pksc.13579_M_gt_{2}Msol'

a_e = np.array([ 1.e-51, 2e-51, 2.236e-51, 2.75e-51, 3e-51 ])
halo_cats = [ files.format( '25', '1e-51'    , '1.129e11' ) ,
              files.format( '25', '2e-51'    , '1.129e11' ) ,
              files.format( '25', '2.236e-51', '1.129e11' ) ,
              files.format( '25', '2.75e-51' , '1.129e11' ) ,
              files.format( '25', '3e-51'    , '1.129e11' ) ,
              files.format( '25', '4e-51'    , '1.129e11' ) ]

for i in range(len(halo_cats)):

    in_cat = open( halo_cats[i] , 'rb' )

    # Read catalogue header
    N_halos    = np.fromfile( in_merge, dtype=np.int32,count=1 )[0]
    RTHmax     = np.fromfile( in_merge, dtype=np.float32,count=1 )[0]
    zin        = np.fromfile( in_merge, dtype=np.float32,count=1 )[0]
    m_halo_max = halo_mass(RTHmax)
    
    # All catalogue parameters are 32-bit (or 4-byte) integers or floats, thus
    size_bytes = os.path.getsize(merge_file)
    N_columns = int( (size_bytes-12)/4/N_halos ) # is the number of columns
    catalogue = np.reshape( np.fromfile(in_merge,dtype=np.float32,count=N_columns*N_halos) , (N_halos,N_columns) )

    x,y,z = catalogue[:,:3].T
    m_halo = halo_mass(catalogue[:,6]) # Msol
    del(catalogue) # Free up memory




# Discretization
m = np.logspace( np.log10(m_halo_min*.5), np.log10(m_halo_max*1.05), 100 )
n_m = np.zeros(len(m))
for i in range(N_halos):
    i_upper  = np.searchsorted(m,m_halo[i],side='left')
    i_interp = i_upper - (m[i_upper]-m_halo[i])/(m[i_upper]-m[i_upper-1])
    if i_upper==0:
        n_m[i_upper]+=1
    elif i_upper>0 and i_upper<len(m):
        n_m[i_upper]   += 1-(i_upper-i_interp)
        n_m[i_upper-1] +=    i_upper-i_interp
n_m = n_m/boxsize**3

# Plot n(m)
fig,ax = plt.subplots(nrows=1,ncols=1)
ax.set_yscale('log')
ax.set_xscale('log')

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

ax.plot(m, n_m, label=r'Peak-Patch')
ax.set_xlabel(r'$M/M_{{\odot}}$')
ax.set_ylabel(r'$dn/dM$ $[$Mpc$^{-3}]$')


ax.legend()
fig.set_size_inches(14,6)
fig.savefig(run_dir+'/halo_mass_function.png',dpi=100)



