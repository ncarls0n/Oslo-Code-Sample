import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.integrate as integrate

# Usage input redshift range:
#    python3 make_a_of_chi.py min max [variable]

# Constants
c = 299792.458 # km/s
s_per_yr = 365.2422*24*3600

# Cosmologcal parameters
h       = 0.6735
Omega_X = 0.2645
Omega_b = 0.0493
Omega_m = Omega_X+Omega_b

# Hubble's constant H_0 = h * (100 km/s/Mpc)
H_0 = h * 100 # km/s/Mpc

# d\chi / dz in Mpc
def dchi_dz(z):
    return c / ( H_0*np.sqrt( Omega_m*(z+1)**3 + 1 - Omega_m ) )

# Read from command line
if len(sys.argv) not in (3,4):
    raise SyntaxError('Incorrect number of command line options.')
elif len(sys.argv)==4:
    variable=str(sys.argv[3])
    if variable in 'xXchiChiCHI':
        variable = 'x'
        xmin     = float(str(sys.argv[1]))
        xmax     = float(str(sys.argv[2]))
else:
    variable = 'z'
    zmin     = float(str(sys.argv[1]))
    zmax     = float(str(sys.argv[2]))

# Redshift range
if variable=='z':
    z   = np.logspace( np.log10(zmin), np.log10(zmax), 1000 )
    chi = np.zeros(len(z))
    for j in range(len(z)):
        chi[j] = integrate.quad( dchi_dz, 0, z[j] )[0]

# Scale factor
a = 1./(z+1.)


# Age as a function of redshift
t = z*0
for j in range(len(z)):
    t[j] = integrate.quad( lambda x: dchi_dz(x)/(1+x), 0, z[j] )[0]
t = t/s_per_yr*1e-9 # in Gyrs

#table = np.fromfile( a_of_chi_file, dtype=np.float32, count=-1 )
#table = np.reshape( table, ( 4, int(len(table)/4) ), order='F' )
#chi = table[0] # in Mpc
#a   = table[1]
#t   = table[2] / s_per_yr * 1e-9 # in Gyrs
#z   = table[3]
#del(table)

def f(q):
    return ('{:.4}'.format(q)).ljust(12)

print( 'index  chi [Mpc]   a           t [Gyr]     z'       )
print( '0      ',f(chi[0]),   f(a[0]),    f(t[0]),    f(z[0])  ,sep='')
print( '-1     ',f(chi[-1]),  f(a[-1]),   f(t[-1]),   f(z[-1]) ,sep='')

c = 299792458 / 1.495978707e11 * np.pi / 648000 / 1e6 * s_per_yr * 1e9 # Mpc/Gyr
index = np.arange(len(a))
def t2chi(t_i):
    i = np.interp(t_i,t,index)
    return np.interp(i,index,chi)
def chi2t(chi_i):
    i = np.interp(chi_i,chi,index)
    return np.interp(i,index,t)

def a2z(a):
    return a**-1-1
def z2a(z):
    return 1/(z+1)


fig,ax = plt.subplots()
ax.plot( chi,a,label='a_of_chi_file' )

# Plot Planck 2018 matter dominated spectrum
def hubble(z):
    return h*100*np.sqrt(Omega_m*(1+z)**3 + 1-Omega_m)
def drdz(z):
    return 299792.458 / hubble(z)
zinterp = np.linspace(0,5,10000)
rinterp = np.cumsum(drdz(zinterp)*(zinterp[1]-zinterp[0]))
rinterp -= rinterp[0]
ax.plot( rinterp, 1/(1+zinterp), ls='--',
    label=r'$Planck$ 2018 $\Lambda-m$-dominated:'+'\n'+r'$z(\chi)=H_0 \int dz \sqrt{ \Omega_{\Lambda,0} + \Omega_{m,0} (1+z)^3 }$'+'\n'+r'$\Omega_{\Lambda,0}=1-\Omega_{m,0}=0.6862 , h=0.6735$' 
    )



ax.set_xlabel(r'comoving distance $\chi$ (Mpc)')
ax.set_ylabel(r'scale factor $a$')

secxax = ax.secondary_xaxis('top', functions=(chi2t,t2chi))
secxax.set_xlabel(r'light-travel time $t$ (Gyr)')
secyax = ax.secondary_yaxis('right', functions=(a2z,z2a))
secyax.set_ylabel(r'redshift $z$')

plt.legend()
plt.show()

#"""
