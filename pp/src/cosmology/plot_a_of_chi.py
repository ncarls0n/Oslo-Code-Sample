import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) == 1:
    a_of_chi_file = 'a_of_chi_table.dat'
elif len(sys.argv) == 2:
    a_of_chi_file = sys.argv[1]
else:
    raise SyntaxError('Incorrect number of command line options.')

s_per_yr = 365.2422*24*3600

table = np.fromfile( a_of_chi_file, dtype=np.float32, count=-1 )
table = np.reshape( table, ( 4, int(len(table)/4) ), order='F' )

chi = table[0] # in Mpc
a   = table[1]
t   = table[2] / s_per_yr * 1e-9 # in Gyrs
z   = table[3]

del(table)

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



ax.plot( chi,a,label=a_of_chi_file )

# Plot Planck 2018 matter dominated spectrum
h       = .6735
omega_m = .2645+.0493
def hubble(z):
    return h*100*np.sqrt(omega_m*(1+z)**3 + 1-omega_m)
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
