import numpy as np
import matplotlib.pyplot as plt

alpha = 1.076104949214213
a     = np.exp(alpha)
H     = 3.434759510783179

p_light = np.loadtxt('pchi_dk.1H_n500_lambdachi0_m1e-2.dat',usecols=(-17,-6))
p_int   = np.loadtxt('pchi_dk.1H_n500_lambdachi0_m1.dat',   usecols=(-17,-6))
p_mid   = np.loadtxt('pchi_dk.1H_n500_lambdachi0_m10.dat',  usecols=(-17,-6))
p_welt  = np.loadtxt('pchi_dk.1H_n500_lambdachi0_m100.dat', usecols=(-17,-6))
p_heavy = np.loadtxt('pchi_dk.1H_n500_lambdachi0_m1e3.dat', usecols=(-17,-6))

fig,ax = plt.subplots(1)
ax.plot( p_light[:,0]/a/H , (p_light[:,0]/a/H)**3 / 2*np.pi**2 *p_light[:,1] , label=r'$m_\chi^2 \sim \frac{H^2}{1000}$' )
ax.plot( p_int[:,0]/a/H , (p_int[:,0]/a/H)**3 / 2*np.pi**2 *p_int[:,1] , label=r'$m_\chi^2 \sim \frac{H^2}{10}$' )
ax.plot( p_mid[:,0]/a/H , (p_mid[:,0]/a/H)**3 / 2*np.pi**2 *p_mid[:,1] , label=r'$m_\chi^2 \sim H^2$')
ax.plot( p_welt[:,0]/a/H , (p_welt[:,0]/a/H)**3 / 2*np.pi**2 *p_welt[:,1] , label=r'$m_\chi^2 \sim 10H^2$')
ax.plot( p_heavy[:,0]/a/H , (p_heavy[:,0]/a/H)**3 / 2*np.pi**2 *p_heavy[:,1] , label=r'$m_\chi^2 \sim 100 H^2$' )

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel(r'$\frac{k}{aH}$')
ax.set_ylabel(r'$\mathcal{P}_{\chi\chi}$')

plt.legend()
plt.show()
