import numpy as np
import matplotlib.pyplot as plt

fig,ax = plt.subplots(1)

for i in [#'2.5701e6','1.28505e6',
          '1.28505e6','1e6','7.7105e5','5.1402e5','2.5701e5','1.28505e5','6.42525e4','2.5701e4','2.5701e3','0']:
    p = np.loadtxt('Pchi_lambdachi{0}.dat'.format(i))
    ax.plot( p[:,0] , p[:,1], label=r'$\lambda_\chi={0}\times10^{{{1}}}$'.
        format( i[:i.find('e')] , i[i.find('e')+1:] ))

ax.set_xlabel(r'$k$ $[$Mpc$^{-1}]$')
ax.set_ylabel(r'$P_{{\chi\chi}}(k)$ $[$Mpc$^3]$')
ax.set_xscale('log')
ax.set_yscale('log')

plt.legend()
plt.show()
