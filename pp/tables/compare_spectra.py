import numpy as np
import matplotlib.pyplot as plt

file1 = 'planck2018_powerspectrum.dat'
file2 = 'planck18_intermittent.dat'

data1 = np.loadtxt(file1, delimiter=' ')
data2 = np.loadtxt(file2, delimiter=' ')

# Compare P(k)
fig1,ax1 = plt.subplots()
ax1.plot( data1[:,0], data1[:,1], label=file1, ls='-'  )
ax1.plot( data2[:,0], data2[:,1], label=file2, ls='--' )
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel(r'$k$')
ax1.set_ylabel(r'$P(k)$')
ax1.legend()

# Compare T(k)
fig2,ax2 = plt.subplots()
ax2.plot( data1[:,0], data1[:,2], label=file1, ls='-'  )
ax2.plot( data2[:,0], data2[:,2], label=file2, ls='--' )
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlabel(r'$k$')
ax2.set_ylabel(r'$T(k)$')
ax2.legend()

plt.show()
