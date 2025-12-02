import numpy as np

a=np.loadtxt('deltaN-LUT-1.875')
a[:,0]=a[:,0]*2*np.pi*1e5
np.savetxt('deltaN-LUT-1.875_s4000n2000.dat',a,fmt=['%22.15e','%22.15e','%22.15e','%7.9f'])
