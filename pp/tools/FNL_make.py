import numpy as np
import pylab as py


sigma = 5.0e-7
Af  = 4.0e-4
w   = 3.#*sigma
piv = 12.#*sigma

chi  = np.linspace(0.0,piv,2000)

zeta = Af * np.exp(-(chi-piv)**2/2/w**2) 

np.savetxt("FNL_spike_w3_piv12.dat",np.transpose([chi,zeta]),fmt='%1.4e')



