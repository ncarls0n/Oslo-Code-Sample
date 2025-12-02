import numpy as np
from os.path import exists

dats = [
  'Pchi_lambdachi2.5701e6.dat',
  'Pchi_lambdachi2.5701e7.dat',
  'Pchi_lambdachi5.1402e5.dat',
  'Pchi_lambdachi6.42525e4.dat',
  'Pchi_lambdachi7.7105e5.dat',
  'Pchi_lambdachi1.28505e5.dat',
  'Pchi_lambdachi1.28505e6.dat',
  'Pchi_lambdachi1.6e3.dat',
  'Pchi_lambdachi1.6e4.dat',
  'Pchi_lambdachi1.6e5.dat',
  'Pchi_lambdachi1.6e6.dat',
  'Pchi_lambdachi1.6e7.dat',
  'Pchi_lambdachi2.5701e3.dat',
  'Pchi_lambdachi2.5701e4.dat',
  'Pchi_lambdachi2.5701e5.dat',
  'Pchi_lambdachi2.5701e5_nk1000.dat',
  'Pchi_lambdachi2.5701e5_nk10000.dat',
  'Pchi_lambdachi2.5701e5_nk300.dat',
  'Pchi_lambdachi2.5701e5_nk5000.dat'  ]

for file_in in dats:
    if not exists(file_in):
        continue
    matrix = np.loadtxt(file_in)
    np.savetxt(file_in,matrix,fmt='%.4e')
