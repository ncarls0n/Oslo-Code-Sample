import os
import numpy as np

os.system('rm -rf ../temp')
os.system('cp -r ../instability ../temp')
gotosrc='cd ../temp;'

def kj(k0,kn,j,n):
    return k0*(kn/k0)**(j/n)

m_chi   = 1.
phi_w   = 0.12547
m_tach  = 1.25e3
a_e     = 1.e-40
next    = 4000
boxsize = 4000

# Additional fixed parameters
vev = 0.1       # Vacuum expectation value
phi_p = 8.49953 # Inflaton at extreme of instability
lambda_chi = m_tach*vev**-2
lcode = 2.6259e-52 # Mpc/(10^5 reduced planck masses)
nk   =5#500
k0   = 2*np.pi/a_e/boxsize*lcode
kn   = (next+2)/2*k0
dk   = kj(k0,kn,-1,nk)
kmax = kj(k0,kn,nk+1,nk)

os.system('echo Making chi power spectrum')
s=gotosrc+"sed 's/{0}/{1}/g' {2} > temp;mv temp {2}"
os.system(s.format('M_CHI_REPLACE'     ,m_chi      ,'potential.f90'  ))
os.system(s.format('PHI_W_REPLACE'     ,phi_w      ,'potential.f90'  ))
os.system(s.format('LAMBDA_CHI_REPLACE',lambda_chi ,'potential.f90'  ))
os.system(s.format('PHI_INIT_REPLACE'  ,phi_p-phi_w,'evolve_corr.f90'))
os.system(s.format('PHI_FIN_REPLACE'   ,phi_p+phi_w,'evolve_corr.f90'))
os.system(s.format('NK_REPLACE'        ,nk         ,'evolve_corr.f90'))
os.system(s.format('DK_REPLACE'        ,dk         ,'evolve_corr.f90'))
os.system(s.format('KMAX_REPLACE'      ,kmax       ,'evolve_corr.f90'))

os.system(gotosrc+'module load intel;make clean -f Makefile_corr;make '
    +'-f Makefile_corr;./corr_test')

os.system('echo')
os.system('echo bears')
os.system('echo')

# Read in k and P_{\chi\chi}(k)
src_dir = gotosrc[3:-1]+'/'
p_chichi = np.loadtxt(src_dir+'corr.out',usecols(-17,-6),skiprows=1)
p_chichi[:,0] *= lcode**-1 # set k to units of Mpc^-1

# Save power as Fortran-ordered 32-bit floats in unformatted binary file
p_chichi.T.astype(np.float32).tofile(src_dir+'../tables/p_chichi.dat')
#np.savetxt(src_dir+'/corr.dat',p,fmt='%.5e')

# Clean up correlation matrix 
os.system(gotosrc+'rm -rf corr.out;make clean -f Makefile_corr')


