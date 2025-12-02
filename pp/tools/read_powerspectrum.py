import numpy as np
import pyfftw
import sys,gc
from scipy.special import spherical_jn

###########################################################################
#                                                                         #
#                    READ POWER SPECTRUM FROM 3D FIELD                    #
#                                                                         #
# This script reads the power spectrum from a Peak Patch format initial   #
# conditions field (i.e. files of the form `<run_dir>/fields/Fvec_<...>`, #
# or any of the other possible outputs of the subroutine in hpkvd.f90     #
# RandomField_Output).                                                    #
#                                                                         #
# USAGE (if Peak Patch is in direcotory `<...>/peak-patch/`):             #
#     cd <run_dir>                                                        #
#     python3 <...>/peak-patch/tools/read_powerspectrum.py fields/<field> #
#         n nbuff s                                                       #
#                                                                         #
# where n is the side length of the field in lattice units (including     #
# buffers), nbuff is buffer thickness in lattice units, and s is the side #
# length of the field in Mpc (not including buffers).                     #
#                                                                         #
# Note that as this is a python script, it is slow and should only be     #
# used for small fields. Furthermore, this script is not designed for     #
# parallel runs; Peak Patch IC fields consist of n^3 32-bit (4 byte)      #
# floating-point real numbers, meaning the total RAM load of this script  #
# will be greater than 4*n^3 bytes. It is recomended to use the Fortran   #
# script read_power.f90 instead for larger fiedls, however, the RAM usage #
# will be similar, so a parallelization scheme must be used if the field  #
# requires more memory than a single CPU has.                             #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

filein  = sys.argv[1]
n       = int(sys.argv[2])
nbuff   = int(sys.argv[3])
sbox    = float(sys.argv[4])

if len(sys.argv)>5:
    outfile = sys.argv[5].strip()
    if outfile[-4:] != '.txt':
        outfile=outfile+'.txt'
else:
    outfile = 'powerspectrum.txt'

neff    = n-2*nbuff
dk      = 2*np.pi/sbox
kmax   = neff*dk
nyquist = int(neff/2)+1

# Distance (lattice units) to furthest cell when you fold at the Nyquist
klen  = int(3**.5*(nyquist-1))

# Allocate wavenumber and power spectrum arrays
kbins = np.arange(dk,dk*(klen+2),dk) #np.logspace( np.log10(dk), np.log10(kmax), neff )
P_k   = np.zeros(len(kbins))
P_0   = np.zeros(len(P_k  ))

# Read field - buffers
field = np.zeros((neff,neff,neff))
for i in range(neff):
    f2     = np.array([])
    offset = 4 * (n**2 * (nbuff+i) + n*nbuff + nbuff)
    for j in range(neff):
        f1 = np.fromfile( filein, dtype=np.float32, count=neff, offset = offset+j*4*n )
        f2 = np.concatenate( (f2,f1) )
    field[:,:,i] = np.reshape( f2, (neff,neff), order='F' ) 

# Initialize FFT
f_rfftwf   = pyfftw.empty_aligned( (neff,neff,neff),      dtype='float32'   )
fk_rfftwf  = pyfftw.empty_aligned( (neff,neff,neff//2+1), dtype='complex64' )
rfftwf_obj = pyfftw.FFTW(f_rfftwf, fk_rfftwf, axes=(0,1,2), direction='FFTW_FORWARD')
for i in range(len(f_rfftwf)):
    for j in range(len(f_rfftwf)):
        for k in range(len(f_rfftwf)):
            f_rfftwf[i,j,k] = field[i,j,k]

# Free up memory
del(field);gc.collect()

# Perform FFT, set result to f_k
f_k = rfftwf_obj()


#logkmin = np.log10(dk)
#logkmax = np.log10(kmax)
#dlogk   = (logkmax-logkmin)/neff

for i in range(neff):
    if i<nyquist: ii = i
    else:         ii = nyquist-i

    for j in range(neff):
        if j<nyquist: jj = j
        else:         jj = nyquist-j

        for k in range(neff):
            if k<nyquist: kk = k
            else:         kk = nyquist-k

            p = np.sqrt(ii**2+jj**2+kk**2)
            l = int(np.floor(p))
            c = (1.-np.array([ l-p , l-p+1 ])**2)**2

            P_k[l:l+2] = P_k[l:l+2] + c * np.abs(f_k[i,j,kk])**2
            P_0[l:l+2] = P_0[l:l+2] + c

            #if wavenum==0.:
            #    pk=0
            #else:
            #    d_bin = ( np.log10(wavenum)-(logkmin+(i_bin-1)*dlogk) )/dlogk
            #    pk

# Free up memory
del(f_k);gc.collect()

for i in range(klen+2):
    if P_0[i] != 0.:
        P_k[i] = P_k[i] / P_0[i]
P_k = (neff/sbox)**3*P_k

np.savetxt(outfile, np.array([kbins,P_k]).T)

## Real-space top-hat kernel
#def W(k,R):
#    3*(k*R)**-1*spherical_jn(1,k*R)
#
## Read sigma_8
#def dsigma8_dk(k):
#    return 
#
#np.sum( (2*np.pi**2)**-1 * W(kbins,8)**2 * kbins**2 * P_k )
#
