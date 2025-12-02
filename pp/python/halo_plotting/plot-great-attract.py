# plot-great-attract.py
#
# Nathan J. Carlson
# September 9, 2021
#
# Optimized to run on Python 3.8.2
#
# 
# This script locates "Great Attractors", wells in the primordial gravi-
# tational potential scalar field Phi that draw in matter from surrounding
# regions of space. Great Attractors are located using two methods: (1)
# minima of the initial Phi field smoothed at various filter scales, and
# (2) using theoretical methods to extrapolate its position from shear
# tensor for each peak. Peak patches are located and displaced using linear
# and quadratic Lagrangian perturbations (1LPT & 2LPT), so the great at-
# tractor model acts as a third correction.

import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import signal
from mayavi import mlab # Mayavi free 3D plotting module

# RUN INSTRUCTIONS
# 
# This script is run in terminal. To do so, you need to pass it as
# arguments the sidelength of cubic simulation volume as a float
# `<side_length>`, the thickness (in number of cells) of the cubic shell
# that makes up the buffer in the density field `<buffer_size>` as an
# integer the merged peak patch output file `<merged_peakfile>` (the final
# result of a peak patch run), and the initial density field
# `<density_field>`. The initial density field is output to a directory
# `fields` when the parameter `ioutfield` is set to 1 or 2 in the
# parameter file.

# The terminal command is thus
# 
#     cd <path_to_peak-patch>
#     python3.8 <py_dir>/plot-great-attract.py <run_dir>
# 
# For example, this might look like
# 
#     cd <path_to_peak-patch>
#     python3.8 python/halo_plotting/plot-great-attract.py runs/ns128

if len(sys.argv) != 2:
    # Error message displayed if you pass the wrong number of arguments.
    print('You\'ve passed me the wrong number of arguments at the command '
        +'line.\nCommand line should take the form\ncd <peak_patch_directo'
        +'ry>\npython3.8 python/halo_plotting/plot-great-attract.py runs/n'
        +'s128')
    sys.exit(2)

# Reading from command line prompts
run_dir = str(sys.argv[1])   # The directory of the run


###########################################################################
### Execute useful lines from parameter file                            ###
###########################################################################

# Open parameter file
with open(run_dir+'/param/param.params') as f:
    params = [i.strip() for i in f.readlines()]

# Define variables
# seed               initial random field seed number (an integer)
# boxsize [h^-1 Mpc] sidelength of simulation volume
# nmesh              number of cells in one tile (including buffers)
# nbuff              number of cells in buffer
# ntile              number
# sigma8             stdev of Gaussian random field smoothed at 8 h^-1 Mpc
# Omx                Omega_CDM, energy density fraction of CDM
# OmB                Omega_B, energy density fraction of baryonic matter
# Omvac              Omega_Lambda, energy denstiy fraction of DE
# h                  "little h", H_0/100 km/s/Mpc
# ns                 n_s, spectral index
# run_name           string used in output files
# short_name         shorter string used in output files
# maximum_redshift   z_max, redshift of primordial fields
# global_redshift    z_0, redshift of Eulerian peaks
for line in params:
    if ( line[:4]=='seed'  or line[:7]=='boxsize' or line[:5]=='nmesh'  or
         line[:5]=='nbuff' or line[:5]=='ntile'   or line[:6]=='sigma8' or
         line[:3]=='Omx'   or line[:3]=='OmB'     or line[:5]=='Omvac'  or
         line[:2]=='h '    or line[:2]=='ns'      or
         line[:8]=='run_name'          or line[:10]=='short_name'       or
         line[:16]=='maximum_redshift' or line[:15]=='global_redshift'  ):
       exec(line)

# Computes the LambdaCDM scale factor a(t)
scale_factor = (1 + maximum_redshift - global_redshift)**-1

# Additional parameters
nlattice = int( (nmesh-2*nbuff) * ntile + 2*nbuff )
neff     = int( nlattice - 2*nbuff )
Omm      = Omx+OmB        # Omega_m, fraction of energy that can cluster
rhocrit  = 2.775e11*h**2  # critical energy density 3H^2/8piG [Msol Mpc^-3]
rho      = rhocrit*Omm    # average CDM density [Msol Mpc^-3]
deltavir = 200            # <something defining virial collapse I think>
outnum   = 33             # number of columns in <merged_peak_file>
G        = 4.517e-48      # gravitational constant [Mpc^3 Msol^-1 s^-2]
# Note here that we use mass units of solar masses Msol, spatial units of
# megaparsecs Mpc, and time units of seconds s.


###########################################################################
### Reading in field & peak catalogue data                              ###
###########################################################################

# Pointer to merged peak patch catalogue
in_catalogue = open( run_dir + '/output/' + run_name + '_merge.pksc.' +
                     str(seed) , 'rb' )

# Pointer to primordial overdensity field
in_delta     = open( run_dir + '/fields/Fvec_' + run_name , 'rb' )

# Reading data from merged peak patch catalogue pointer
Non       = np.fromfile(in_catalogue,dtype=np.int32,count=1)[0]
RTHmax    = np.fromfile(in_catalogue,dtype=np.float32,count=1)[0]
zin       = np.fromfile(in_catalogue,dtype=np.float32,count=1)[0]
catalogue = np.fromfile(in_catalogue,dtype=np.float32,count=outnum*Non)
catalogue = np.reshape(catalogue,(Non,outnum))
# The merged peak catalogue has three header parameters: `Non`, the number
# of halos found, `RTHmax`, the largest top-hat filter at which peaks were
# found, `zin` initial redshift. The remainder of the file contains a 1D 
# array that we then reshape into an XY-column matrix `peakdata`. Each row
# of this XY-column matrix contains data describing a unique halo.

# Reading data from primordial overdensity field pointer
delta = np.fromfile(in_delta,dtype=np.float32,count=-1)
delta = np.reshape(delta, (nlattice,nlattice,nlattice), order='F')
delta = delta[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]

# Defines X,Y,Z as meshgrid
s,n = boxsize,neff
edges = np.linspace( -s/2 , s/2 , n+1 )
X,Y,Z = np.meshgrid(edges,edges,edges,indexing='ij')

# Next, for the sake of clarity, we copy the raw halo data from the XY
# columns of `peakdata` to XY better-labelled NumPy arrays:
x      = catalogue[:,0]  # x,y,z: components of the final (Eulerian) halo
y      = catalogue[:,1]  #     position vector [h^-1 Mpc]
z      = catalogue[:,2]  # vx,vy,vz: components of the final (Eulerian) halo
vx     = catalogue[:,3]  #     velocity vector [km/s/Mpc]
vy     = catalogue[:,4]  # Rth: top-hat filter radius halo was found at
vz     = catalogue[:,5]  #     (roughly the radius of the halo) [h^-1 Mpc]
Rth    = catalogue[:,6]  # xL,yL,zL: components of the initial (Lagrangian)
xL     = catalogue[:,7]  #     halo position vector [h^-1 Mpc]
yL     = catalogue[:,8]  # Fpk: linearly extrapolated initial overdensity
zL     = catalogue[:,9]  #     field at peak centre, also called Fcollv in
Fpk    = catalogue[:,10] #     fortran codes
e_v    = catalogue[:,11] # e_v: ellipticity, e_v=(lambda_3-lambda1)/2Fpk
p_v    = catalogue[:,12] # p_v: prolaticity, p_v=1/2-3lambda_2/2Fpk
epsxx  = catalogue[:,13] #     where e_v >= 0, -e_v <= p_v <= e_v
epsyy  = catalogue[:,14] # eps^ij: the i,jth component of the symmetric
epszz  = catalogue[:,15] #     strain tensor at the corresponding peak
epsyz  = catalogue[:,16] #     
epsxz  = catalogue[:,17] #     
epsxy  = catalogue[:,18] # 
#F_d2   = catalogue[:,19]  # 
#zform   = catalogue[:,20] # This one =-1 for most peaks for some reason
#gradx   = catalogue[:,21] #
#grady   = catalogue[:,22] # Don't need this stuff yet. Will have to dig
#gradz   = catalogue[:,23] # through hpkvd.f90 etc to figure out what any
#gradfx  = catalogue[:,24] # of these actually means.
#gradfy  = catalogue[:,25] #
#gradfz  = catalogue[:,26] #
#Rfclv   = catalogue[:,27] #
#FcollvRf= catalogue[:,28] #
#F_d2Rf  = catalogue[:,29] #
#gradrfx = catalogue[:,30] #
#gradrfy = catalogue[:,31] #
#gradrfz = catalogue[:,32] #

# Let's free up some hard drive space
del(catalogue)

print('The extent of the filter radii used was ',np.min(Rth), np.max(Rth) )

# Eigenvalues, lambda3 >= lambda2 >= lambda1
lambda_3 = Fpk/3 * (1+3*e_v+p_v) #
lambda_2 = Fpk/3 * (1-2*p_v)     # Fpk = sum_i lambda_i = -Tr(eps)
lambda_1 = Fpk/3 * (1-3*e_v+p_v) #

# The anisotropic part of the eigenvalues lambda`_i
lambdan3 = lambda_3 - Fpk/3
lambdan2 = lambda_2 - Fpk/3
lambdan1 = lambda_1 - Fpk/3

# Semi-axes of the initial Lagrangian ellipsoids
R_1 = Rth/2 * np.exp( -lambdan1 )
R_2 = Rth/2 * np.exp( -lambdan2 )
R_3 = Rth/2 * np.exp( -lambdan3 )
# Where R_3 <= R_2 <= R_1, R_1 is the last axis to collapse.


# Read in 1LPT parameters
params1LPT = run_dir+'/output/1LPT/'+run_name+'_merge.pksc.'+str(seed)
from read1LPT import get1LPTparams
( x1LPT  , y1LPT  , z1LPT, # 1LPT Eulerian peak positions
  vx1LPT , vy1LPT , vz1LPT # 1LPT Eulerian peak velocities
) = get1LPTparams(params1LPT)


# Finally, we can determine the approximate mass of each halo based on its
# radius (which is greater than or equal to the largest top-hat radius at
# which it is found to have collased `Rth`) as well as the critical density
# of cold dark matter `rho`:
M = 4*np.pi/3*Rth**3*rho # halo masses [Msol]


###########################################################################
### Theoretical great attractor locations                               ###
###########################################################################

# 1LPT peak position (Eulerian space)
r_pk1LPT = np.array([ [x1LPT[i],y1LPT[i],z1LPT[i]] for i in range(Non) ])

# Peak position (Eulerian space)
r_pk = np.array([ [ x[i] , y[i] , z[i] ] for i in range(Non) ])

# Peak position (Lagrangian pace)
q_pk = np.array([ [ xL[i] , yL[i] , zL[i] ] for i in range(Non) ])

# Peak displacement
s_pk = r_pk - q_pk

# Strain tensor defined as eps^ij = (ds^i/dq^j+ds^j/dq^i)/2
eps = -np.array([
    [[ epsxx[i] , epsxy[i] , epsxz[i] ],
     [ epsxy[i] , epsyy[i] , epsyz[i] ],
     [ epsxz[i] , epsyz[i] , epszz[i] ]] for i in range(Non) ])

# # Anisotropic strain tensor
# eps_aniso = eps + np.array([
#     [[ Fpk[i]/3 ,       0. ,       0. ],
#      [       0. , Fpk[i]/3 ,       0. ],
#       [       0. ,       0. , Fpk[i]/3 ]] for i in range(Non) ])

# Inverse strain tensor
epsinv = np.array([ np.linalg.inv(eps[i]) for i in range(Non) ])

# Lagrangian Great Attractor position
q_GA = q_pk - np.array([ epsinv[i].dot(s_pk[i]) for i in range(Non) ])

# Check q_GA is actually an attractor
sb_lt_x = s_pk + np.array([ eps[i].dot(q_GA[i] -
    np.array([boxsize/neff/2,0,0])-q_pk[i]) for i in range(Non) ])
sb_gt_x = s_pk + np.array([ eps[i].dot(q_GA[i] +
    np.array([boxsize/neff/2,0,0])-q_pk[i]) for i in range(Non) ])
sb_lt_y = s_pk + np.array([ eps[i].dot(q_GA[i] -
    np.array([0,boxsize/neff/2,0])-q_pk[i]) for i in range(Non) ])
sb_gt_y = s_pk + np.array([ eps[i].dot(q_GA[i] +
    np.array([0,boxsize/neff/2,0])-q_pk[i]) for i in range(Non) ])
sb_lt_z = s_pk + np.array([ eps[i].dot(q_GA[i] -
    np.array([0,0,boxsize/neff/2])-q_pk[i]) for i in range(Non) ])
sb_gt_z = s_pk + np.array([ eps[i].dot(q_GA[i] +
    np.array([0,0,boxsize/neff/2])-q_pk[i]) for i in range(Non) ])

flags = np.intersect1d(
    np.intersect1d(
        np.intersect1d(
            np.argwhere( sb_lt_x[:,0]>0 ), np.argwhere( sb_gt_x[:,0]<0 )),
        np.intersect1d(
            np.argwhere( sb_lt_y[:,1]>0 ), np.argwhere( sb_gt_y[:,1]<0 ))),
    np.intersect1d(
        np.argwhere( sb_lt_z[:,2]>0 ), np.argwhere( sb_gt_z[:,2]<0 )))


###########################################################################
### Eigenvectors of strain tensor                                       ###
###########################################################################

# Find eigenvalues lambda_pk and eigenvectors eig_v_pk of strain tensor
eigs_pk   = [ np.linalg.eig(-eps[i]) for i in range(Non) ]
lambda_pk = np.array([ eigs_pk[i][0] for i in range(Non) ])
eigvec_pk = np.array([ eigs_pk[i][1] for i in range(Non) ]); del(eigs_pk)

# Arrange so that lambda_1 < lambda_2 < lambda_3
for i in range(Non):
    args = np.argsort(lambda_pk[i])
    lambda_pk[i] = lambda_pk[i,args]
    eigvec_pk[i] = eigvec_pk[i,args]

# Rotation matrix i.e. P eps P^-1 = diag(lambda_1,lambda_2,lambda_3)
P = np.array([ np.array( [eigvec_pk[i,0],eigvec_pk[i,1],eigvec_pk[i,2]] ).T
               for i in range(Non) ])

# Euler angles for P = R_zxz
beta  = np.array([ np.arccos(P[i,2,2]) for i in range(Non) ])
alpha = np.array([ np.arcsin(P[i,0,2]/np.sin(beta[i])) % (2*np.pi)
                   for i in range(Non) ])
gamma = np.array([ np.arcsin(P[i,2,0]/np.sin(beta[i])) % (2*np.pi)
                   for i in range(Non) ])

'''
for i in range(40,43):
    print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _')
    print('alpha=',round(alpha[i],2),',  beta=',round(beta[i],2),',  gamma=',round(gamma[i],2))
    print()
    print(P[i])
    print()
    rotmatx = np.array(
        [[ np.cos(alpha[i])*np.cos(gamma[i]) - np.sin(alpha[i])*np.cos(beta[i])*np.sin(gamma[i])  ,
          -np.cos(alpha[i])*np.sin(gamma[i]) - np.sin(alpha[i])*np.cos(beta[i])*np.cos(gamma[i])  ,
           np.sin(alpha[i])*np.sin(beta[i])] ,
         [ np.sin(alpha[i])*np.cos(gamma[i]) + np.cos(alpha[i])*np.cos(beta[i])*np.sin(gamma[i])  ,
          -np.sin(alpha[i])*np.sin(gamma[i]) + np.cos(alpha[i])*np.cos(beta[i])*np.cos(gamma[i])  ,
          -np.cos(alpha[i])*np.sin(beta[i])] ,
         [ np.sin(beta[i])*np.sin(gamma[i])  , np.sin(beta[i])*np.cos(gamma[i]) , np.cos(beta[i]) ]])
    print( rotmatx )
    print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _')
'''
# # Inverse of P and the diagonalized eps matrix
# P_inv = np.array([ np.linalg.inv(P[i]) for i in range(Non) ])
# D = np.array([ (P[i].dot(eps[i])).dot(P[i].T) for i in range(Non) ])

# Unit eigenvectors
unit_eigvec_pk = np.array([ [
    eigvec_pk[i,0]/np.linalg.norm(eigvec_pk[i,0]) ,
    eigvec_pk[i,1]/np.linalg.norm(eigvec_pk[i,1]) ,
    eigvec_pk[i,2]/np.linalg.norm(eigvec_pk[i,2]) ] for i in range(Non) ])
# unit_eigvec_pk[:,0] points in the direction of the major axis
# unit_eigvec_pk[:,1] points in the direction of the intermediate axis
# unit_eigvec_pk[:,2] points in the direction of the minor axis

# Ellipsoid semi-axes
semiax_pk = np.array([
    [ np.exp( -lambda_1[i]+Fpk[i]/3 ) * Rth[i] ,
      np.exp( -lambda_2[i]+Fpk[i]/3 ) * Rth[i] ,
      np.exp( -lambda_3[i]+Fpk[i]/3 ) * Rth[i] ] for i in range(Non) ])
# semiax_pk[:,0] is the semimajor axis
# semiax_pk[:,1] is the semiintermediate axis
# semiax_pk[:,2] is the semiminor axis

#print(lambda_1[30],    lambda_2[30],    lambda_3[30]   )
#print(lambda_pk[30,0], lambda_pk[30,1], lambda_pk[30,2])

# Find peak displacement in terms of principle axes
s_pk_pr = np.array([
    [ s_pk[i].dot(unit_eigvec_pk[i,0]) ,
      s_pk[i].dot(unit_eigvec_pk[i,1]) ,
      s_pk[i].dot(unit_eigvec_pk[i,2]) ] for i in range(Non) ])


###########################################################################
### Fourier transforms                                                  ###
###########################################################################

# Fourier phase k [h Mpc^-1]
kscale     = scipy.fft.fftfreq(neff,boxsize/neff)
k_for_plot = scipy.fft.fftshift(kscale) #=np.linspace(-n/s/2, (n-2)/s/2, n)
kscale[np.argwhere(kscale==0)] = np.inf
kX,kY,kZ   = np.meshgrid(kscale,kscale,kscale,indexing='ij')
# The k=(0,0,0) phase component of deltak corresponds to an additive term
# in delta, and the DFT is invariant under additive factors, so we can just
# set this term of the sum over k to zero. Since we take
# iFFT(k^-2 FFT(delta)), this is accomplished by replacing k=0 with k=inf.

# Fourier space overdensity field
deltak          = scipy.fft.fftn(delta)               # FFT(delta)
# deltak_for_plot = np.real(scipy.fft.fftshift(deltak)) # shifts domain
# inverse_deltak  = np.real(scipy.fft.ifftn(deltak))    # iFFT(FFT(delta))
''' For reasons related to symmetries of phase space, the scipy function
    fftfreq places 0 as the first entry of the an array, followed by the
    positive frequencies, followed by the negative frequencies:
    k = fftfreq(n,s/n) = [ 0.   ,  1/s      ,  2/s ,  ...  ,  (n-1)/2s, 
                           -n/2s,  -(n-1)/2s,  ... ,   -2/s,  -1/s    ])
    this isn't very conveninet if we want to plot it, so for the array
    deltak_for_plot, I use the scipy function fftshift to rearange this
    numpy array into a more intuitive form:
    fftshift(k) = [ -n/2s,  -(n-1)/2s,  ...,  -2/s,  -1/s    ,
                    0.   ,  1/s      ,  2/s,  ... ,  (n-1)/2s])         '''

# Gravitational potential field Phi [Mpc^2 s^-2]
Phi = scipy.fft.ifftn( deltak / (kX**2+kY**2+kZ**2))
Phi = -4*np.pi * G * rho * np.real(Phi)
# In Lagrangian space, we drop the scale factor a(t)=1

# Restore kX,kY,kZ
kscale[np.argwhere(kscale==np.inf)] = 0.
kX,kY,kZ = np.meshgrid(kscale,kscale,kscale,indexing='ij')

# # Fourier transform of Phi
# Phik = scipy.fft.fftn(Phi)
# 
# # delta from Phi0^-1 * iFFT(-k^2 FFT(Phik))
# delta_from_Phi = ( (4*np.pi*G*rho)**(-1) *
#     np.real(scipy.fft.ifftn( -(kX**2+kY**2+kZ**2) * Phik )) )
# 
# Gradient of Phi given by iFFT(k FFT(Phik)) where k is a vector
# dPhi = np.array([ -1.j * np.real( scipy.fft.ifftn(kX*Phik) ) ,
#                   -1.j * np.real( scipy.fft.ifftn(kY*Phik) ) ,
#                   -1.j * np.real( scipy.fft.ifftn(kZ*Phik) ) ])
# dPhi_mag = np.sqrt(np.abs( dPhi[0]**2 + dPhi[1]**2 + dPhi[2]**2 ))


###########################################################################
### Smoothers                                                           ###
###########################################################################

def W_G(qX,qY,qZ,R_G):
    return (2*np.pi)**-1.5 * R_G**-3 * np.exp(
        -(qX**2+qY**2+qZ**2) / (2*R_G**2) )

def W_th(qX,qY,qZ,R_th):
    return 3/(4*np.pi) * R_th**-3 * np.heaviside(
        R_th - np.sqrt(qX**2+qY**2+qZ**2) , 0 )


###########################################################################
### Locating extrema of primordial fields                               ###
###########################################################################

def minima_and_maxima(Phi, n, R=1, out='all'):
    '''
    '''
    # Gradient of Phi field
    dPhi = np.gradient( Phi )
    dPhi_mag = np.sqrt( dPhi[0]**2 + dPhi[1]**2 + dPhi[2]**2 )
    
    Phi_local_minima   = np.array([ [], [], [] ])
    Phi_local_min_mags = np.array([])
    Phi_local_maxima   = np.array([ [], [], [] ])
    Phi_local_max_mags = np.array([])
    dPhi_threshold   = np.min(dPhi_mag) + np.max(dPhi_mag)/5

    for a in range(int(R/2),n,R):
        for b in range(int(R/2),n,R):
            for c in range(int(R/2),n,R):
                if dPhi_mag[a,b,c] < dPhi_threshold:
                    if ( Phi[a,b,c] > Phi[ (a+R)%n, b,       c       ]and
                         Phi[a,b,c] > Phi[ a-R,     b,       c       ]and
                         Phi[a,b,c] > Phi[ a,       (b+R)%n, c       ]and
                         Phi[a,b,c] > Phi[ a,       b-R,     c       ]and
                         Phi[a,b,c] > Phi[ a,       b,       (c+R)%n ]and
                         Phi[a,b,c] > Phi[ a,       b,       c-R     ]):

                        Phi_local_maxima = np.array([
                            np.append( Phi_local_maxima[0],a ),
                            np.append( Phi_local_maxima[1],b ),
                            np.append( Phi_local_maxima[2],c )])
                        Phi_local_max_mags = np.append(
                            Phi_local_max_mags, Phi[a,b,c] )

                    elif ( Phi[a,b,c] < Phi[ (a+R)%n, b,       c       ]and
                           Phi[a,b,c] < Phi[ a-R,     b,       c       ]and
                           Phi[a,b,c] < Phi[ a,       (b+R)%n, c       ]and
                           Phi[a,b,c] < Phi[ a,       b-R,     c       ]and
                           Phi[a,b,c] < Phi[ a,       b,       (c+R)%n ]and
                           Phi[a,b,c] < Phi[ a,       b,       c-R     ]):

                        Phi_local_minima = np.array([
                            np.append( Phi_local_minima[0],a ),
                            np.append( Phi_local_minima[1],b ),
                            np.append( Phi_local_minima[2],c )])
                        Phi_local_min_mags = np.append(
                            Phi_local_min_mags, Phi[a,b,c] )
    if out=='min':
        return Phi_local_minima
    elif out=='minmag':
        return Phi_local_minima,Phi_local_min_mags
    elif out=='max':
        return Phi_local_maxima
    elif out=='maxmag':
        return Phi_local_maxima,Phi_local_max_mags
    else:
        return( Phi_local_minima , Phi_local_min_mags ,
                Phi_local_maxima , Phi_local_max_mags )


###########################################################################
### Plotting functions                                                  ###
###########################################################################

def plot_density_hist(field, X,Y, label1, label2, show=True):
    """ field, X and Y are NumPy arrays of dimension n,n,n representing
        respectively a density field and the x and y components of a
        meshgrid describing it's dependent variables.
        label1 and label2 are lists of strings.

        Usage:
        fig0,axs0=plot_density_hist( delta, X,Y,
            (r'a)', r'$x$ [$h^{-1}$Mpc]', r'$y$ [$h^{-1}$Mpc]', r'$\delta$ or something'),
            (r'b)', r'$\bar{\rho}\delta = \rho-\bar{\rho}$', r'counts per overdensity bin'),
            show=False)
    """    
    fig,axs = plt.subplots(nrows=1,ncols=2)

    # Plotting slice of overdensity field
    ax = axs[0].pcolormesh( X[:,:,0], Y[:,:,0], field[:,:,25], cmap='viridis')
    axs[0].set_title(  label1[0] )
    axs[0].set_xlabel( label1[1] )
    axs[0].set_ylabel( label1[2] )
    ax_pos = make_axes_locatable(axs[0])
    ax_colorbar_pos = ax_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(ax, cax=ax_colorbar_pos)
    ax_colorbar_pos.yaxis.set_label_position('right')
    ax_colorbar_pos.set_ylabel( label1[3], labelpad=15, rotation=270)
    axs[0].set_aspect(1)

    # Plotting overdensity histogram
    bins = int( (np.max(field) - np.min(field)) *      # Number of bins using
        len(field) / (2*scipy.stats.iqr(field)) + .5 ) # Freedman-Draconis rule
    counts,bin_edges = np.histogram(field, bins=bins)
    bin_centres      = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2
    axs[1].hist( bin_edges[:-1], bin_edges, weights=counts )
    axs[1].set_title(  label2[0] )
    axs[1].set_xlabel( label2[1] )
    axs[1].set_ylabel( label2[2] )

    if show==True:
        plt.show()

    else:
        return fig,axs


def plot_great_attractors( q_pk, r_pk, r_pk1LPT, Rth, Phi, PhiL ,q_GA,
        s, n, xspan=[0,neff], yspan=[0,neff], zspan=[neff//2,neff//2+1],
        R_label=r'$R_{lattice}$', display=False):
    ''' INPUTS/OUTPUTS:
        'PhiL' is a 4D NumPy array containing the positions of each of the
        local maxima it is of the same form as 'Phi_local_maxima' and 
        'Phi_local_minima' defined in the main body of the code above. The
        input parameters 'xspan', 'yspan', and 'zspan' are lists of two
        integers and 'display' is a boolian parameter. If display=True,
        there is no output, and the function simply displays the figure
        generated. If display=False, then 'fig' (of the class
        'matplotlib.figure.Figure') is returned. The remaining input
        parameters share the class of like-named variables in the main code
        body. 
        
        DESCRIPTION:
        This function plots the Phi peaks (found using a sequential
        nearest-neighbours search) and peak patches that lie within the
        spatial bounds laid out by 'xspan', 'yspan' and 'zspan' super-
        imposed on the Phi slice defined by the middle of this region (i.e.
        if zspan[1]-zspan[0]<xspan[1]-xspan[0],yspan[1]-yspan[0] then the
        slice at Phi[x,y,z=(zspan[1]-zspan[0])/2] is plotted).

        USAGE:
        This function can be used to directly display a plot (as described
        above) to the screen, e.g. with
        >>> plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,PhiL=Phi_local_maxima)
        
        Or to return a matplotlib.figure.Figure class object, which allows
        us to display multiple figures at once, e.g. with
        >>> figure4 = plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,
            PhiL=Phi_local_maxima,xspan=[0,n],yspan=[0,n],
            zspan=[n//2,n//2+1],display=False)
        >>> plt.show()
    '''
    from matplotlib.lines import Line2D

    # Convert Phi from Mpc^2 s^-2 to km^2 s^-2
    Phi = 9.5214061369e38*Phi

    # Discarding peaks that lie outside of the slice
    w = np.intersect1d(
          np.intersect1d(
            np.intersect1d( np.argwhere( q_pk[:,0]>=-s/2+s/n*xspan[0] ) ,
                            np.argwhere( q_pk[:,0]< -s/2+s/n*xspan[1] ) ),
            np.intersect1d( np.argwhere( q_pk[:,1]>=-s/2+s/n*yspan[0] ) ,
                            np.argwhere( q_pk[:,1]< -s/2+s/n*yspan[1] ) )),
          np.intersect1d( np.argwhere( q_pk[:,2]>=-s/2+s/n*zspan[0] ),
                          np.argwhere( q_pk[:,2]< -s/2+s/n*zspan[1] )))
    q_pk     = q_pk[w]
    r_pk     = r_pk[w]
    q_GA     = q_GA[w]
    Rth      = Rth[w].reshape(len(w))
    r_pk1LPT = r_pk1LPT[w]
    
    # Discarding Phi minima that lie outside the slice
    w = np.intersect1d(
            np.intersect1d(
                np.intersect1d( np.argwhere( PhiL[0,:] >= xspan[0] ),
                                np.argwhere( PhiL[0,:] <  xspan[1] )),
                np.intersect1d( np.argwhere( PhiL[1,:] >= yspan[0] ),
                                np.argwhere( PhiL[1,:] <  yspan[1] ))),
            np.intersect1d( np.argwhere( PhiL[2,:] >= zspan[0] ),
                            np.argwhere( PhiL[2,:] <  zspan[1] )))
    PhiL = -s/2 + s/n * PhiL[:,w]
    
    # Dependent variable used to plot circular peak patches
    theta = np.linspace(0, 2*np.pi, 100, endpoint=True)
    
    # Define bin edges for Phi plot
    edges = np.linspace( -s/2 , s/2 , n+1 )
    X,Y,Z = np.meshgrid(edges,edges,edges,indexing='ij')
 
    # Creates figure
    fig,axs = plt.subplots(1, figsize=(12,10))

    # Determines what plane the volume slice is in and defines figure
    # labels accordingly
    if ( zspan[1]-zspan[0]<=xspan[1]-xspan[0] and
         zspan[1]-zspan[0]<=yspan[1]-yspan[0] ):
        the_plane  = r'$(x,y)$'
        Phi_cells  = (zspan[0]+zspan[1]) // 2
        Phi_slice  = -s/2+s/n*Phi_cells
        the_xlim   = [-s/2+s/n*xspan[0],-s/2+s/n*xspan[1]]
        the_ylim   = [-s/2+s/n*yspan[0],-s/2+s/n*yspan[1]]
        the_zlim   = [-s/2+s/n*zspan[0],-s/2+s/n*zspan[1]]
        #the_title  = (r'$\Phi(x,y,z=$'+str(Phi_slice)+
        #    r'$~h^{-1}\mathrm{Mpc}),~\mathrm{peaks~at}~z\in[$'+
        #    str(round(the_zlim[0],2))+r'$~h^{-1}\mathrm{Mpc},$'+
        #    str(round(the_zlim[1],2))+r'$~h^{-1}\mathrm{Mpc})$')
        the_xlabel = r'$x ~ [h^{-1} \mathrm{Mpc}]$'
        the_ylabel = r'$y ~ [h^{-1} \mathrm{Mpc}]$'
        cbar_label = (r'$\Phi(x,y,z=$'+str(round(Phi_slice,2))+
            r'$~h^{-1}\mathrm{Mpc})~[\mathrm{km}^2\,\mathrm{s}^{-2}]$')
        legend_string = ( r'$z\in[$'+str(round(the_zlim[0],2))+
            r'$~h^{-1}\mathrm{Mpc},$'+str(round(the_zlim[1],2))+
            r'$~h^{-1}\mathrm{Mpc})$')

        # Plots corresponding slice of Phi
        ax = axs.pcolormesh( X[:,:,0], Y[:,:,0], Phi[:,:,Phi_cells],
                             cmap='viridis')

    elif ( xspan[1]-xspan[0]<=yspan[1]-yspan[0] and
           xspan[1]-xspan[0]<=zspan[1]-zspan[0] ):
        the_plane  = 'r$(y,z)$'
        Phi_cells  = (xspan[0]+xspan[1]) // 2
        Phi_slice  = -s/2+s/n*Phi_cells
        the_xlim   = [-s/2+s/n*yspan[0],-s/2+s/n*yspan[1]]
        the_ylim   = [-s/2+s/n*zspan[0],-s/2+s/n*zspan[1]]
        the_zlim   = [-s/2+s/n*xspan[0],-s/2+s/n*xspan[1]]
        #the_title  = (r'$\Phi(x=$'+str(Phi_slice)+
        #    r'$~h^{-1}\mathrm{Mpc},y,z),~\mathrm{peaks~at}~x\in[$'+
        #    str(round(the_zlim[0],2))+r'$~h^{-1}\mathrm{Mpc},$'+
        #    str(round(the_zlim[1],2))+r'$~h^{-1}\mathrm{Mpc})$')
        the_xlabel = r'$y ~ [h^{-1} \mathrm{Mpc}]$'
        the_ylabel = r'$z ~ [h^{-1} \mathrm{Mpc}]$'
        cbar_label = (r'$\Phi(x=$'+str(round(Phi_slice,2))+
            r'$~h^{-1}\mathrm{Mpc},y,z)$')
        legend_string = ( r'$x\in[$'+str(round(the_zlim[0],2))+
            r'$~h^{-1}\mathrm{Mpc},$'+str(round(the_zlim[1],2))+
            r'$~h^{-1}\mathrm{Mpc})$')

        # Plots corresponding slice of Phi
        ax = axs.pcolormesh( Y[0,:,:], Z[0,:,:], Phi[Phi_cells,:,:],
                             cmap='viridis')

    else:
        the_plane  = r'$(x,z)$'
        Phi_cells  = (yspan[0]+yspan[1]) // 2
        Phi_slice  = -s/2+s/n*Phi_cells
        the_xlim   = [-s/2+s/n*xspan[0],-s/2+s/n*xspan[1]]
        the_ylim   = [-s/2+s/n*zspan[0],-s/2+s/n*zspan[1]]
        the_zlim   = [-s/2+s/n*yspan[0],-s/2+s/n*yspan[1]]
        #the_title  = (r'$\Phi(x,y=$'+str(Phi_slice)+
        #    r'$~h^{-1}\mathrm{Mpc},z),~\mathrm{peaks~at}~y\in[$'+
        #    str(round(the_zlim[0],2))+r'$~h^{-1}\mathrm{Mpc},$'+
        #    str(round(the_zlim[1],2))+r'$~h^{-1}\mathrm{Mpc})$')
        the_xlabel = r'$x ~ [h^{-1} \mathrm{Mpc}]$'
        the_ylabel = r'$z ~ [h^{-1} \mathrm{Mpc}]$'
        cbar_label = (r'$\Phi(x,y=$'+str(round(Phi_slice,2))+
            r'$~h^{-1}\mathrm{Mpc},z)$')
        legend_string = ( r'$y\in[$'+str(round(the_zlim[0],2))+
            r'$~h^{-1}\mathrm{Mpc},$'+str(round(the_zlim[1],2))+
            r'$~h^{-1}\mathrm{Mpc})$')

        # Plots corresponding slice of Phi
        ax = axs.pcolormesh( X[:,0,:], Z[:,0,:],Phi[:,Phi_cells,:],
                             cmap='viridis')

    # Setting plot formatting
    #axs.set_title(the_title, fontsize=16)
    axs.set_xlabel(the_xlabel, fontsize=14)
    axs.set_ylabel(the_ylabel, fontsize=14)
    ax_pos = make_axes_locatable(axs)
    ax_colorbar_pos = ax_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(ax, cax=ax_colorbar_pos)
    ax_colorbar_pos.yaxis.set_label_position('right')
    ax_colorbar_pos.set_ylabel(cbar_label, labelpad=15, rotation=270,
                               fontsize=14)
    axs.set_xlim(the_xlim)
    axs.set_ylim(the_ylim)
    if the_xlim[1]-the_xlim[0] == the_ylim[1]-the_ylim[0]:
        axs.set_aspect(1)

    for a in range(len( q_pk[:,0] )):
        # Plots peak patches as circles centred at xL,yL,zL with radius Rth
        mpp, = axs.plot( Rth[a]*np.cos(theta)+r_pk[a,0], 
                         Rth[a]*np.sin(theta)+r_pk[a,1], 
                         color='white')

        # Plot 1LPT peak displacements
        lpt = axs.arrow( q_pk[a,0], q_pk[a,1],
            r_pk1LPT[a,0]-q_pk[a,0], r_pk1LPT[a,1]-q_pk[a,1], color='w',
            ls=':')#head_width=1, head_length=1, length_includes_head=True )

        # Plot 2LPT peak displacements
        lpt2 = axs.arrow( r_pk1LPT[a,0], r_pk1LPT[a,1],
            r_pk[a,0]-r_pk1LPT[a,0], r_pk[a,1]-r_pk1LPT[a,1], color='w',
            )#head_width=1, head_length=1, length_includes_head=True )

        # Plot theoretical great attractors as arrows
        if lambda_1[a] >= 0:
            ga = axs.arrow( q_pk[a,0], q_pk[a,1],
                q_GA[a,0]-q_pk[a,0], q_GA[a,1]-q_pk[a,1], color='k',ls='-.',
                head_width=1, head_length=1, length_includes_head=True )

        # # Plot theoretical great attractors as points
        # GA, = axs.plot( q_GA[a,0], q_GA[a,1], ls='none', marker='x',
        #     mfc='none', mec='k' )

    # Plots minima of the Phi field
    axs.plot( PhiL[0]+s/2/n, PhiL[1]+s/2/n, ls='none', marker='x',
        mfc='white', mec='red')#, label=r'Minima of $\Phi(\mathbf{x})$' )

    # Legend
    axs.legend(loc='center', bbox_to_anchor=(0., 1., 1., .15), ncol=3,
        facecolor='grey', 
        title=(r'For a slice '+legend_string+' projected onto the '+
            the_plane+r' plane'),
        handles=(
            Line2D([0],[0],ls='none',marker='o',mfc='none',mec='w',
                markersize=10,label=r'Peak patches'),
            Line2D([0],[0],ls='none',marker='x',mec='r',
                label=r'Minima of $\Phi(\mathbf{q},$'+R_label+r'$)$'),
            Line2D([0],[0],ls=':',color='w',
                label=r'$\mathbf{s}_{pk}^{(1)}$'),
            Line2D([0],[0],ls='-',color='w',
                label=r'$\mathbf{s}_{pk}^{(2)}$'),
            Line2D([0],[0],ls='-.',color='k',
                label=r'$\mathbf{q}_{GA}$') ), fancybox=True)

    # Displays plot or returns figure object
    if display==True:
         plt.show()
    else:
         return fig,axs




def plot_great_attractor_prog_( q_pk, r_pk, r_pk1LPT, Rth, Phi, q_GA, 
        s, n, xspan=[0,neff], yspan=[0,neff], zspan=[neff//2,neff//2+1],
        display=False):
    from matplotlib.lines import Line2D

    # Convert Phi from Mpc^2 s^-2 to km^2 s^-2
    Phi = 9.5214061369e38*Phi

    # Discarding peaks that lie outside of the slice
    w = np.intersect1d(
          np.intersect1d(
            np.intersect1d( np.argwhere( q_pk[:,0]>=-s/2+s/n*xspan[0] ) ,
                            np.argwhere( q_pk[:,0]< -s/2+s/n*xspan[1] ) ),
            np.intersect1d( np.argwhere( q_pk[:,1]>=-s/2+s/n*yspan[0] ) ,
                            np.argwhere( q_pk[:,1]< -s/2+s/n*yspan[1] ) )),
          np.intersect1d( np.argwhere( q_pk[:,2]>=-s/2+s/n*zspan[0] ),
                          np.argwhere( q_pk[:,2]< -s/2+s/n*zspan[1] )))
    q_pk     = q_pk[w]
    r_pk     = r_pk[w]
    q_GA     = q_GA[w]
    Rth      = Rth[w].reshape(len(w))
    r_pk1LPT = r_pk1LPT[w]
    
    # Dependent variable used to plot circular peak patches
    theta = np.linspace(0, 2*np.pi, 100, endpoint=True)
    
    # Define bin edges for Phi plot
    edges = np.linspace( -s/2 , s/2 , n+1 )
    X,Y,Z = np.meshgrid(edges,edges,edges,indexing='ij')
 
    R_th_scales = np.linspace(0,2.5*RTHmax,20)
    for q in range(20):
        if q==0:
            Phi_W = Phi
            R_label = (r'$R_{\mathrm{lattice}} = $'+str(round(
                boxsize/neff/2,2 ))+r' $h^{-1}$Mpc')
        else:
            Phi_W = signal.convolve( Phi, W_th( qX,qY,qZ,R_th_scales[q] ),
                mode='same')
            R_label = (r'$R_{th} = $'+str(round( R_th_scales[q],2 ))+
                r' $h^{-1}$Mpc')

        fig,axs = plt.subplots(1,1, figsize=(12,12))

        PhiL = minima_and_maxima(Phi_W, neff, R=1, out='min')
        PhiM = minima_and_maxima(Phi_W, neff, R=1, out='max')

        # Discarding Phi minima that lie outside the slice
        w = np.intersect1d(
                np.intersect1d(
                    np.intersect1d( np.argwhere( PhiL[0,:] >= xspan[0] ),
                                     np.argwhere( PhiL[0,:] <  xspan[1] )),
                    np.intersect1d( np.argwhere( PhiL[1,:] >= yspan[0] ),
                                    np.argwhere( PhiL[1,:] <  yspan[1] ))),
                np.intersect1d( np.argwhere( PhiL[2,:] >= zspan[0] ),
                                np.argwhere( PhiL[2,:] <  zspan[1] )))
        PhiL = -s/2 + s/n * PhiL[:,w]

        # Discarding Phi maxima that lie outside the slice
        w = np.intersect1d(
                np.intersect1d(
                    np.intersect1d( np.argwhere( PhiM[0,:] >= xspan[0] ),
                                     np.argwhere( PhiM[0,:] <  xspan[1] )),
                    np.intersect1d( np.argwhere( PhiM[1,:] >= yspan[0] ),
                                    np.argwhere( PhiM[1,:] <  yspan[1] ))),
                np.intersect1d( np.argwhere( PhiM[2,:] >= zspan[0] ),
                                np.argwhere( PhiM[2,:] <  zspan[1] )))
        PhiM = -s/2 + s/n * PhiM[:,w]

        # Determines what plane the volume slice is in and defines figure
        # labels accordingly
        the_plane  = r'$(x,y)$'
        Phi_cells  = (zspan[0]+zspan[1]) // 2
        Phi_slice  = -s/2+s/n*Phi_cells
        the_xlim   = [-s/2+s/n*xspan[0],-s/2+s/n*xspan[1]]
        the_ylim   = [-s/2+s/n*yspan[0],-s/2+s/n*yspan[1]]
        the_zlim   = [-s/2+s/n*zspan[0],-s/2+s/n*zspan[1]]
        the_xlabel = r'$x ~ [h^{-1} \mathrm{Mpc}]$'
        the_ylabel = r'$y ~ [h^{-1} \mathrm{Mpc}]$'
        cbar_label = (r'$\Phi(x,y,z=$'+str(round(Phi_slice,2))+
             r'$~h^{-1}\mathrm{Mpc})~[\mathrm{km}^2\,\mathrm{s}^{-2}]$')

        # Plots corresponding slice of Phi
        ax = axs.pcolormesh( X[:,:,0], Y[:,:,0], Phi_W[:,:,Phi_cells],
            cmap='viridis')

        for a in range(len( q_pk[:,0] )):
            # Plots peak patches as circles centred at xL,yL,zL with radius Rth
            mpp, = axs.plot( Rth[a]*np.cos(theta)+r_pk[a,0], 
                         Rth[a]*np.sin(theta)+r_pk[a,1], 
                         color='white')

            # Plot 1LPT peak displacements
            lpt = axs.arrow( q_pk[a,0], q_pk[a,1],
                r_pk1LPT[a,0]-q_pk[a,0], r_pk1LPT[a,1]-q_pk[a,1], color='w',
                ls=':')#head_width=1, head_length=1, length_includes_head=True )

            # Plot 2LPT peak displacements
            lpt2 = axs.arrow( r_pk1LPT[a,0], r_pk1LPT[a,1],
                r_pk[a,0]-r_pk1LPT[a,0], r_pk[a,1]-r_pk1LPT[a,1], color='w',
                )#head_width=1, head_length=1, length_includes_head=True )

            # Plot theoretical great attractors as arrows
            if (a in flags) and (lambda_1[a]>=0):
                if np.linalg.norm(q_GA[a]-r_pk[a]) < boxsize/2:
                    ga = axs.arrow( r_pk[a,0], r_pk[a,1],
                        q_GA[a,0]-r_pk[a,0], q_GA[a,1]-r_pk[a,1], color='k',ls='-.',
                        head_width=1, head_length=1, length_includes_head=True )

        # Plots minima of the Phi field
        axs.plot( PhiL[0]+s/2/n, PhiL[1]+s/2/n, ls='none', marker='x',
            mfc='w', mec='r')#, label=r'Minima of $\Phi(\mathbf{x})$' )

        # Plots maxima of the Phi field
        axs.plot( PhiM[0]+s/2/n, PhiM[1]+s/2/n, ls='none', marker='+',
            mfc='w', mec='violet')#, label=r'Minima of $\Phi(\mathbf{x})$' )

        # Setting plot formatting
        axs.set_title(R_label)
        axs.set_xlabel(the_xlabel, fontsize=14)
        axs.set_ylabel(the_ylabel, fontsize=14)
        ax_pos = make_axes_locatable(axs)
        ax_colorbar_pos = ax_pos.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(ax, cax=ax_colorbar_pos)
        ax_colorbar_pos.yaxis.set_label_position('right')
        ax_colorbar_pos.set_ylabel(cbar_label, labelpad=15, rotation=270,
                fontsize=14)
        axs.set_xlim(the_xlim)
        axs.set_ylim(the_ylim)
        if the_xlim[1]-the_xlim[0] == the_ylim[1]-the_ylim[0]:
            axs.set_aspect(1)
        
        plt.show()












def plot_great_attractor_prog( q_pk, r_pk, r_pk1LPT, Rth, Phi, q_GA, 
        s, n, xspan=[0,neff], yspan=[0,neff], zspan=[neff//2,neff//2+1],
        display=False):
    ''' INPUTS/OUTPUTS:
        'PhiL' is a 4D NumPy array containing the positions of each of the
        local maxima it is of the same form as 'Phi_local_maxima' and 
        'Phi_local_minima' defined in the main body of the code above. The
        input parameters 'xspan', 'yspan', and 'zspan' are lists of two
        integers and 'display' is a boolian parameter. If display=True,
        there is no output, and the function simply displays the figure
        generated. If display=False, then 'fig' (of the class
        'matplotlib.figure.Figure') is returned. The remaining input
        parameters share the class of like-named variables in the main code
        body. 
        
        DESCRIPTION:
        This function plots the Phi peaks (found using a sequential
        nearest-neighbours search) and peak patches that lie within the
        spatial bounds laid out by 'xspan', 'yspan' and 'zspan' super-
        imposed on the Phi slice defined by the middle of this region (i.e.
        if zspan[1]-zspan[0]<xspan[1]-xspan[0],yspan[1]-yspan[0] then the
        slice at Phi[x,y,z=(zspan[1]-zspan[0])/2] is plotted).

        USAGE:
        This function can be used to directly display a plot (as described
        above) to the screen, e.g. with
        >>> plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,PhiL=Phi_local_maxima)
        
        Or to return a matplotlib.figure.Figure class object, which allows
        us to display multiple figures at once, e.g. with
        >>> figure4 = plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,
            PhiL=Phi_local_maxima,xspan=[0,n],yspan=[0,n],
            zspan=[n//2,n//2+1],display=False)
        >>> plt.show()
    '''
    from matplotlib.lines import Line2D

    # Convert Phi from Mpc^2 s^-2 to km^2 s^-2
    Phi = 9.5214061369e38*Phi

    # Discarding peaks that lie outside of the slice
    w = np.intersect1d(
          np.intersect1d(
            np.intersect1d( np.argwhere( q_pk[:,0]>=-s/2+s/n*xspan[0] ) ,
                            np.argwhere( q_pk[:,0]< -s/2+s/n*xspan[1] ) ),
            np.intersect1d( np.argwhere( q_pk[:,1]>=-s/2+s/n*yspan[0] ) ,
                            np.argwhere( q_pk[:,1]< -s/2+s/n*yspan[1] ) )),
          np.intersect1d( np.argwhere( q_pk[:,2]>=-s/2+s/n*zspan[0] ),
                          np.argwhere( q_pk[:,2]< -s/2+s/n*zspan[1] )))
    q_pk     = q_pk[w]
    r_pk     = r_pk[w]
    q_GA     = q_GA[w]
    Rth      = Rth[w].reshape(len(w))
    r_pk1LPT = r_pk1LPT[w]
    
    # Dependent variable used to plot circular peak patches
    theta = np.linspace(0, 2*np.pi, 100, endpoint=True)
    
    # Define bin edges for Phi plot
    edges = np.linspace( -s/2 , s/2 , n+1 )
    X,Y,Z = np.meshgrid(edges,edges,edges,indexing='ij')
 
    # Creates figure
    fig,axs = plt.subplots(2,2, figsize=(12,12))

    R_th_scales = np.linspace(0,2*RTHmax,4)
    for q in range(4):
        if q==0:
            Phi_W = Phi
            R_label = (r'$R_{\mathrm{lattice}} = $'+str(round( boxsize/neff/2,2 ))+
                r' $h^{-1}$Mpc')
        else:
            Phi_W = signal.convolve( Phi, W_th( qX,qY,qZ,R_th_scales[q] ),
                mode='same')
            R_label = (r'$R_{th} = $'+str(round( R_th_scales[q],2 ))+
                r' $h^{-1}$Mpc')

        PhiL = minima_and_maxima(Phi_W, neff, R=1, out='min')
        PhiM = minima_and_maxima(Phi_W, neff, R=1, out='max')

        # Discarding Phi minima that lie outside the slice
        w = np.intersect1d(
                np.intersect1d(
                    np.intersect1d( np.argwhere( PhiL[0,:] >= xspan[0] ),
                                     np.argwhere( PhiL[0,:] <  xspan[1] )),
                    np.intersect1d( np.argwhere( PhiL[1,:] >= yspan[0] ),
                                    np.argwhere( PhiL[1,:] <  yspan[1] ))),
                np.intersect1d( np.argwhere( PhiL[2,:] >= zspan[0] ),
                                np.argwhere( PhiL[2,:] <  zspan[1] )))
        PhiL = -s/2 + s/n * PhiL[:,w]

        # Discarding Phi maxima that lie outside the slice
        w = np.intersect1d(
                np.intersect1d(
                    np.intersect1d( np.argwhere( PhiM[0,:] >= xspan[0] ),
                                     np.argwhere( PhiM[0,:] <  xspan[1] )),
                    np.intersect1d( np.argwhere( PhiM[1,:] >= yspan[0] ),
                                    np.argwhere( PhiM[1,:] <  yspan[1] ))),
                np.intersect1d( np.argwhere( PhiM[2,:] >= zspan[0] ),
                                np.argwhere( PhiM[2,:] <  zspan[1] )))
        PhiM = -s/2 + s/n * PhiM[:,w]

        # Determines what plane the volume slice is in and defines figure
        # labels accordingly
        if ( zspan[1]-zspan[0]<=xspan[1]-xspan[0] and
             zspan[1]-zspan[0]<=yspan[1]-yspan[0] ):
            the_plane  = r'$(x,y)$'
            Phi_cells  = (zspan[0]+zspan[1]) // 2
            Phi_slice  = -s/2+s/n*Phi_cells
            the_xlim   = [-s/2+s/n*xspan[0],-s/2+s/n*xspan[1]]
            the_ylim   = [-s/2+s/n*yspan[0],-s/2+s/n*yspan[1]]
            the_zlim   = [-s/2+s/n*zspan[0],-s/2+s/n*zspan[1]]
            the_xlabel = r'$x ~ [h^{-1} \mathrm{Mpc}]$'
            the_ylabel = r'$y ~ [h^{-1} \mathrm{Mpc}]$'
            cbar_label = (r'$\Phi(x,y,z=$'+str(round(Phi_slice,2))+
                 r'$~h^{-1}\mathrm{Mpc})~[\mathrm{km}^2\,\mathrm{s}^{-2}]$')
            legend_string = ( r'$z\in[$'+str(round(the_zlim[0],2))+
                r'$~h^{-1}\mathrm{Mpc},$'+str(round(the_zlim[1],2))+
                r'$~h^{-1}\mathrm{Mpc})$')

            # Plots corresponding slice of Phi
            ax = axs[q//2,q%2].pcolormesh( X[:,:,0], Y[:,:,0], Phi_W[:,:,Phi_cells],
                cmap='viridis')

        elif ( xspan[1]-xspan[0]<=yspan[1]-yspan[0] and
               xspan[1]-xspan[0]<=zspan[1]-zspan[0] ):
            the_plane  = 'r$(y,z)$'
            Phi_cells  = (xspan[0]+xspan[1]) // 2
            Phi_slice  = -s/2+s/n*Phi_cells
            the_xlim   = [-s/2+s/n*yspan[0],-s/2+s/n*yspan[1]]
            the_ylim   = [-s/2+s/n*zspan[0],-s/2+s/n*zspan[1]]
            the_zlim   = [-s/2+s/n*xspan[0],-s/2+s/n*xspan[1]]
            the_xlabel = r'$y ~ [h^{-1} \mathrm{Mpc}]$'
            the_ylabel = r'$z ~ [h^{-1} \mathrm{Mpc}]$'
            cbar_label = (r'$\Phi(x=$'+str(round(Phi_slice,2))+
                r'$~h^{-1}\mathrm{Mpc},y,z)$')
            legend_string = ( r'$x\in[$'+str(round(the_zlim[0],2))+
                r'$~h^{-1}\mathrm{Mpc},$'+str(round(the_zlim[1],2))+
                r'$~h^{-1}\mathrm{Mpc})$')

            # Plots corresponding slice of Phi
            ax = axs[q//2,q%2].pcolormesh( Y[0,:,:], Z[0,:,:], Phi_W[Phi_cells,:,:],
                cmap='viridis')

        else:
            the_plane  = r'$(x,z)$'
            Phi_cells  = (yspan[0]+yspan[1]) // 2
            Phi_slice  = -s/2+s/n*Phi_cells
            the_xlim   = [-s/2+s/n*xspan[0],-s/2+s/n*xspan[1]]
            the_ylim   = [-s/2+s/n*zspan[0],-s/2+s/n*zspan[1]]
            the_zlim   = [-s/2+s/n*yspan[0],-s/2+s/n*yspan[1]]
            the_xlabel = r'$x ~ [h^{-1} \mathrm{Mpc}]$'
            the_ylabel = r'$z ~ [h^{-1} \mathrm{Mpc}]$'
            cbar_label = (r'$\Phi(x,y=$'+str(round(Phi_slice,2))+
                r'$~h^{-1}\mathrm{Mpc},z)$')
            legend_string = ( r'$y\in[$'+str(round(the_zlim[0],2))+
                r'$~h^{-1}\mathrm{Mpc},$'+str(round(the_zlim[1],2))+
                r'$~h^{-1}\mathrm{Mpc})$')

            # Plots corresponding slice of Phi
            ax = axs[q//2,q%2].pcolormesh( X[:,0,:], Z[:,0,:],Phi_W[:,Phi_cells,:],
                cmap='viridis')

        for a in range(len( q_pk[:,0] )):
            # Plots peak patches as circles centred at xL,yL,zL with radius Rth
            mpp, = axs[q//2,q%2].plot( Rth[a]*np.cos(theta)+r_pk[a,0], 
                         Rth[a]*np.sin(theta)+r_pk[a,1], 
                         color='white')

            # Plot 1LPT peak displacements
            lpt = axs[q//2,q%2].arrow( q_pk[a,0], q_pk[a,1],
                r_pk1LPT[a,0]-q_pk[a,0], r_pk1LPT[a,1]-q_pk[a,1], color='w',
                ls=':')#head_width=1, head_length=1, length_includes_head=True )

            # Plot 2LPT peak displacements
            lpt2 = axs[q//2,q%2].arrow( r_pk1LPT[a,0], r_pk1LPT[a,1],
                r_pk[a,0]-r_pk1LPT[a,0], r_pk[a,1]-r_pk1LPT[a,1], color='w',
                )#head_width=1, head_length=1, length_includes_head=True )

            # Plot theoretical great attractors as arrows
            if (a in flags) and (lambda_1[a]>=0):
                if np.linalg.norm(q_GA[a]-r_pk[a]) < boxsize/2:
                    ga = axs[q//2,q%2].arrow( r_pk[a,0], r_pk[a,1],
                        q_GA[a,0]-r_pk[a,0], q_GA[a,1]-r_pk[a,1], color='k',ls='-.',
                        head_width=1, head_length=1, length_includes_head=True )

            # Plots minima of the Phi field
            axs[q//2,q%2].plot( PhiL[0]+s/2/n, PhiL[1]+s/2/n, ls='none', marker='x',
                mfc='w', mec='r')#, label=r'Minima of $\Phi(\mathbf{x})$' )

            # Plots maxima of the Phi field
            axs[q//2,q%2].plot( PhiM[0]+s/2/n, PhiM[1]+s/2/n, ls='none', marker='+',
                mfc='w', mec='violet')#, label=r'Minima of $\Phi(\mathbf{x})$' )

        # Setting plot formatting
        if q//2 != 0:
            axs[q//2,q%2].set_xlabel(the_xlabel, fontsize=14)
        if q%2 == 0:
            axs[q//2,q%2].set_ylabel(the_ylabel, fontsize=14)
        ax_pos = make_axes_locatable(axs[q//2,q%2])
        ax_colorbar_pos = ax_pos.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(ax, cax=ax_colorbar_pos)
        ax_colorbar_pos.yaxis.set_label_position('right')
        if q%2 != 0:
            ax_colorbar_pos.set_ylabel(cbar_label, labelpad=15, rotation=270,
                fontsize=14)
        axs[q//2,q%2].set_xlim(the_xlim)
        axs[q//2,q%2].set_ylim(the_ylim)
        if the_xlim[1]-the_xlim[0] == the_ylim[1]-the_ylim[0]:
            axs[q//2,q%2].set_aspect(1)

        # Legend
        axs[q//2,q%2].legend(loc='center', bbox_to_anchor=(0.03, 1., 1., .2), ncol=3,
            facecolor='grey', 
            handles=(
                Line2D([0],[0],ls='none',marker='o',mfc='none',mec='w',
                    markersize=10,label=r'Peak patches'),
                Line2D([0],[0],ls='none',marker='x',mec='r',
                    label=r'$\Phi_{\mathrm{min}}(\mathbf{q},$'+R_label+r'$)$'),
                Line2D([0],[0],ls='-.',color='k',
                    label=r'$\mathbf{q}_{GA}$'),
                Line2D([0],[0],ls='none',marker='+',mec='violet',
                    label=r'$\Phi_{\mathrm{max}}$'),
                Line2D([0],[0],ls=':',color='w',
                    label=r'$\mathbf{s}_{pk}^{(1)}$'),
                Line2D([0],[0],ls='-',color='w',
                    label=r'$\mathbf{s}_{pk}^{(2)}$') ), fancybox=True)
    fig.subplots_adjust(top=0.88)
    fig.suptitle=(r'For a slice '+legend_string+' projected onto the '+
        the_plane+r' plane')

    # Displays plot or returns figure object
    if display==True:
        plt.show()
    else:
        return fig,axs













###########################################################################
### Plotting density slices and histograms                              ###
###########################################################################

fig0,axs0=plot_density_hist( delta, X,Y,
    (r'a)', r'$x$ [$h^{-1}$Mpc]', r'$y$ [$h^{-1}$Mpc]', r'$\delta_L$'),#r'$\rho-\bar{\rho}$'),
    (r'b)', r'$\delta_L$', r'counts per overdensity bin'),#r'$\bar{\rho}\delta = \rho-\bar{\rho}$', r'counts per overdensity bin'),
    show=False)

# Test to see what a run with f_{NL}>0 should look like
f_NL    = 10
avg     = np.mean(delta**2)
avgk    = np.mean( np.abs(deltak)**2 )

#deltaNL = delta + f_NL*(delta**2+avg)
deltakNL= deltak + f_NL*( np.abs(deltak)**2+avgk )
deltaNL = scipy.fft.ifftn( deltakNL ).real

fig01,axs01 = plot_density_hist( deltaNL, X,Y,
    (r'a)', r'$x$ [$h^{-1}$Mpc]', r'$y$ [$h^{-1}$Mpc]', r'$\delta_L$$'),#r'$y$ [$h^{-1}$Mpc]', r'$\rho-\bar{\rho}$'),
    (r'b)', r'$\delta_L$', r'counts per overdensity bin'),#r'$\bar{\rho}\delta = \rho-\bar{\rho}$', r'counts per overdensity bin'),
    show=False)

# Check if fNL field exist and plot
if os.path.isfile( run_dir + '/fields/Fvec_fNL_' + run_name ):
    in_deltafNL = open( run_dir + '/fields/Fvec_fNL_' + run_name , 'rb' )
    deltafNL = np.fromfile(in_deltafNL,dtype=np.float32,count=-1)
    deltafNL = np.reshape(deltafNL, (nlattice,nlattice,nlattice), order='F')
    deltafNL = deltafNL[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]

    fig1,axs1=plot_density_hist( deltafNL, X,Y,
        (r'a)', r'$x$ [$h^{-1}$Mpc]', r'$y$ [$h^{-1}$Mpc]', r'$\delta_{NL}$'),#r'$\rho_{fNL}-\bar{\rho}_{fNL}$'),
        (r'b)', r'$\delta_{NL}$', r'counts per overdensity bin'),#r'$\bar{\rho}_{NL}\delta_{fNL} = \rho_{fNL}-\bar{\rho}_{fNL}$', r'counts per overdensity bin'),
        show=False)

# Check if zetag field exists and plot
if os.path.isfile( run_dir + '/fields/zetag_' + run_name ):
    in_zetag = open( run_dir + '/fields/zetag_' + run_name , 'rb' )
    zetag = np.fromfile(in_zetag,dtype=np.float32,count=-1)
    zetag = np.reshape(zetag, (nlattice,nlattice,nlattice), order='F')
    zetag = zetag[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]

    fig2,axs2=plot_density_hist( zetag, X,Y,
        (r'a)', r'$x$ [$h^{-1}$Mpc]', r'$y$ [$h^{-1}$Mpc]', r'$\zeta$ or something'),
        (r'b)', r'$\zeta$ or something', r'counts per $\zeta$ bin'),
        show=False)

plt.show()















###########################################################################
### Making plots                                                        ###
###########################################################################

# Lagrangian space scale
q_axes   = np.linspace(0,boxsize,neff,endpoint=False)-boxsize*(1-1/neff)/2
qX,qY,qZ = np.meshgrid(q_axes,q_axes,q_axes,indexing='ij')

# Angle between 2LPT displacement and GA displacement
"""
s_qGA_theta = np.arccos( np.array([ s_pk[i].dot(q_GA[i]-q_pk[i]) /
    ( np.linalg.norm(s_pk[i]) * np.linalg.norm(q_GA[i]-q_pk[i]))
    for i in range(Non) ]))

fig_theta,ax_theta = plt.subplots(1,2,figsize=(12.,5.))

# Plot of angle vs smoothing radius
ax_theta[0].plot( Rth, s_qGA_theta, ls='none', marker='.' )
ax_theta[0].set_xlabel(r'$R_{th,c}$ $[h^{-1}$ Mpc$]$')
ax_theta[0].set_ylabel(r'$\theta_{GA,2LPT} = \mathrm{arccos}\left(\frac{\mathbf{s}_c\cdot(\mathbf{q}_{GA,c}-\mathbf{q}_c)}{|\mathbf{s}_c||\mathbf{q}_{GA,c}-\mathbf{q}_c|}\right)$ [rad]')

# Plot of angle vs lambda_v,i
ax_theta[1].plot( lambda_3, s_qGA_theta, ls='none', marker='.',
    label=r'$\lambda_{c,3}$')
ax_theta[1].plot( lambda_2, s_qGA_theta, ls='none', marker='.',
    label=r'$\lambda_{c,2}$')
ax_theta[1].plot( lambda_1, s_qGA_theta, ls='none', marker='.',
    label=r'$\lambda_{c,1}$')
ax_theta[1].set_xlabel(r'$\lambda_{c,i}$')
ax_theta[1].legend()

fig_taylor_eval,ax_taylor_eval = plt.subplots(1,2,figsize=(12.,5.))
qGA_minus_qpk_norm = np.array([
    np.linalg.norm( q_GA[i]-q_pk[i] ) for i in range(Non) ])
w = np.argwhere( lambda_1 >= 0 )

ax_taylor_eval[0].plot(s_qGA_theta, qGA_minus_qpk_norm/Rth, ls='none', marker='.' )
#ax_taylor_eval[0].plot(s_qGA_theta[w] , qGA_minus_qpk_norm[w]/Rth[w], ls='none', marker='.' )
ax_taylor_eval[0].set_xlabel(r'$\theta_{GA,2LPT} = \mathrm{arccos}\left(\frac{\mathbf{s}_c\cdot(\mathbf{q}_{GA,c}-\mathbf{q}_c)}{|\mathbf{s}_c||\mathbf{q}_{GA,c}-\mathbf{q}_c|}\right)$ [rad]')
ax_taylor_eval[0].set_ylabel(r'$\frac{| \mathbf{q}_{GA,c}-\mathbf{q}_{c} |}{R_{th,c}}$')
ax_taylor_eval[0].set_yscale('log')

ax_taylor_eval[1].plot(lambda_3, qGA_minus_qpk_norm/Rth, ls='none', marker='.',
    label=r'$\lambda_{c,3}$')
ax_taylor_eval[1].plot(lambda_2, qGA_minus_qpk_norm/Rth, ls='none', marker='.',
    label=r'$\lambda_{c,2}$')
ax_taylor_eval[1].plot(lambda_1, qGA_minus_qpk_norm/Rth, ls='none', marker='.',
    label=r'$\lambda_{c,1}$')
ax_taylor_eval[1].set_xlabel(r'$\lambda_{c,i}$')
ax_taylor_eval[1].set_yscale('log')
ax_taylor_eval[1].legend()
#ax_taylor_eval[0].set_ylim([0,128])
#ax_taylor_eval[1].set_ylim([0,128])
plt.show()
#"""

# Mag GA distance vs Rth and lambda
"""
qGA_minus_qc = np.array([ np.linalg.norm( q_GA[i]-q_pk[i] ) for i in range(Non) ])

fig_norm,ax_norm = plt.subplots(1,2,figsize=(12.,5.))
# Plot of angle vs smoothing radius
ax_norm[0].plot( Rth, qGA_minus_qc, ls='none', marker='.' )
ax_norm[0].set_xlabel(r'$R_{th,c}$ $[h^{-1}$ Mpc$]$')
ax_norm[0].set_ylabel(r'$|\mathbf{q}_{GA,c}-\mathbf{q}_c|$ $[h^{-1}$ Mpc$]$')

# Plot of angle vs lambda_v,i
ax_norm[1].plot( lambda_3, qGA_minus_qc, ls='none', marker='.', ms=2,
    label=r'$\lambda_{c,3}$')
ax_norm[1].plot( lambda_2, qGA_minus_qc, ls='none', marker='.', ms=2,
    label=r'$\lambda_{c,2}$')
ax_norm[1].plot( lambda_1, qGA_minus_qc, ls='none', marker='.', ms=2,
    label=r'$\lambda_{c,1}$')
ax_norm[1].set_xlabel(r'$\lambda_{c,i}$')
ax_norm[1].legend()

ax_norm[0].set_ylim([0, 128])
ax_norm[1].set_ylim([0, 128])
plt.show()
"""

# Lambda vs Rth
"""
fig_R_lambda,ax_R_lambda = plt.subplots(1,2)
ax_R_lambda[0].plot(lambda_1,Rth,ls='none',marker='.',label=r'$i=1$')
ax_R_lambda[0].plot(lambda_2,Rth,ls='none',marker='.',label=r'$i=2$')
ax_R_lambda[0].plot(lambda_3,Rth,ls='none',marker='.',label=r'$i=3$')

w = np.argwhere(lambda_1 >= 0)
ax_R_lambda[1].plot(lambda_1[w],Rth[w],ls='none',marker='.',label=r'$i=1$')
ax_R_lambda[1].plot(lambda_2[w],Rth[w],ls='none',marker='.',label=r'$i=2$')
ax_R_lambda[1].plot(lambda_3[w],Rth[w],ls='none',marker='.',label=r'$i=3$')

the_xlims = ax_R_lambda[0].get_xlim()
the_ylims = ax_R_lambda[0].get_ylim()
ax_R_lambda[1].set_xlim(the_xlims)
ax_R_lambda[1].set_ylim(the_ylims)

ax_R_lambda[0].set_ylabel( r'$R_{th,c}$ $[h^{-1}$ Mpc$]$' )
ax_R_lambda[0].set_xlabel( r'$\lambda_{c,i}$' )
ax_R_lambda[1].set_xlabel( r'$\lambda_{c,i}$' )
ax_R_lambda[0].legend()
ax_R_lambda[1].legend()
plt.show()
#"""

# Plot of principle axis displacement vs eigenvalues
"""
fig_principle,ax_principle = plt.subplots(3,1,figsize=(12.,10.))

ax_principle[0].plot(
    lambda_pk[:,0], s_pk_pr[:,0]/lambda_pk[:,0], ls='none', 
    marker='.', mec='b', mfc='b', label=r'$i=1$' )
ax_principle[1].plot(
    lambda_pk[:,1], s_pk_pr[:,1]/lambda_pk[:,1], ls='none', 
    marker='.', mec='r', mfc='r', label=r'$i=2$' )
ax_principle[2].plot(
    lambda_pk[:,2], s_pk_pr[:,2]/lambda_pk[:,2], ls='none',
    marker='.', mec='purple', mfc='purple', label=r'$i=3$' )
#plt.title( r'Principle axis peak displacement ${s^\prime}_c^i$ over eigenvalue $\lambda_{c,i}$' )
ax_principle[2].set_xlabel( r'$\lambda_{c,i}$' )
ax_principle[0].set_ylabel(
    r'${s^\prime}_c^1 (\lambda_{c,1})^{-1}$ $[h^{-1}$ Mpc$]$' )
ax_principle[1].set_ylabel(
    r'${s^\prime}_c^2 (\lambda_{c,2})^{-1}$ $[h^{-1}$ Mpc$]$' )
ax_principle[2].set_ylabel(
    r'${s^\prime}_c^3 (\lambda_{c,3})^{-1}$ $[h^{-1}$ Mpc$]$' )
ax_principle[0].legend()
ax_principle[1].legend()
ax_principle[2].legend()
plt.show()
#"""

# Phi scanning through q_z
"""#Phi_W  = signal.convolve( Phi, W_th( qX,qY,qZ,10 ), mode='same')
Phi_W  = Phi
Phi_local_minima = minima_and_maxima(Phi_W, neff, out='min')

fig=[[],[],[],[]]
for j in range(4):

    fig[j] = plot_great_attractors(q_pk, r_pk, r_pk1LPT, Rth,
        Phi_W, Phi_local_minima, q_GA, boxsize, neff,
        xspan=[0,neff], yspan=[0,neff],
        zspan=[j*12,(j+1)*12])[0]
#"""

# Phi at various filter scales
"""
figs=[[],[],[],[],[]]
R_th_scales = np.linspace(0,2*RTHmax,5)
# globalmin   = ( (np.argmin(Phi)//(len(Phi)**2))%len(Phi) ,
#                 (np.argmin(Phi)//len(Phi))%len(Phi)      ,
#                 np.argmin(Phi)%len(Phi)                  )
for j in range(5):
    if j==0:
        Phi_W = Phi
        Phi_local_minima = minima_and_maxima(Phi_W, neff, R=1, out='min')
        print('Found ', len(Phi_local_minima[0]) ,
            'minima at R_th = lattice scale')
        the_label = (r'$R_{th} = $'+str(round( boxsize/neff/2,2 ))+
            r' $h^{-1}$Mpc')
    else:
        Phi_W  = signal.convolve( Phi, W_th( qX,qY,qZ,R_th_scales[j] ),
            mode='same')
        Phi_local_minima = minima_and_maxima(Phi_W, neff, 1,
            #int( 2*neff*(R_th_scales[j])**(1/3)/boxsize )
            out='min')
        print('Found ', len(Phi_local_minima[0]) , 'minima at R_th = ' ,
            R_th_scales[j])
        the_label = (r'$R_{th} = $'+str(round( R_th_scales[j],2 ))+
            r' $h^{-1}$Mpc')

    figs[j] = plot_great_attractors(q_pk, r_pk, r_pk1LPT, Rth,
        Phi_W, Phi_local_minima, q_GA, boxsize, neff,
        xspan=[neff//4,neff], yspan=[neff//8,neff-neff//8], zspan=[8,15],
        R_label=the_label)[0]
plt.show()
#"""

# Phi at various filter scales in one subplot
"""
zrange=[6,17]
print('z in [', round(zrange[0]/neff*boxsize-boxsize/2,2), ', ', round(zrange[1]/neff*boxsize-boxsize/2,2), ')')
fig_GAs,axs_GAs = plot_great_attractor_prog(q_pk, r_pk, r_pk1LPT,
        Rth*scale_factor*100, Phi, q_GA, boxsize, neff,
        #xspan=[neff//4,neff], yspan=[neff//8,neff-neff//8], 
        xspan=[0,neff], yspan=[0,neff], 
        zspan=zrange)

plt.show()
#"""
"""
print(np.average(Fpk))


zrange=[6,17]
print('z in [', round(zrange[0]/neff*boxsize-boxsize/2,2), ', ', round(zrange[1]/neff*boxsize-boxsize/2,2), ')')
plot_great_attractor_prog_(q_pk, r_pk, r_pk1LPT,
        Rth*scale_factor*100, Phi, q_GA, boxsize, neff,
        xspan=[0,neff], yspan=[0,neff], 
        zspan=zrange)
"""


# Number of Phi minima as a function of R_th and R_g
"""
R_th_scales = np.linspace(boxsize/neff, 2*RTHmax, 150)
R_Gg_scales = np.linspace(boxsize/neff/2, 2*RTHmax, 150)
N_minima    = np.zeros( np.shape(R_th_scales) )
N_mins_g    = np.zeros( np.shape(R_th_scales) )
for j in range( len(R_th_scales) ):
    Phi_W = signal.convolve( Phi, W_th( qX,qY,qZ,R_th_scales[j] ), mode='same')
    Phi_local_minima = minima_and_maxima(Phi_W, neff, 1, out='min')
    Phi_G = signal.convolve( Phi, W_G( qX,qY,qZ,R_Gg_scales[j]/2 ), mode='same')
    Phi_minima_gauss = minima_and_maxima(Phi_G, neff, 1, out='min')
    N_minima[j] = len(Phi_local_minima[0])
    N_mins_g[j] = len(Phi_minima_gauss[0])

figure,axis = plt.subplots(1)
axis.plot( R_th_scales, N_minima, label=r'$R_{th} = R_W$' )
axis.plot( R_Gg_scales, N_mins_g, label=r'$R_{gauss} = 2R_W$' )
axis.set_ylabel(r'Local $\Phi$ minima found $N_{\mathrm{minima}}$')
axis.set_xlabel(r'Filter scale $R_{W}$ $[h^{-1}$Mpc$]$')
plt.legend()
plt.show()
#"""




"""






def plot_peak_displacements(xL,yL,zL,x,y,z,display=True):
    ''' INPUTS/OUTPUTS:
        All inputs except 'display' are share the class of like-named
        variables in the main code body. If display=True, there is no
        output, and the function simply displays the figure generated. If
        display=False, then 'vectorplotobj' (a Mayavi plot object of class
        'mayavi.modules.vectors.Vectors') is returned.
        
        DESCRIPTION:
        Plots the vector displacement of the peak patches identified and
        passed to the code at command line with the pointer
        <merged_peak_file>.
        The free plotting library Mayavi is used to plot vectors indicating
        the displacement of the mass peaks from the initial (Lagrangian)
        position xL,yL,zL and the final (Eulerian) position x,y,z. Vectors
        are shown with tails at xL,yL,zL and heads at x,y,z.
        
        USAGE:
        This function can be used to directly display a plot (as described
        above) to the screen, e.g. with
        >>> plot_peak_displacements(xL,yL,zL,x,y,z)
        
        Or to return a 'mayavi.modules.vectors.Vectors' class object, which
        allows us to display multiple figures at once, e.g. with
        >>> maya1 = plot_peak_displacements(xL,yL,zL,x,y,z,display=False)
        >>> mlab.show()
    '''
    # Plots/returns plot of displacement vectors of each  peak
    vectorplotobj = mlab.quiver3d(xL,yL,zL,x-xL,y-yL,z-zL)
    if display==True:
        mlab.show()
    else:
        return vectorplotobj


def plot_peaks(xL,yL,zL,Rth,display=True):
    ''' INPUTS/OUTPUTS:
        All inputs except 'display' share the class of like-named variables
        in the main code body. If display=True, there is no output, and the
        function simply displays the figure generated. If display=False,
        then 'pointplotobj' (a Mayavi plot object of class
        'mayavi.modules.glyph.Glyph') is returned.
        
        DESCRIPTION:
        Plots the mass peaks identified by a peak patch run and passed to
        the code at command line with the pointer <merged_peak_file>. The
        free plotting library Mayavi is used to plot peaks as spheres with
        the radius (in Lagrangian space) approximated by the largest top
        hat filter radius 'Rth' at which the peak was identified (in
        reality, the radius is greater than or equal to 'Rth'). The
        (Lagrangian space) positions of the N peaks are, in terms of
        components, the N elements of xL,yL,zL.
        
        USAGE:
        This function can be used to directly display a plot (as described
        above) to the screen, e.g. with
        >>> plot_peaks(xL,yL,zL,Rth)
        
        Or to return a 'mayavi.modules.glyph.Glyph' class object, which
        allows us to display multiple figures at once, e.g. with
        >>> maya2 = plot_peaks(xL,yL,zL,Rth,display=False)
        >>> mlab.show()
    '''
    # Plots/returns plot of peaks
    pointplotobj = mlab.points3d(x,y,z,Rth,colormap='viridis')
    if display==True:
        mlab.show()
    else:
        return pointplotobj





Phi_W = signal.convolve( Phi, W_th(qX,qY,qZ,RTHmax), mode='same')
PhiL = minima_and_maxima(Phi_W, neff, R=1, out='min')
#PhiM = minima_and_maxima(Phi_W, neff, R=1, out='max')

q_GA_4plot = q_GA[flags]
x_4plot,y_4plot,z_4plot = x[flags] , y[flags] , z[flags]
xL4plot,yL4plot,zL4plot = xL[flags], yL[flags], zL[flags]
q_4plot = np.array([ [xL4plot[i],yL4plot[i],zL4plot[i]] for i in range(len(q_GA_4plot)) ])
r_4plot = np.array([ [x_4plot[i],y_4plot[i],z_4plot[i]] for i in range(len(q_GA_4plot)) ])

q_GA_norm = np.array([ np.linalg.norm( q_GA_4plot[a]-r_4plot[a] ) for a in range(len(q_GA_4plot)) ])

w = np.argwhere( q_GA_norm < 64. )
q_GA_4plot = q_GA_4plot[w]
q_4plot    = q_4plot[w]
r_4plot    = r_4plot[w]

mlab.quiver3d(xL,yL,zL,x-xL,x-yL,x-zL, colormap='Blues')
mlab.quiver3d(q_4plot[:,0,0],q_4plot[:,0,1],q_4plot[:,0,2],
              q_GA_4plot[:,0,0]-q_4plot[:,0,0],
              q_GA_4plot[:,0,1]-q_4plot[:,0,1],
              q_GA_4plot[:,0,2]-q_4plot[:,0,2], colormap='Reds')

mlab.points3d( ( PhiL[0]-neff/2 )*boxsize/neff ,
               ( PhiL[1]-neff/2 )*boxsize/neff ,
               ( PhiL[2]-neff/2 )*boxsize/neff ,
               scale_factor=1., color=(0,0,0) )

GAs_in_3D = mlab.gcf()
mlab.figure( figure=GAs_in_3D , bgcolor=(.5,.5,.5) )
v = mlab.view()

for i in range(-15,10,5):
    mlab.view( azimuth    = v[0]+i ,
               elevation  = 90     ,
               distance   = v[2]   ,
               focalpoint = v[3]   )

    mlab.savefig('afig'+str(i)+'.png',size=(800,800),figure=GAs_in_3D)

mlab.show()
#"""
"""
#mlab.quiver3d(xL,yL,zL,x-xL,y-yL,z-zL)

#mlab.points3d(x,y,z,Rth,colormap='viridis')
mlab.points3d(xL,yL,zL,Rth,colormap='inferno')
mlab.quiver3d(xL,yL,zL,x-xL,x-yL,x-zL, colormap='inferno')
peaks = mlab.gcf()
mlab.figure(figure=peaks, bgcolor=(0,0,0))
v = mlab.view()

for i in range(0,360,5):
    v2 = (v[0]+10. , v[1], v[2], v[3])

    mlab.view( azimuth    = v[0]+i ,
               elevation  = v[1]   ,
               distance   = v[2]   ,
               focalpoint = v2[3]  )

    mlab.savefig('fig'+str(i)+'.png',size=(800,800),figure=peaks)

#mlab.show()



# Plot the peak pathces in the simmulation volume
#maya1 = plot_peaks(x,y,z,Rth,display=False)

# Display Mayavi figures
#mlab.show()

#xLscale = q_axes


#plot_peak_displacements(xL,yL,zL,x,y,z,display=True)
"""

def plot_ellipsoids(q_pk,R_1,R_2,R_3,Rth,rotmatx):
    '''
    Input parameters:
    q_pk: NumPy array of length N, Lagrangian position of ellipsoid centres.
    R_1, R_2, R_3: NumPy arrays of length N with semiaxes of each ellipsoid
        in Lagrangian space where R_3 <= R_2 <= R_1, R_1 is the last axis
        to collapse.
    Rth: NumPy array of length N of top-hat filter scale at which peak was
        found.
    rotmatx: NumPy array with dimension (N,3,3), the rotation matrix
        describing the orientation of each of the N ellipsoids.

    where N is a common integer, corresponding to the number of ellipsoids
    to be plotted.

    Output:
    None, prints a Mayavi plot to screen.

    Given relevant parameters, this function plots N ellipsoids in space.
    '''
    import matplotlib.cm
    cmap = matplotlib.cm.get_cmap('viridis')
    Rth0 = (Rth - np.min(Rth)) / (np.max(Rth) - np.min(Rth))

    dtheta, dphi = np.pi / 25. , np.pi / 25.
    [theta, phi] = np.mgrid[ 0 : np.pi + dtheta*1.5 : dtheta ,
                             0 : 2*np.pi + dphi*1.5 : dphi   ]
    
    colorscale = (Rth - np.min(Rth)) / (np.max(Rth) - np.min(Rth))

    for i in range(len(R_1)):
        x = R_1[i] * np.sin(theta) * np.cos(phi)
        y = R_2[i] * np.sin(theta) * np.sin(phi)
        z = R_3[i] * np.cos(theta)
        
        X = np.array([ [x[j],y[j],z[j]] for j in range(len(phi)) ])
        Y = np.array([ rotmatx[i].dot(X[j]) for j in range(len(phi)) ])

        mlab.mesh( Y[:,0] + float(q_pk[i,0]) ,
                   Y[:,1] + float(q_pk[i,1]) ,
                   Y[:,2] + float(q_pk[i,2]) , color=cmap(Rth0[i])[:3] )
        
        del(x,y,z,X,Y)
    
    #mlab.axes()
    mlab.view( azimuth    = 45.0,
               elevation  = 54.7,
               distance   = 465.3,
               focalpoint = [-2.13211303, -0.48169468,  0.83177272] )

    mlab.show()


plot_ellipsoids(q_pk, R_1, R_2, R_3, Rth, P)

# Check overlap
"""
void_count, halo_count, over_count = 0,0,0
for i in range(1000):
    rand_coords = boxsize * np.random.rand(3)
    flag = False
    flaG = False
    for j in range(len(Rth)):
        if np.linalg.norm(rand_coords-q_pk[j,:].T) < 2*Rth[j]:
            if flag==True:
                flaG = True
            flag = True
    if flaG == True:
        over_count += 1
    elif flag == True:
        halo_count += 1
    else:
        void_count += 1

print('void:  ',void_count)
print('halo:  ',halo_count)
print('over:  ',over_count)
print('total: ',void_count+halo_count+over_count)
"""

'''
def test_mesh():
    dtheta, dphi = np.pi / 250.0, np.pi / 250.0
    [theta, phi] = np.mgrid[ 0 : np.pi + dtheta*1.5 : dtheta ,
                             0 : 2*np.pi + dphi*1.5 : dphi   ]
    a = 3.
    b = 2.
    c = 1.
    x = a * np.sin(theta) * np.cos(phi) 
    y = b * np.sin(theta) * np.sin(phi)
    z = c * np.cos(theta)
    
    # Euler angles
    alpha, beta, gamma = np.pi/4, np.pi/4, np.pi/4

    Rotation_matrix = (
        np.array([[ np.cos(gamma) , -np.sin(gamma) , 0 ], 
                  [ np.sin(gamma) ,  np.cos(gamma) , 0 ], 
                  [ 0             ,  0             , 1 ]]).dot(

        np.array([[ 1 , 0            ,  0            ],   
                  [ 0 , np.cos(beta) , -np.sin(beta) ],
                  [ 0 , np.sin(beta) ,  np.cos(beta) ]])).dot(

        np.array([[ np.cos(alpha) , -np.sin(alpha) , 0 ],
                  [ np.sin(alpha) ,  np.cos(alpha) , 0 ],
                  [ 0             ,  0             , 1 ]])
        ) )

    X = np.array([ [x[i],y[i],z[i]] for i in range(len(phi)) ])
    Y = np.array([ Rotation_matrix.dot(X[i]) for i in range(len(phi)) ])

    mlab.mesh(Y[:,0]+3, Y[:,1]+3, Y[:,2]+3, colormap="viridis")
    mlab.mesh(x, y, z, colormap="viridis")
    mlab.axes()
    mlab.show()

def test_mesh2():
    dtheta, dphi = np.pi / 250.0, np.pi / 250.0
    [theta, phi] = np.mgrid[ 0 : np.pi + dtheta*1.5 : dtheta ,
                             0 : 2*np.pi + dphi*1.5 : dphi   ]
    a = 3. 
    b = 2. 
    c = 1. 
    r = ( (np.sin(theta)*np.cos(phi)/a)**2 +
          (np.sin(theta)*np.sin(phi)/b)**2 +
          (np.cos(theta)/c)**2 )**-.5
    x = r * np.sin(theta+np.pi) * np.cos(phi) 
    y = r * np.sin(theta+np.pi) * np.sin(phi)
    z = r * np.cos(theta+np.pi)
    
    mlab.mesh(x, y, z, colormap="viridis")
    mlab.show()

test_mesh()
'''



'''
# Euler angles
alpha, beta, gamma = np.pi/2, 0., 0.

Rotation_matrix = (
    np.array([[ np.cos(gamma) , -np.sin(gamma) , 0 ],
              [ np.sin(gamma) ,  np.cos(gamma) , 0 ],
              [ 0             ,  0             , 1 ]]).dot(

    np.array([[ 1 , 0            ,  0            ],
              [ 0 , np.cos(beta) , -np.sin(beta) ],
              [ 0 , np.sin(beta) ,  np.cos(beta) ]])).dot(

    np.array([[ np.cos(alpha) , -np.sin(alpha) , 0 ],
              [ np.sin(alpha) ,  np.cos(alpha) , 0 ],
              [ 0             ,  0             , 1 ]])
    ) )






# Plot Ellipsoidal patches
#"""

# Ellipoid solid angle
dtheta, dphi = np.pi / 250.0, np.pi / 250.0
[theta, phi] = np.mgrid[ 0 : np.pi + dtheta*1.5 : dtheta ,
                         0 : 2*np.pi + dphi*1.5 : dphi   ]

# Ellipsoide semiaxes
a     = .5 * np.exp( - lambda_1 + Fpk/3 ) * Rth
b     = .5 * np.exp( - lambda_2 + Fpk/3 ) * Rth
c     = .5 * np.exp( - lambda_3 + Fpk/3 ) * Rth

# Orienting ellipsoids
ell_x      = a * np.sin(theta) * np.cos(phi)
ell_y      = b * np.sin(theta) * np.sin(phi)
ell_z      = c * np.cos(theta)
ellipsoids = np.array([ [ell_x[i],ell_y[i],ell_z[i]]
    for i in range(len(phi)) ])
ellipsoids = np.array([ P.dot(


for i in range(Non):
    
    

'''







