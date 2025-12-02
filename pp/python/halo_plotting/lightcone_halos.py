# plot-great-attract.py
#
# Nathan J. Carlson
# September 9, 2021
#
# Optimized to run on Python 3.8.2
#

import numpy as np
import matplotlib.pyplot as plt
import os

nmesh = 236
nbuff = 20
ntile = 10
boxsize = 4000.
Omx = 0.2645
OmB = 0.0493
h = 0.6735

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

cat0 = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/merges/4000Mpc_n236_nb20_nt10_merge.pksc.13579_fnl{0}'
cat = [ cat0.format('1e4') ,
        cat0.format('5e3') ,
        cat0.format('1e3') ,
        cat0.format('5e2') ,
        cat0.format('0')   ]

label= [ 'fNL=1e4'  ,
         'fNL=5e3'  ,
         'fNL=1e3'  ,
         'fNL=5e2'  ,
         'Gaussian' ]

colours = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']

fig = plt.figure()
ax  = fig.add_subplot(projection='polar')

for i in range(len(cat)):
# Read size of halo catalogue in bytes
    size_bytes = os.path.getsize(cat[i])

    # Reading data from merged peak patch catalogue pointer
    Non       = np.fromfile(cat[i],dtype=np.int32,count=1)[0]
    RTHmax    = np.fromfile(cat[i],dtype=np.float32,count=1)[0]
    zin       = np.fromfile(cat[i],dtype=np.float32,count=1)[0]
    outnum = int( (size_bytes-12)/4/Non )
    catalogue = np.fromfile(cat[i],dtype=np.float32,count=outnum*Non)
    catalogue = np.reshape(catalogue,(Non,outnum))

    # Defines X,Y,Z as meshgrid
    s,n = boxsize,neff
    edges = np.linspace( -s/2 , s/2 , n+1 )

    # Next, for the sake of clarity, we copy the raw halo data from the XY
    # columns of `peakdata` to XY better-labelled NumPy arrays:
    x      = catalogue[:,0]  # x,y,z: components of the final (Eulerian) halo
    y      = catalogue[:,1]  #     position vector [h^-1 Mpc]
    z      = catalogue[:,2]  # vx,vy,vz: components of the final (Eulerian) halo
    #vx     = catalogue[:,3]  #     velocity vector [km/s/Mpc]
    #vy     = catalogue[:,4]  # Rth: top-hat filter radius halo was found at
    #vz     = catalogue[:,5]  #     (roughly the radius of the halo) [h^-1 Mpc]
    Rth    = catalogue[:,6]  # xL,yL,zL: components of the initial (Lagrangian)
    #xL     = catalogue[:,7]  #     halo position vector [h^-1 Mpc]
    #yL     = catalogue[:,8]  # Fpk: linearly extrapolated initial overdensity
    #zL     = catalogue[:,9]  #     field at peak centre, also called Fcollv in

    # Let's free up some hard drive space
    del(catalogue)
    # M = 4*np.pi/3*Rth**3*rho # halo masses [Msol]

    #for j in range(Non):
    for j in range(5):
        if z[j]<-10 or z[j]>10: continue
        r = np.sqrt(x[j]**2+y[j]**2)
        if r > s/2: continue
        if x[j]==0: theta = np.pi/2
        else      : theta = np.arctan(np.float32(y[j]/x[j]))
        if theta<j*len(cat) or theta>=j*len(cat)+90/len(cat): continue
        ax.scatter( theta, r, c=colours[i], s=np.pi*Rth[j]**2, alpha=0.125 )

        print(x[j], y[j], theta, r)

        if j%10000==0:
            print('Halo {0} of {1} done for field {2} of {3}'.format(j+1,Non,i+1,len(cat)))

ax.set_thetamin(0)
ax.set_thetamax(90)
fig.savefig('lightcone_halos.png',dpi=100)
