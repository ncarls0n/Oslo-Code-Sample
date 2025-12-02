# plot-great-attract.py
#
# Nathan J. Carlson
# September 9, 2021
#
# Optimized to run on Python 3.8.2
#

# TO RUN ME:
# 
#     python lightcone_halos2.0.p
# 

import numpy as np
from matplotlib.patches import Circle as C
from matplotlib.collections import PatchCollection as PC
import matplotlib.pyplot as plt
import os
import astropy.units as u
from astropy.cosmology import Planck18, z_at_value
import astropy.cosmology.units as cu

#nongauss = '7_m40_test' # Peak Patch non-Gaussianity model
nongauss = '0' # Peak Patch non-Gaussianity model
nmesh = 580
nbuff = 40
ntile = 2
boxsize = 1000.
Omx = 0.2645
OmB = 0.0493
h = 0.6735

nlattice = int( (nmesh-2*nbuff) * ntile + 2*nbuff )
neff     = int( nlattice - 2*nbuff )
Omm      = Omx+OmB        # Omega_m, fraction of energy that can cluster
rhocrit  = 2.775e11*h**2  # critical energy density 3H^2/8piG [Msol Mpc^-3]
rho      = rhocrit*Omm    # average CDM density [Msol Mpc^-3]
deltavir = 200            # <something defining virial collapse I think>
G        = 4.517e-48      # gravitational constant [Mpc^3 Msol^-1 s^-2]
# Note here that we use mass units of solar masses Msol, spatial units of
# megaparsecs Mpc, and time units of seconds s.

#colours = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
run_dir = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.02.05_SBsuite/ng{0}'.format(nongauss)
runs = ['cenz0500Mpc' , 'cenz1500Mpc' , 'cenz2500Mpc' , 'cenz3500Mpc' , 'cenz4500Mpc' , 'cenz5500Mpc' , 'cenz6500Mpc']
seed = [ 13579        ,  23456        ,  31317        ,  42297        ,  59841        ,  66600        ,  76543       ]
cenz = np.linspace( 500., 6500., 7 )

cat  = [ '{0}/{1}/output/{2}Mpc_n{3}_nb{4}_nt{5}_merge.pksc.{6}'.format(
    run_dir, runs[j], int(boxsize), nmesh, nbuff, ntile, seed[j] )
    for j in range(len(runs)) ]

# os.system('ls {0}'.format(cat[3]))

#fig = plt.figure( figsize=(14,2) )
#ax  = fig.add_subplot(projection='polar')
fig,ax = plt.subplots( figsize=(14,2) )
theta = np.linspace(0,2*np.pi,100)
#set(ax, 'Units', 'points' )
# for i in range(len(cat)):
# Read size of halo catalogue in bytes
for i in range(len(cat)):
    size_bytes = os.path.getsize(cat[i])

    in_cat = open( cat[i] , 'rb' )

    # Reading data from merged peak patch catalogue pointer
    Non       = np.fromfile(in_cat,dtype=np.int32,count=1)[0]
    RTHmax    = np.fromfile(in_cat,dtype=np.float32,count=1)[0]
    zin       = np.fromfile(in_cat,dtype=np.float32,count=1)[0]
    outnum = int( (size_bytes-12)/4/Non ) # typically 11 for Peak Patch circa 2023
    catalogue = np.fromfile(in_cat,dtype=np.float32,count=outnum*Non)
    catalogue = np.reshape(catalogue,(Non,outnum))

    # Defines X,Y,Z as meshgrid
    s,n = boxsize,neff
    edges = np.linspace( -s/2 , s/2 , n+1 )

    x    = catalogue[:,0]
    y    = catalogue[:,1]
    z    = catalogue[:,2]
    Rvir = catalogue[:,6] / (200.**(1./3.)) * 2 # twice virial radius

    # Make y slice
    yslice = np.argwhere( np.abs(y)<25. )
    x      = x[yslice][:,0]
    z      = z[yslice][:,0]
    Rvir   = Rvir[yslice][:,0]

    # Trim overlaps between adjacent boxes
    zoverlap = np.argwhere( np.abs(z-cenz[i])<=cenz[0] )
    z = z[zoverlap][:,0]
    x = x[zoverlap][:,0]
    Rvir = Rvir[zoverlap][:,0]

    #def t(a): return print(a.min(),a.mean(),a.max())
    #t(z)
    #t(x)

    # Let's free up some hard drive space
    del(catalogue,y,yslice,zoverlap)
    #M = 4*np.pi/3*Rth**3*rho # halo masses [Msol]

    #ax.plot( z,x,ls='none',marker='.')

    ax.plot( z,x,ls='none',marker='.',markersize=1,mfc='k',mec='k',alpha=0.05)


    #for j in range(1000):
        #ax.add_patch(plt.Circle( (z[j],x[j]), Rvir[j], edgecolor='b', facecolor='k', lw=1 ))
        #ax.plot( z[j] + Rvir[j]*np.cos(theta) , x[j] + Rvir[j]*np.sin(theta) )

    #patches = []
    #for j in range(10):
    #    c = C((z[j],x[j]),Rvir[j])
    #    patches.append(c)
    #circle = [ plt.Circle( (z[j],x[j]), Rvir[j], facecolor='k', lw=0 ) for j in range(10) ]
    #c = mpl.collections.PatchCollection(circle)

    #p = PC(patches, alpha=0.4)
    #ax.add_collection(p)

    #for j in range(10):
    #    ax.add_collection( mpl.collections.PatchCollection( plt.Circle((z[j],x[j]), Rvir[j], facecolor='k', lw=0 ) ))

    #ax.add_collection(c)
    #r         = (x**2+z**2)**.5
    #theta     = np.angle( z+x*1.0j )
    #theta_max = np.arctan( (cenz[-1]+cenz[0])/cenz[0] )
    #ax.scatter( theta, r ) #, c=colours[i])#, s=np.pi*Rth[j]**2, alpha=0.125 )

# ax.set_thetamin(-theta_max)
# ax.set_thetamax(+theta_max)

#ax.axis('equal')
ax.set_xlim( 0, cenz[-1]+cenz[0] )
ax.set_ylim( -cenz[0] , cenz[0]  )

ax.set_xlabel(r'comoving distance from observer $\chi$ (Mpc)')
ax.set_ylabel(r'$\chi$ (Mpc)')


ax2 = ax.twiny()
ax2.set_xticks( ax.get_xticks() )
ax2.set_xbound( ax.get_xbound() )



ax2.set_xticklabels( [0]+[ round(z_at_value(Planck18.comoving_distance, chi*u.Mpc ).value,2) for chi in ax.get_xticks()[1:] ]  )
ax2.set_xlabel(r'redshift $z$')

for i in range(len(cenz)):
    ax.plot( [ cenz[i]+cenz[0] , cenz[i]+cenz[0] ] , [-cenz[0],cenz[0]] , color='k' , lw=1 )



fig.savefig('{0}/lightcone_halos.pdf'.format(run_dir),dpi=200,bbox_inches='tight')


print('{0}/lightcone_halos.pdf'.format(run_dir))
#'''










# 108 #keep only a certain number of halos
# 109 if rankflag<0:
# 110     nhalos = -rankflag
# 111     dm = np.argsort(Rvir)[-nhalos:]
# 112     x=x[dm]
# 113     y=y[dm]
# 114     Rvir=Rvir[dm]
# 115     M=M[dm]

'''
133 for i in range(len(Rvir)-1):
134     ax.add_patch(py.Circle((x[i],y[i]),Rvir[i], edgecolor='b',facecolor='none',lw=1))
135 
136 #titletext=simtype+"\n"+L+" x "+L+" x "+W+" Mpc/h\n"+"z = "+Z
137 titletext=L+" x "+L+" x "+W+" Mpc\n"+"z = "+Z
138 font = {'horizontalalignment' :'left',
139         'fontsize'            : 26,
140         'multialignment'      :'left'}
141 
142 ax.set_xlabel("x [Mpc]",fontsize=26)
143 ax.set_ylabel("y [Mpc]",fontsize=26)
144 ax.tick_params(axis='both',which='major',labelsize=24)
145 t=py.title(titletext,font)
146 t.set_x(0.0)
147 py.xlim(bx*fb,bx*(1-fb))
148 py.ylim(bx*fb,bx*(1-fb))
149 
150 py.savefig(outfile+'.pdf')

'''

