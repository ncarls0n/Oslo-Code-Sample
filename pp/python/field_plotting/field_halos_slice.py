import numpy as np
import matplotlib as mpl
import subprocess,os

s = str(subprocess.check_output('hostname'))
if s[:5]+'####'+s[9:] == "b'nia####.scinet.local\n'":
    mpl.use('Agg')

import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy
import sys
import os
try:
    plt.style.use('nate')
except OSError:
    plt.style.use('default')

#fieldfile = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.04.29_LIM_production_runs/ng0_cenz7150Mpc_C/fields/Fvec_500Mpc_n432_nb74_nt3'
fieldfile = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.04.29_LIM_production_runs/ng7_cenz7150_C_mlambda20_sigmag/fields/Fvec_500Mpc_n432_nb74_nt3'
halofile  = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.04.29_LIM_production_runs/ng7_cenz7150_C_mlambda20_sigmag/output/500Mpc_n432_nb74_nt3_merge.pksc.13579_M_gt_1.129e+11Msol'

# Read catalogue file header
in_cat  = open( halofile, 'rb' )
N_halos = np.fromfile( in_cat, dtype=np.int32  , count=1 )[0]
RTHmax  = np.fromfile( in_cat, dtype=np.float32, count=1 )[0]
zin     = np.fromfile( in_cat, dtype=np.float32, count=1 )[0]

# Determine number of columns in halo catalogue table
size_bytes = os.path.getsize(halofile)
cols       = int( (size_bytes-12)/4/N_halos )

# Read halo catalogue
catalogue = np.reshape( np.fromfile( in_cat, dtype=np.float32, count=N_halos*cols ), (N_halos,cols) )

# Keep only halos whos Langrangian position is within 1 filter radius of yL=0
lowpass   = np.where( catalogue[:,4] <  catalogue[:,6] )[0] ; catalogue = catalogue[lowpass ,:]
highpass  = np.where( catalogue[:,4] > -catalogue[:,6] )[0] ; catalogue = catalogue[highpass,:]


nmesh   = 432
nbuff   = 74
boxsize = 500.0
ntile   = 3
cenz    = 7150.



# Additional parameters
nlattice = int( (nmesh-2*nbuff) * ntile + 2*nbuff )
neff     = int( nlattice - 2*nbuff )

# x, y, & z axes for field plots
edges = np.linspace( -boxsize/2 , boxsize/2 , neff+1 )
X,Y = np.meshgrid(edges,edges,indexing='ij')
Y += cenz
cut  = [int(3*neff/8),int(5*neff/8)]
xcut = [ X[ cut[0],    0  ] , X[ cut[1],    0  ] ]
zcut = [ Y[     0 ,cut[0] ] , Y[     0 ,cut[1] ] ]
lowpass   = np.where( catalogue[:,7] < xcut[1] )[0] ; catalogue = catalogue[lowpass ,:]
highpass  = np.where( catalogue[:,7] > xcut[0] )[0] ; catalogue = catalogue[highpass,:]
lowpass   = np.where( catalogue[:,9] < zcut[1] )[0] ; catalogue = catalogue[lowpass ,:]
highpass  = np.where( catalogue[:,9] > zcut[0] )[0] ; catalogue = catalogue[highpass,:]
x,y,z=catalogue[:,:3].T ; R,xL,yL,zL = catalogue[:,6:10].T
del(catalogue)

def read_field_slice(fieldfile,plane,intercept,nlattice,nbuff):
    # fieldfile is a string containing the field binary file
    # plane = 'xy', 'xz', or 'yz' and defines the plane
    # intercept is the intercept with plane and axis perpendicular to it
    #     (-nlattice+nbuff < intercept < nlattice-nbuff)
    # nlattice is full resolution of the box including buffers
    # nbuff is the buffer thickness
    #
    # Returns 2D NumPy array with specified slice of the 3D field

    field = np.array([])
    neff  = nlattice-2*nbuff   # Lattice dimension excluding buffers
    y_offset = 4 * nlattice    # byte offset to go from (i,j,k)->(i,j+1,k)
    z_offset = 4 * nlattice**2 # byte offset to go from (i,j,k)->(i,j,k+1)

    # For the x-y plane f(x,y,z=intercept)->field[i,j]
    if 'z' not in plane.lower():
        offset   = 4 * (nlattice**2*(nbuff+intercept)+nlattice*nbuff+nbuff)
        for j in range(neff):
            f = np.fromfile( fieldfile, dtype=np.float32,
                             count  = nlattice-2*nbuff,
                             offset = offset + j*y_offset )
            field = np.concatenate((field,f))

    # For the x-z plane f(x,y=intercept,z)->field[i,k]
    elif 'y' not in plane.lower():
        offset   = 4 * (nlattice**2*nbuff+nlattice*(nbuff+intercept)+nbuff)
        for k in range(neff):
            f = np.fromfile( fieldfile, dtype=np.float32,
                             count=nlattice-2*nbuff,
                             offset=offset + k*z_offset   )
            field = np.concatenate( (field,f) )

    # For the y-z plane f(x=intercept,y,z)->field[j,k]
    elif 'x' not in plane.lower():
        offset   = 4 * (nlattice**2*nbuff+nlattice*nbuff+nbuff+intercept)
        for k in range(neff):
            for j in range(neff):
                f = np.fromfile( fieldfile, dtype=np.float32, count=1,
                                 offset = offset + j*y_offset + k*z_offset )
                field = np.concatenate( (field,f) )

    return np.reshape( field, (neff,neff), order='F' )

# Function to setup a field slice subplot
def slice_subplot(fig,ax,fieldslice,y0,y1,title,xtitle,ytitle):
    fig.colorbar( ax.pcolormesh(y0,y1,fieldslice,cmap='viridis', vmin=-20, vmax=20),
        cax=make_axes_locatable(ax).
        append_axes('right', size='5%', pad=0.05)   )
    ax.set_aspect(1)                ; ax.set_xlabel(xtitle,fontsize=16)
    ax.set_title(title,fontsize=20) ; ax.set_ylabel(ytitle,fontsize=16)
    return fig,ax


fig,ax = plt.subplots( nrows=1, ncols=1 )
xy_slice = read_field_slice(fieldfile,'xz',0,nlattice,nbuff)
fig,ax = slice_subplot( fig, ax, xy_slice[cut[0]:cut[1],cut[0]:cut[1]], X[cut[0]:cut[1],cut[0]:cut[1]], Y[cut[0]:cut[1],cut[0]:cut[1]],
    r'$\delta(x,y=0,z)$' ,r'$x$ [Mpc]',r'$z$ [Mpc]')

theta = np.linspace( 0, 2*np.pi, 100 )
for i in range(len(x)):
    ax.plot( xL[i]+R[i]*np.cos(theta) , zL[i]+R[i]*np.sin(theta) , ls='-', color='w' )

for i in range(len(x)):
    ax.arrow( xL[i],zL[i],x[i]-xL[i],z[i]-zL[i],color='w')

ax.set_xlim(xcut)
ax.set_ylim(zcut)

fig.set_size_inches(14,12)
#fig.savefig(os.path.dirname(fieldfile)+'/delta.png',bbox_inches='tight')
fig.savefig('/cita/d/www/home/njcarlson/cifar_meeting/varying_lambdav2H-2/33.92997_halos_delta.png',bbox_inches='tight')
#        r'$\delta_{L,nG}(x,y,z=0)$'
#        r'$\delta_{L,nG}(x,y=0,z)$'
#        r'$\delta_{L,nG}(x=0,y,z)$'
