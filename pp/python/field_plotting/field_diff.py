import numpy as np
import matplotlib as mpl
import subprocess

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

q=[ ['1'    , '8.48249' ],
    ['2.236','19.08561' ],
    ['2.5'  ,'33.92997' ],
    ['2.75' ,'53.01558' ],
    ['3'    ,'76.34243' ],
    ['4'    ,''         ] ]
s='/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.04.29_LIM_production_runs/ng7_cenz7150Mpc_C_mlambda25_a{0}e-51/fields/{1}_500Mpc_n432_nb74_nt3'
j=5
fieldfile1 = s.format(q[j][0],'Fvec')
fieldfile2 = s.format(q[j][0],'rhog')
outfile = '/cita/d/www/home/njcarlson/cifar_meeting/varying_ae/{0}e-51_sigma8planck_deltang.png'.format(q[j][0])
nmesh = 432
nbuff = 74
boxsize = 500.0
ntile = 3

# Additional parameters
nlattice = int( (nmesh-2*nbuff) * ntile + 2*nbuff )
neff     = int( nlattice - 2*nbuff )

# x, y, & z axes for field plots
edges = np.linspace( -boxsize/2 , boxsize/2 , neff+1 )
X,Y = np.meshgrid(edges,edges,indexing='ij')


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
    fig.colorbar( ax.pcolormesh(y0,y1,fieldslice,cmap='viridis',vmin=-10,vmax=40),
        cax=make_axes_locatable(ax).
        append_axes('right', size='5%', pad=0.05)   )
    ax.set_aspect(1)                ; ax.set_xlabel(xtitle,fontsize=16)
    ax.set_title(title,fontsize=20) ; ax.set_ylabel(ytitle,fontsize=16)
    return fig,ax


fig,ax = plt.subplots( nrows=1, ncols=1 )
xy_slice = read_field_slice(fieldfile1,'xy',0,nlattice,nbuff) - read_field_slice(fieldfile2,'xy',0,nlattice,nbuff)
fig,ax = slice_subplot( fig, ax, xy_slice[:,:], X[:,:], Y[:,:],
    r'$\delta(x,y,z=0)-\delta_G(x,y,z=0)$' ,r'$x$ [Mpc]',r'$y$ [Mpc]')
fig.set_size_inches(14,12)
fig.savefig(outfile,bbox_inches='tight')

#        r'$\delta_{L,nG}(x,y,z=0)$'
#        r'$\delta_{L,nG}(x,y=0,z)$'
#        r'$\delta_{L,nG}(x=0,y,z)$'
