import numpy as np
import matplotlib as mpl
import subprocess

s = str(subprocess.check_output('hostname'))
if s[:5]+'####'+s[9:] == "b'nia####.scinet.local\n'":
    mpl.use('Agg')

import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy
import sys
import os
#try:
#    plt.style.use('nate')
#except OSError:
#    plt.style.use('default')

"""
USAGE:

python3 python/plotting/run_checks/check_fields.py <run> <field>

where <field> = "delta", "deltag", "deltang", "zeta", "zetag", "zetang",
"chi" or "all" (assumes "all" if no value for <field> given).
"""

fieldfile1 = sys.argv[1]
fieldfile2 = sys.argv[2]
#fieldfile3 = sys.argv[3]
#fieldfile4 = sys.argv[4]

if len(sys.argv)>5:
    title1 = sys.argv[5]
    title2 = sys.argv[6]
    title3 = sys.argv[7]
    title4 = sys.argv[8]
else:
    title1 = r'$\zeta_G(\mathbf{q})$'
    title2 = r'$\zeta(\mathbf{q})$, $f_{nG}=10^3$'
    title3 = r'$\zeta(\mathbf{q})$, $f_{nG}=3.75\times10^3$'
    title4 = r'$\zeta(\mathbf{q})$, $f_{nG}=5\times10^3$'

boxsize  = 4000.
nmesh    = 236
nbuff    = 20
ntile    = 10
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


# def read_maxmin(fieldfile,nlattice,nbuff):
#     # fieldfile is a string containing the field binary file
#     # nlattice is full resolution of the box including buffers
#     # nbuff is the buffer thickness
#     # 
#     # Returns the max and min of the field (excluding buffers)
#     # 
#     # Calls read_field_slice()
# 
# 
#     fieldmax , fieldmin = -np.inf , np.inf
#     neff = nlattice - 2*nbuff
# 
#     for i in range(neff):
#         field_slice = read_field_slice(fieldfile,'yz',i,nlattice,nbuff)
#         slicemax , slicemin = np.max(field_slice) , np.min(field_slice)
#         fieldmax            = np.max([fieldmax,slicemax])
#         fieldmin            = np.min([fieldmin,slicemin])
#         del(field_slice)
# 
#     return fieldmax,fieldmin
# 
# 
# def read_hist(fieldfile,nlattice,nbuff):
#     # fieldfile is a string containing the field binary file
#     # nlattice is full resolution of the box including buffers
#     # nbuff is the buffer thickness
#     #    
#     # Returns NumPy arrays representing bin edges and bin weights for a
#     # histogram of the field file
#     # 
#     # Calls read_maxmin() and read_field_slice()
# 
#     neff = nlattice - 2*nbuff
#     fieldmax,fieldmin = read_maxmin( fieldfile, nlattice, nbuff )
# 
#     if fieldmin == fieldmax:
#         bin_edges = np.linspace( fieldmin-1e-6, fieldmax+1e-6, 101 )
# 
#     else:
#         # Estimate interquartile range to calculate number of bins
#         field = read_field_slice(fieldfile,'x',0,nlattice,nbuff)
#         IQR   = scipy.stats.iqr(field) ; del(field)
# 
#         # number of bins using Freedman-Draconis rule
#         num_bins  = int( (fieldmax-fieldmin) * neff / 2 / IQR +.5 )
#         bin_edges = np.linspace(fieldmin,fieldmax,num_bins)
#         counts    = np.zeros(( len(bin_edges)-1 ))
# 
#     for i in range(neff):
#         fieldslice = read_field_slice(fieldfile,'x',i,nlattice,nbuff)
#         counts_i,junk = np.histogram( fieldslice, bins=bin_edges )
#         counts += counts_i
#         del(fieldslice,counts_i)
# 
#     return bin_edges,counts


# Function to setup a field slice subplot
def slice_subplot(fig,ax,fieldslice,y0,y1,title,xtitle,ytitle):
    fig.colorbar( ax.pcolormesh(y0,y1,fieldslice,cmap='viridis'),
                  cax=make_axes_locatable(ax).
                      append_axes('right', size='5%', pad=0.05)   )
    ax.set_aspect(1)    ; ax.set_xlabel(xtitle)
    ax.set_title(title) ; ax.set_ylabel(ytitle)
    return fig,ax


fig,axs = plt.subplots(1,4)

field1 = read_field_slice(fieldfile1,'xy',int(neff/2),nlattice,nbuff)
fig,axs[0] = slice_subplot( fig, axs[0], field1, X,Y,title1,r'$x$ [Mpc]',r'$y$ [Mpc]')

field2 = read_field_slice(fieldfile2,'xy',int(neff/2),nlattice,nbuff)
fig,axs[1] = slice_subplot( fig, axs[1], field1+field2, X,Y,title2,r'$x$ [Mpc]',r' ')
del(field2)
axs[1].axes.yaxis.set_ticklabels([])

#field3 = read_field_slice(fieldfile3,'xy',int(neff/2),nlattice,nbuff)
#fig,axs[2] = slice_subplot( fig, axs[2], field1+field3, X,Y,title3,r'$x$ [Mpc]',r' ')
#del(field3)
#axs[2].axes.yaxis.set_ticklabels([])
#
#field4 = read_field_slice(fieldfile4,'xy',int(neff/2),nlattice,nbuff)
#fig,axs[3] = slice_subplot( fig, axs[3], field1+field4, X,Y,title4,r'$x$ [Mpc]',r' ')
#del(field4)
#axs[3].axes.yaxis.set_ticklabels([])

fig.set_size_inches(12,2.75)
fig.savefig('delta.png',dpi=500)

