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

"""
USAGE:

python3 python/plotting/run_checks/check_fields.py <run> <field>

where <field> = "delta", "deltag", "deltang", "zeta", "zetag", "zetang",
"chi" or "all" (assumes "all" if no value for <field> given).
"""

if len(sys.argv)<2 or len(sys.argv)>3:
    # Error message displayed if you pass the wrong number of arguments.
    raise ValueError('You\'ve passed me the wrong number of arguments at t'
        +'he command line. Command\nline should take the form\n    cd <pea'
        +'k_patch_directory>\n    python3.8 python/plotting/run_checks/che'
        +'ck_fields2.0.py <run> <field>\nwhere <field> can be "delta", "de'
        +'ltag" or "deltang" for the overdensity\nfield, its Gaussian or n'
        +'on-Gaussian component, or similarly "zeta", "zetag"\nor "zetang"'
        +' for the superhorizon hypersurface invariant of curvature\npertu'
        +'rbations, "chi" for model-specific transverse inflationary field'
        +', or\n"all" for all present fields. If no value for <field> is g'
        +'iven, "all" is\nassumed.')

if len(sys.argv)==2:
    fieldlist = 'all'
else:
    fieldlist = str(sys.argv[2]).lower()

# Reading from command line prompts
run_dir = str(sys.argv[1])+'/'   # The directory of the run

# Make plots subirectory if it doesn't already exist
if not os.path.isdir('{0}plots'.format(run_dir)):
    os.system('cd {0};mkdir plots'.format(run_dir))

###########################################################################
### Execute useful lines from parameter file                            ###
###########################################################################

# Open parameter file
with open(run_dir+'param/param.params') as f:
    params = [i.strip() for i in f.readlines()]

# Define variables
# seed               initial random field seed number (an integer)
# boxsize [h^-1 Mpc] sidelength of simulation volume
# nmesh              number of cells in one tile (including buffers)
# nbuff              number of cells in buffer
# ntile              number
# sigma8             stdev of Gaussian random field smoothed at 8 h^-1 Mpc
# h                  "little h", H_0/100 km/s/Mpc
# ns                 n_s, spectral index
# fNL                Amplitude of non-Gaussianity
# run_name           string used in output files
# short_name         shorter string used in output files
for line in params:
    if ( line[:4]=='seed'  or line[:7]=='boxsize' or line[:5]=='nmesh'  or
         line[:5]=='nbuff' or line[:5]=='ntile'   or line[:6]=='sigma8' or
         line[:2]=='h '    or line[:2]=='ns'      or line[:3]=='fNL'    or  
         line[:8]=='run_name'          or line[:10]=='short_name'       ):
       exec(line)

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


def read_maxmin(fieldfile,nlattice,nbuff):
    # fieldfile is a string containing the field binary file
    # nlattice is full resolution of the box including buffers
    # nbuff is the buffer thickness
    # 
    # Returns the max and min of the field (excluding buffers)
    # 
    # Calls read_field_slice()


    fieldmax , fieldmin = -np.inf , np.inf
    neff = nlattice - 2*nbuff

    for i in range(neff):
        field_slice = read_field_slice(fieldfile,'yz',i,nlattice,nbuff)
        slicemax , slicemin = np.max(field_slice) , np.min(field_slice)
        fieldmax            = np.max([fieldmax,slicemax])
        fieldmin            = np.min([fieldmin,slicemin])
        del(field_slice)

    return fieldmax,fieldmin


def read_hist(fieldfile,nlattice,nbuff):
    # fieldfile is a string containing the field binary file
    # nlattice is full resolution of the box including buffers
    # nbuff is the buffer thickness
    #    
    # Returns NumPy arrays representing bin edges and bin weights for a
    # histogram of the field file
    # 
    # Calls read_maxmin() and read_field_slice()

    neff = nlattice - 2*nbuff
    fieldmax,fieldmin = read_maxmin( fieldfile, nlattice, nbuff )

    if fieldmin == fieldmax:
        bin_edges = np.linspace( fieldmin-1e-6, fieldmax+1e-6, 101 )

    else:
        # Estimate interquartile range to calculate number of bins
        field = read_field_slice(fieldfile,'x',0,nlattice,nbuff)
        IQR   = scipy.stats.iqr(field) ; del(field)

        # number of bins using Freedman-Draconis rule
        num_bins  = int( (fieldmax-fieldmin) * neff / 2 / IQR +.5 )
        bin_edges = np.linspace(fieldmin,fieldmax,num_bins)
        counts    = np.zeros(( len(bin_edges)-1 ))

    for i in range(neff):
        fieldslice = read_field_slice(fieldfile,'x',i,nlattice,nbuff)
        counts_i,junk = np.histogram( fieldslice, bins=bin_edges )
        counts += counts_i
        del(fieldslice,counts_i)

    return bin_edges,counts


# Function to setup a field slice subplot
def slice_subplot(fig,ax,fieldslice,y0,y1,title,xtitle,ytitle):
    fig.colorbar( ax.pcolormesh(y0,y1,fieldslice,cmap='viridis'),
                  cax=make_axes_locatable(ax).
                      append_axes('right', size='5%', pad=0.05)   )
    ax.set_aspect(1)    ; ax.set_xlabel(xtitle)
    ax.set_title(title) ; ax.set_ylabel(ytitle)
    return fig,ax


# Plotting function that plots slices of the 3D field and histogram
def plot_field_slices(fieldfile,x0=X,x1=Y,nlattice=nlattice,nbuff=nbuff,
    xy_slice_title= r'$\delta_L(x,y,z=0)$',
    xz_slice_title= r'$\delta_L(x,y=0,z)$',
    yz_slice_title= r'$\delta_L(x=0,y,z)$', hist_x_label = r'$\delta_L$',
    x_intercept   = int(neff/2), y_intercept   = int(neff/2),
    z_intercept   = int(neff/2), return_slices=False, include_hist=True  ):   
 
    fig,axs = plt.subplots( nrows=2, ncols=2 )

    # Plot histogram
    if include_hist==True:
        bin_edges,counts = read_hist(fieldfile,nlattice,nbuff)
        axs[0,1].hist(  bin_edges[:-1], bin_edges, weights=counts )
        axs[0,1].set_xlabel(hist_x_label)
        axs[0,1].set_ylabel(r'counts')
        del(bin_edges,counts) # Conserve RAM

    # Plot colormesh slices in each plane of the 3D field
    xy_slice = read_field_slice(fieldfile,'xy',z_intercept,nlattice,nbuff)
    fig,axs[0,0] = slice_subplot( fig, axs[0,0], xy_slice, x0, x1,
                                  xy_slice_title,r'$x$ [Mpc]',r'$y$ [Mpc]')
    if return_slices == False: del(xy_slice) # Conserve RAM

    xz_slice = read_field_slice(fieldfile,'xz',y_intercept,nlattice,nbuff)
    fig,axs[1,0] = slice_subplot( fig, axs[1,0], xz_slice, x0, x1,
                                  xz_slice_title,r'$x$ [Mpc]',r'$z$ [Mpc]')
    if return_slices == False: del(xz_slice) # Conserve RAM

    yz_slice = read_field_slice(fieldfile,'yz',x_intercept,nlattice,nbuff)
    fig,axs[1,1] = slice_subplot( fig, axs[1,1], yz_slice, x0, x1,
                                  yz_slice_title,r'$y$ [Mpc]',r'$z$ [Mpc]')
    
    if return_slices == False:
        return fig,axs
    else:
        return fig,axs,xy_slice,xz_slice,yz_slice


# If there's a delta field
a,b,c,d,e  = run_dir, int(boxsize), nmesh, nbuff, ntile
deltafile  = '{0}fields/Fvec_{1}Mpc_n{2}_nb{3}_nt{4}'.format(a,b,c,d,e)
deltaWfile = '{0}fields/Fvec_smooth_{1}Mpc_n{2}_nb{3}_nt{4}'.format(a,b,c,d,e)
Ldeltafile = '{0}fields/Fvec_fNL_{1}Mpc_n{2}_nb{3}_nt{4}'.format(a,b,c,d,e)
deltagfile = '{0}fields/rhog_{1}Mpc_n{2}_nb{3}_nt{4}'.format(a,b,c,d,e)
zetagfile  = '{0}fields/zetag_{1}Mpc_n{2}_nb{3}_nt{4}'.format(a,b,c,d,e)
zetafile   = '{0}fields/zetang_{1}Mpc_n{2}_nb{3}_nt{4}'.format(a,b,c,d,e)
chifile    = '{0}fields/chi_{1}Mpc_n{2}_nb{3}_nt{4}'.format(a,b,c,d,e)

#############################
###   PLOT DELTA FIELDS   ###
#############################

# Overdensity and its non-Gaussian part
if fieldlist=='deltang' or ( fieldlist=='all' and os.path.isfile(deltafile)
                             and os.path.isfile(deltagfile)              ):

    print('Plotting delta fields:\n{0}\n{1}\nand their difference'
        .format( deltafile, deltagfile ) )

    fig1,axs1,delta_xy,delta_xz,delta_yz    = plot_field_slices(deltafile,
                                                        return_slices=True)
    fig2,axs2,deltag_xy,deltag_xz,deltag_yz = plot_field_slices(deltagfile,
        xy_slice_title= r'$\delta_{L,G}(x,y,z=0)$',
        xz_slice_title= r'$\delta_{L,G}(x,y=0,z)$',
        yz_slice_title= r'$\delta_{L,G}(x=0,y,z)$',
        hist_x_label  = r'$\delta_{L,G}$', return_slices=True )

    delta_xy = delta_xy - deltag_xy; del(deltag_xy)
    delta_xz = delta_xz - deltag_xz; del(deltag_xz)
    delta_yz = delta_yz - deltag_yz; del(deltag_yz)

    fig3,axs3 = plt.subplots( nrows=2, ncols=2 )
    fig3,axs3[0,0] = slice_subplot( fig3,axs3[0,0],delta_xy,X,Y,
        r'$\delta_{L,nG}(x,y,z=0)$',r'$x$ [Mpc]',r'$y$ [Mpc]')
    del(delta_xy)
    fig3,axs3[1,0] = slice_subplot( fig3,axs3[1,0],delta_xz,X,Y,
        r'$\delta_{L,nG}(x,y=0,z)$',r'$x$ [Mpc]',r'$z$ [Mpc]')
    del(delta_xz)
    fig3,axs3[1,1] = slice_subplot( fig3,axs3[1,1],delta_yz,X,Y,
        r'$\delta_{L,nG}(x=0,y,z)$',r'$y$ [Mpc]',r'$z$ [Mpc]')
    del(delta_yz)

    fig1.set_size_inches(14,12)
    fig1.savefig('{0}plots/delta.png'.format(run_dir),dpi=100)
    fig2.set_size_inches(14,12)
    fig2.savefig('{0}plots/deltag.png'.format(run_dir),dpi=100)
    fig3.set_size_inches(14,12)
    fig3.savefig('{0}plots/deltang.png'.format(run_dir),dpi=100)

# Plot overdensity with both Gaussian & non-Gaussian components
elif fieldlist=='delta' or ( fieldlist=='all' 
                             and os.path.isfile(deltafile) ):
    print('Plotting delta field:\n{0}'.format(deltafile) )
    fig1,axs1 = plot_field_slices( deltafile )
    fig1.set_size_inches(14,12)
    fig1.savefig('{0}plots/delta.png'.format(run_dir),dpi=100)

# Plot Gaussian component of overdensity field
elif fieldlist=='deltag' or ( fieldlist=='all' 
                              and os.path.isfile(deltagfile) ):
    print('Plotting Gaussian delta field:\n{0}'.format(deltagfile) )
    fig1,axs1 = plot_field_slices( deltagfile,
        xy_slice_title= r'$\delta_{L,G}(x,y,z=0)$',
        xz_slice_title= r'$\delta_{L,G}(x,y=0,z)$',
        yz_slice_title= r'$\delta_{L,G}(x=0,y,z)$',
        hist_x_label  = r'$\delta_{L,G}$',include_hist=False )
    fig1.set_size_inches(14,12)
    fig1.savefig('{0}plots/deltag.png'.format(run_dir),dpi=100)


#############################
###   PLOT ZETA FIELDS    ###
#############################

# zeta and its non-Gaussian part
if fieldlist=='zetang' or ( fieldlist=='all' and os.path.isfile(zetafile)
                             and os.path.isfile(zetagfile)               ):

    print('Plotting zeta fields:\n{0}\n{1}\nand their difference'
        .format( zetafile, zetagfile ) )

    fig4,axs4,zeta_xy,zeta_xz,zeta_yz    = plot_field_slices(zetafile,
                                                        return_slices=True)
    fig5,axs5,zetag_xy,zetag_xz,zetag_yz = plot_field_slices(zetagfile,
        xy_slice_title= r'$\zeta_{L,G}(x,y,z=0)$',
        xz_slice_title= r'$\zeta_{L,G}(x,y=0,z)$',
        yz_slice_title= r'$\zeta_{L,G}(x=0,y,z)$',
        hist_x_label  = r'$\zeta_{L,G}$', return_slices=True )

    zeta_xy = zeta_xy - zetag_xy; del(zetag_xy)
    zeta_xz = zeta_xz - zetag_xz; del(zetag_xz)
    zeta_yz = zeta_yz - zetag_yz; del(zetag_yz)

    fig6,axs6 = plt.subplots( nrows=2, ncols=2 )
    fig6,axs6[0,0] = slice_subplot( fig6,axs6[0,0],zeta_xy,X,Y,
        r'$\zeta_{L,nG}(x,y,z=0)$',r'$x$ [Mpc]',r'$y$ [Mpc]')
    del(zeta_xy)
    fig6,axs6[1,0] = slice_subplot( fig6,axs6[1,0],zeta_xz,X,Y,
        r'$\zeta_{L,nG}(x,y=0,z)$',r'$x$ [Mpc]',r'$z$ [Mpc]')
    del(zeta_xz)
    fig6,axs6[1,1] = slice_subplot( fig6,axs6[1,1],zeta_yz,X,Y,
        r'$\zeta_{L,nG}(x=0,y,z)$',r'$y$ [Mpc]',r'$z$ [Mpc]')
    del(zeta_yz)

    fig4.set_size_inches(14,12)
    fig4.savefig('{0}plots/zeta.png'.format(run_dir),dpi=100)
    fig5.set_size_inches(14,12)
    fig5.savefig('{0}plots/zetag.png'.format(run_dir),dpi=100)
    fig6.set_size_inches(14,12)
    fig6.savefig('{0}plots/zetang.png'.format(run_dir),dpi=100)

# Plot zeta with both Gaussian & non-Gaussian components
elif fieldlist=='zeta' or ( fieldlist=='all' 
                             and os.path.isfile(zetafile) ):
    print('Plotting zeta field:\n{0}'.format(zetafile) )
    fig4,axs4 = plot_field_slices( zetafile,
        xy_slice_title= r'$\zeta_L(x,y,z=0)$',
        xz_slice_title= r'$\zeta_L(x,y=0,z)$',
        yz_slice_title= r'$\zeta_L(x=0,y,z)$',
        hist_x_label  = r'$\zeta_L$'          )
    fig4.set_size_inches(14,12)
    fig4.savefig('{0}plots/zeta.png'.format(run_dir),dpi=100)

# Plot Gaussian component of zeta field
elif fieldlist=='zetag' or ( fieldlist=='all' 
                              and os.path.isfile(zetagfile) ):
    print('Plotting Gaussian zeta field:\n{0}'.format(zetagfile) )
    fig4,axs4 = plot_field_slices( zetagfile,
        xy_slice_title= r'$\zeta_{L,G}(x,y,z=0)$',
        xz_slice_title= r'$\zeta_{L,G}(x,y=0,z)$',
        yz_slice_title= r'$\zeta_{L,G}(x=0,y,z)$',
        hist_x_label  = r'$\zeta_{L,G}$' )
    fig4.set_size_inches(14,12)
    fig4.savefig('{0}plots/zetag.png'.format(run_dir),dpi=100)


#############################
###   PLOT CHI FIELDS     ###
#############################

# Plot zeta with both Gaussian & non-Gaussian components
if fieldlist=='chi' or ( fieldlist=='all' 
                             and os.path.isfile(chifile) ):
    print('Plotting chi field:\n{0}'.format(chifile) )
    fig7,axs7 = plot_field_slices( chifile,
        xy_slice_title= r'$\chi_E(x,y,z=0)$',
        xz_slice_title= r'$\chi_E(x,y=0,z)$',
        yz_slice_title= r'$\chi_E(x=0,y,z)$',
        hist_x_label  = r'$\chi_E$',include_hist=False )
    fig7.set_size_inches(14,12)
    fig7.savefig('{0}plots/chi.png'.format(run_dir),dpi=100)


######################################
###   PLOT Laplace(delta) FIELDS   ###
######################################

# Plot Laplace(delta)
if fieldlist=='laplacedelta' or ( fieldlist=='all' 
                             and os.path.isfile(Ldeltafile) ):
    print('Plotting Laplacian(delta) field:\n{0}'.format(deltafile) )
    fig8,axs8 = plot_field_slices( Ldeltafile,
        xy_slice_title= r'$\nabla^2\delta_L(x,y,z=0)$',
        xz_slice_title= r'$\nabla^2\delta_L(x,y=0,z)$',
        yz_slice_title= r'$\nabla^2\delta_L(x=0,y,z)$',
        hist_x_label  = r'$\nabla^2\delta_L$'          )
    fig8.set_size_inches(14,12)
    fig8.savefig('{0}plots/laplacedelta.png'.format(run_dir),dpi=100)
    

###########################################################################


# plt.show()

