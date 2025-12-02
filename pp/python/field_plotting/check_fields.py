import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy
from scipy import signal
import sys
import os
try:
    plt.style.use('nate')
except OSError:
    plt.style.use('default')

"""
How to make me go:

python3 python/run_checks/check_fields.py runs/<run_name> [field]

"""
if len(sys.argv) != 2:
    # Error message displayed if you pass the wrong number of arguments.
    print('You\'ve passed me the wrong number of arguments at the command '
        +'line.\nCommand line should take the form\ncd <peak_patch_directo'
        +'ry>\npython3.8 python/run_checks/check_fields.py runs/ns128')
    sys.exit(2)

# Reading from command line prompts
run_dir = str(sys.argv[1])+'/'   # The directory of the run

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
# Omx                Omega_CDM, energy density fraction of CDM
# OmB                Omega_B, energy density fraction of baryonic matter
# Omvac              Omega_Lambda, energy denstiy fraction of DE
# h                  "little h", H_0/100 km/s/Mpc
# ns                 n_s, spectral index
# fNL                Amplitude of non-Gaussianity
# run_name           string used in output files
# short_name         shorter string used in output files
# maximum_redshift   z_max, redshift of primordial fields
# global_redshift    z_0, redshift of Eulerian peaks
for line in params:
    if ( line[:4]=='seed'  or line[:7]=='boxsize' or line[:5]=='nmesh'  or
         line[:5]=='nbuff' or line[:5]=='ntile'   or line[:6]=='sigma8' or
         line[:3]=='Omx'   or line[:3]=='OmB'     or line[:5]=='Omvac'  or
         line[:2]=='h '    or line[:2]=='ns'      or line[:3]=='fNL'    or  
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

# x, y, & z axes for field plots
edges = np.linspace( -boxsize/2 , boxsize/2 , neff+1 )
X,Y,Z = np.meshgrid(edges,edges,edges,indexing='ij')

# Plotting function that plots slices of the 3D field as well as a histogram
def plot_field_slices( field, xmesh, ymesh, zmesh,
    z_slice=int(neff/2), y_slice=int(neff/2), x_slice=int(neff/2),
    z_slice_title=r'$\delta_L(x,y,z=0)$',
    y_slice_title=r'$\delta_L(x,y=0,z)$',
    x_slice_title=r'$\delta_L(x=0,y,z)$',
    hist_x_label =r'$\delta_L$'          ):
    
    fig,axs = plt.subplots(nrows=2, ncols=2)
    
    # Plot colormesh slice in the x-y plane of the delta field
    ax00 = axs[0,0].pcolormesh( xmesh[:,:,z_slice], ymesh[:,:,z_slice], field[:,:,z_slice], cmap='viridis')
    ax00_pos = make_axes_locatable(axs[0,0])
    ax00_colorbar_pos = ax00_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(ax00, cax=ax00_colorbar_pos)
    axs[0,0].set_aspect(1)
    axs[0,0].set_title(z_slice_title)
    axs[0,0].set_xlabel(r'$x$ [$h^{-1}$ Mpc]')
    axs[0,0].set_ylabel(r'$y$ [$h^{-1}$ Mpc]')

    # Plot colormesh slice in the x-z plane of the delta field
    ax10 = axs[1,0].pcolormesh( zmesh[:,y_slice,:], xmesh[:,y_slice,:], field[:,y_slice,:], cmap='viridis')
    ax10_pos = make_axes_locatable(axs[1,0])
    ax10_colorbar_pos = ax10_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(ax10, cax=ax10_colorbar_pos)
    axs[1,0].set_aspect(1)
    axs[1,0].set_title(y_slice_title)
    axs[1,0].set_xlabel(r'$x$ [$h^{-1}$ Mpc]')
    axs[1,0].set_ylabel(r'$y$ [$h^{-1}$ Mpc]')

    # Plot colormesh slice in the y-z plane of the delta field
    ax11 = axs[1,1].pcolormesh( ymesh[x_slice,:,:], zmesh[x_slice,:,:], field[x_slice,:,:], cmap='viridis')
    ax11_pos = make_axes_locatable(axs[1,1])
    ax11_colorbar_pos = ax11_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(ax11, cax=ax11_colorbar_pos)
    axs[1,1].set_aspect(1)
    axs[1,1].set_title(x_slice_title)
    axs[1,1].set_xlabel(r'$y$ [$h^{-1}$ Mpc]')
    axs[1,1].set_ylabel(r'$z$ [$h^{-1}$ Mpc]')

    # Plot a histogram of all delta
    max_field = np.max(field)
    min_field = np.min(field)
    
    if max_field==min_field:
        bin_edges = np.linspace(min_field-1e-6, max_field+1e-6, 101)
    else:
        num_bins  = int( (max_field-min_field)             # number of bins using
            * len(field) / (2*scipy.stats.iqr(field)) +.5) # Freedman-Draconis rule
        bin_edges = np.linspace(min_field, max_field, num_bins)
    counts,junk  = np.histogram(field, bins=bin_edges) # Histogram weights
    axs[0,1].hist(  bin_edges[:-1], bin_edges, weights=counts )
    axs[0,1].set_xlabel(hist_x_label)
    axs[0,1].set_ylabel(r'counts')
    
    return fig,axs


# If there's a delta field
deltafile = '{0}fields/Fvec_{1}Mpc_n{2}_nb{3}_nt{4}'.format(
    run_dir, int(boxsize), nmesh, nbuff, ntile )
smoothfile = '{0}fields/Fvec_smooth_{1}Mpc_n{2}_nb{3}_nt{4}'.format(
    run_dir, int(boxsize), nmesh, nbuff, ntile )
deltagfile = '{0}fields/rhog_{1}Mpc_n{2}_nb{3}_nt{4}'.format(
    run_dir, int(boxsize), nmesh, nbuff, ntile )
zetagfile = '{0}fields/zetag_{1}Mpc_n{2}_nb{3}_nt{4}'.format(
    run_dir, int(boxsize), nmesh, nbuff, ntile )
zetangfile = '{0}fields/zetang_{1}Mpc_n{2}_nb{3}_nt{4}'.format(
    run_dir, int(boxsize), nmesh, nbuff, ntile )
chifile = '{0}fields/chi_{1}Mpc_n{2}_nb{3}_nt{4}'.format(
    run_dir, int(boxsize), nmesh, nbuff, ntile )

# Plot what delta fields we have
if os.path.isfile(deltafile) and os.path.isfile(deltagfile):
    # Read in delta(x), the full Gaussian + non-Gaussian overdensity
    in_delta = open( deltafile, 'rb' )
    delta = np.fromfile(in_delta,dtype=np.float32,count=-1)
    delta = np.reshape(delta, (nlattice,nlattice,nlattice), order='F')
    delta = delta[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
    # Read in delta_G(x), the Gaussian part of delta
    in_deltag = open( deltagfile, 'rb' )
    deltag = np.fromfile(in_deltag,dtype=np.float32,count=-1)
    deltag = np.reshape(deltag, (nlattice,nlattice,nlattice), order='F')
    deltag = deltag[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
    # Plot final delta field
    fig1,axs1 = plot_field_slices( delta, X,Y,Z )
    # Plot just the NG part of delta
    fig3,axs3 = plot_field_slices( delta-deltag, X,Y,Z,
         int(neff/2), int(neff/2), int(neff/2),
         r'$\delta_{L,NG}(x,y,z=0)$',
         r'$\delta_{L,NG}(x,y=0,z)$',
         r'$\delta_{L,NG}(x=0,y,z)$',
         r'$\delta_{L,NG}$')
    del(delta)
    # Plot delta_G(x), the Gaussian part of delta
    fig2,axs2 = plot_field_slices( deltag, X,Y,Z,
         int(neff/2), int(neff/2), int(neff/2),
         r'$\delta_{L,G}(x,y,z=0)$',
         r'$\delta_{L,G}(x,y=0,z)$',
         r'$\delta_{L,G}(x=0,y,z)$',
         r'$\delta_{L,G}$')
    del(deltag)
    # Save figures
    fig1.set_size_inches(14,12)
    fig1.savefig('delta.png',dpi=100)
    fig2.set_size_inches(14,12)
    fig2.savefig('deltag.png',dpi=100)
    fig3.set_size_inches(14,12)
    fig3.savefig('deltang.png',dpi=100)

elif os.path.isfile(deltafile):
    # Read in delta(x), the full Gaussian + non-Gaussian overdensity
    in_delta = open( deltafile, 'rb' )
    delta = np.fromfile(in_delta,dtype=np.float32,count=-1)
    delta = np.reshape(delta, (nlattice,nlattice,nlattice), order='F')
    delta = delta[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
    # Plot delta(x)
    fig1,axs1 = plot_field_slices( delta, X,Y,Z )
    del(delta)
    # Save figures
    fig1.set_size_inches(14,12)
    fig1.savefig('delta.png',dpi=100)

elif os.path.isfile(deltagfile):
    print('Plotting zeta fields')
    # Read in delta_G(x), the Gaussian part of delta
    in_deltag = open( deltagfile, 'rb' )
    deltag = np.fromfile(in_deltag,dtype=np.float32,count=-1)
    deltag = np.reshape(deltag, (nlattice,nlattice,nlattice), order='F')
    deltag = deltag[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
    # Plot delta_G(x)
    fig2,axs2 = plot_field_slices( deltag, X,Y,Z,
         int(neff/2), int(neff/2), int(neff/2),
         r'$\delta_{L,G}(x,y,z=0)$',
         r'$\delta_{L,G}(x,y=0,z)$',
         r'$\delta_{L,G}(x=0,y,z)$',
         r'$\delta_{L,G}$')
    del(deltag)
    # Save figures
    fig2.set_size_inches(14,12)
    fig2.savefig('deltag.png',dpi=100)
    

# If there is a smoothed field, read it in and determine sigma_8
if os.path.isfile(smoothfile):
    # Read in the Gaussian overdensity smoothed at 8 Mpc/h
    in_smooth = open( smoothfile, 'rb' )
    smooth = np.fromfile(in_smooth,dtype=np.float32,count=-1)
    smooth = np.reshape(smooth, (nlattice,nlattice,nlattice), order='F')
    sigma8 = ( np.sum(smooth**2) * nlattice**-3 )**.5
    print('sigma_8=',sigma8,',  mean=',np.mean(smooth),',  (w/ buffers)')
    smooth = smooth[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
    # Plot the Gaussian overdensity smoothed at 8 Mpc/h
    fig4,axs4 = plot_field_slices( smooth, X,Y,Z,
        int(neff/2), int(neff/2), int(neff/2),
        r'$\langle\delta_L(x,y,z=0) | R_{TH}=8 h^{-1}$Mpc$\rangle$',
        r'$\langle\delta_L(x,y=0,z) | R_{TH}=8 h^{-1}$Mpc$\rangle$',
        r'$\langle\delta_L(x=0,y,z) | R_{TH}=8 h^{-1}$Mpc$\rangle$',
        r'$\langle\delta_L | R_{TH}=8 h^{-1}$Mpc$\rangle$' )
    del(smooth)
    # Save figures
    fig4.set_size_inches(14,12)
    fig4.savefig('delta_smooth.png',dpi=100)


# If there's a zeta_G(x) field, plot that
if os.path.isfile(zetagfile) and os.path.isfile(zetangfile):
    print('Plotting zeta fields')
    # Read in zeta_G(x), the Gaussian part of zeta
    in_zetag = open( zetagfile, 'rb')
    zetag = np.fromfile(in_zetag,dtype=np.float32,count=-1)
    zetag = np.reshape(zetag, (nlattice,nlattice,nlattice), order='F')
    zetag = zetag[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
    # Read in zeta(x) (the Gaussian + non-Gaussian part)
    in_zetang = open( zetangfile, 'rb')
    zetang = np.fromfile(in_zetang,dtype=np.float32,count=-1)
    zetang = np.reshape(zetang, (nlattice,nlattice,nlattice), order='F')
    zetang = zetang[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
    # Plot zeta(x)-zeta_G(x), the non-Gaussian part of zeta
    fig5,axs5 = plot_field_slices( zetang-zetag, X,Y,Z,
        int(neff/2), int(neff/2), int(neff/2),
        r'$(\zeta-\zeta_G)(x,y,z=0)$',
        r'$(\zeta-\zeta_G)(x,y=0,z)$',
        r'$(\zeta-\zeta_G)(x=0,y,z)$',
        r'$\zeta-\zeta_G$')
    del(zetang)
    # Plot zeta_G(x), the Gaussian part of zeta
    fig6,axs6 = plot_field_slices( zetag, X,Y,Z,
        int(neff/2), int(neff/2), int(neff/2),
        r'$\zeta_G(x,y,z=0)$',
        r'$\zeta_G(x,y=0,z)$',
        r'$\zeta_G(x=0,y,z)$',
        r'$\zeta_G$')
    del(zetag)
    # Save figures
    fig5.set_size_inches(14,12)
    fig5.savefig('zetang.png',dpi=100)
    fig6.set_size_inches(14,12)
    fig6.savefig('zetag.png',dpi=100)
    
# If there's a chi(x) feld, plot that
if os.path.isfile(chifile):
    print('Plotting chi field')
    # Read in chi(x), the transverse inflaton field
    in_chi = open( chifile, 'rb' )
    chi = np.fromfile(in_chi,dtype=np.float32,count=-1)
    chi = np.reshape(chi, (nlattice,nlattice,nlattice), order='F')
    chi = chi[nbuff:-nbuff,nbuff:-nbuff,nbuff:-nbuff]
    # Plot chi(x)
    fig7,axs7 = plot_field_slices( chi, X,Y,Z,
        int(neff/2), int(neff/2), int(neff/2),
        r'$\chi(x,y,z=0)$',
        r'$\chi(x,y=0,z)$',
        r'$\chi(x=0,y,z)$',
        r'$\chi$')
    del(chi)
    # Save figures
    fig7.set_size_inches(14,12)
    fig7.savefig('chi.png',dpi=100)

###########################################################################


# plt.show()

