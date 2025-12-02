# great-attractors.py
# 
# Nathan J. Carlson
# May 19, 2021
# 
# Optimized to run on Python 3.8.2
# 
# This script plots peak displacements and extrema of the gradient of the 
# initial density field to see the effects of great attractors.

import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
#     python3.8 great-attractors.py <side_length> <buffer_size> \
#     <merged_peak_file> <initial_delta_file>
# 
# For example, this might look like
# 
#     cd <path_to_peak-patch>
#     python3.8 python/halo_plotting/great-attractors.py 256. 40 \
#     runs/n256s256/output/256Mpc_n256_nb40_nt1_merge.pksc.13579 \
#     runs/n256s256/fields/Fvec_256Mpc_n256_nb40_nt1
#
# Note: make sure there are no traling spaces following the '\' (which tell
# terminal to continue on the next line), or the code will mess up. 

if len(sys.argv) != 5:
    print('Input Error\nTo run this script in terminal use `python3.8 grea'
        +'t-attractors.py\n<side_length> <buffer_size> <merged_peakfile> <'
        +'density_field>i` where\n<side_length> is the sidelength of the c'
        +'ubic simulation volume in Mpc\nas a float, <buffer_size> is the '
        +'thickness in cells of the buffer (which\nis a cubic shell) <merg'
        +'ed_peakfile> is the final, merged output from a peak\npatch run '
        +'(see peak-patch/readme.md for more information on generating\nme'
        +'rged peak patch output), and <density_field> is the initial dens'
        +'ity field\noutputed by the peak patch run.')
    # Error message displayed if you pass the wrong number of arguments.
    sys.exit(2)

# Reading from command line prompts
s       = float(sys.argv[1])     # cube sidelength [Mpc]
buff    = int(sys.argv[2])       # buffer thickness [cells]
infile  = open(sys.argv[3],'rb') # pointer to <merged_peak_file>
indelta = open(sys.argv[4],'rb') # pointer to <initial_delta_file>

# Constants
omegam   = 0.25           # CDM density fraction
h        = 0.7            # little h
rhocrit  = 2.775e11*h**2  # critical energy density 3H^2/8piG [Msol Mpc^-3]
rho      = rhocrit*omegam # average CDM density [Msol Mpc^-3]
deltavir = 200            # <something defining virial collapse I think>
outnum   = 11             # number of columns in <merged_peak_file>
G        = 4.517e-48      # gravitational constant [Mpc^3 Msol^-1 s^-2]
# Note here that we use mass units of solar masses Msol, spatial units of
# megaparsecs Mpc, and time units of seconds s. 

################################################################
### Reading in data from command line and defining constants ###
################################################################
# Reading data from <merged_peak_file> pointer
Non      = np.fromfile(infile,dtype=np.int32,count=1)[0]           
RTHmax   = np.fromfile(infile,dtype=np.float32,count=1)[0]
zin      = np.fromfile(infile,dtype=np.float32,count=1)[0]
peakdata = np.fromfile(infile,dtype=np.float32,count=outnum*Non)
peakdata = np.reshape(peakdata,(Non,outnum))
# The merged peak file has three header parameters: `Non`, the number of
# halos found, `RTHmax`, the largest top-hat filter at which a peak was
# found, `zin` initial redshift. The remainder of the file contains a 1D 
# array that we then reshape into an 11-column matrix `peakdata`. Each row
# of this 11-column matrix contains data describing a unique halo.

# Next, for the sake of clarity, we copy the raw halo data from the 11
# columns of `peakdata` to 11 better-labelled NumPy arrays:
x      = peakdata[:,0]  #  
y      = peakdata[:,1]  # x,y,z: components of the final (Eulerian) halo
z      = peakdata[:,2]  #     position vector
vx     = peakdata[:,3]  # 
vy     = peakdata[:,4]  # vx,vy,vz: components of the final (Eulerian) halo
vz     = peakdata[:,5]  #     velocity vector
Rth    = peakdata[:,6]  # Rth: top-hat filter radius at which halo was
xL     = peakdata[:,7]  #     found (roughly the radius of the halo)
yL     = peakdata[:,8]  # xL,yL,zL: components of the initial (Lagrangian)
zL     = peakdata[:,9]  #     halo position vector 
deltah = peakdata[:,10] # deltah: halo overdensity

# Finally, we can determine the approximate mass of each halo based on its
# radius (which is greater than or equal to the largest top-hat radius at
# which it is found to have collased `Rth`) as well as the critical density
# of cold dark matter `rho`:
M = 4*np.pi/3*Rth**3*rho # halo masses [Msol]

# Reading data from <initial_peak_file> pointer
delta = np.fromfile(indelta,dtype=np.float32,count=-1)
n     = int(round((len(delta))**(1/3)))
delta = np.reshape(delta, (n,n,n), order='F')
delta = delta[buff:-buff,buff:-buff,buff:-buff]
n     = n-2*buff


def Marianne(delta):
    # Creates figure object with three subplots in separate columns
    fig, axs = plt.subplots(1)

    # Defines bin edges for real-space and phase-space plots
    edges    = np.linspace( -s/2 , s/2 , n+1 )
    X,Y,Z    = np.meshgrid(edges,edges,edges,indexing='ij')
    XX,YY,ZZ = np.meshgrid(edges[:-1],edges[:-1],edges[:-1],indexing='ij')

    # Plots delta in first subplot
    #ax1 = axs[0]
    p1  = axs.pcolormesh( X[:,:,0], Y[:,:,0], -YY[:,:,n//4] + delta[:,:,n//4],
                          cmap='gist_rainbow' )
    axs.set_title(r'Happy Pride!',fontsize=40)
    #ax1.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    #ax1.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    #ax1.set_aspect('equal')
    #ax1_pos = make_axes_locatable(ax1)
    #ax1_colorbar_pos = ax1_pos.append_axes('right', size='5%', pad=0.05)
    #fig.colorbar(p1, cax=ax1_colorbar_pos)
    #ax1_colorbar_pos.yaxis.set_label_position('right')
    #ax1_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
    #                           labelpad=15, rotation=270, fontsize=14)
    return fig

pride = Marianne(delta)
plt.show()

##########################
### Fourier transforms ###
##########################

# Fourier phase k [Mpc^-1]
kscale     = scipy.fft.fftfreq(n,s/n) 
k_for_plot = scipy.fft.fftshift(kscale) #=np.linspace(-n/s/2, (n-2)/s/2, n)
kscale[np.argwhere(kscale==0)] = np.inf
kX,kY,kZ   = np.meshgrid(kscale,kscale,kscale,indexing='ij')
# The k=(0,0,0) phase component of deltak corresponds to an additive term
# in delta, and the DFT is invariant under additive factors, so we can just
# set this term of the sum over k to zero. Since we take
# iFFT(k^-2 FFT(delta)), this is accomplished by replacing k=0 with k=inf.

# Fourier space overdensity field
deltak          = scipy.fft.fftn(delta)               # FFT(delta)
deltak_for_plot = np.real(scipy.fft.fftshift(deltak)) # shifts domain
inverse_deltak  = np.real(scipy.fft.ifftn(deltak))    # iFFT(FFT(delta))
# For reasons related to symmetries of phase space, the scipy function
# fftfreq places 0 as the first entry of the an array, followed by the
# positive frequencies, followed by the negative frequencies:
# k = fftfreq(n,s/n) = [ 0.   ,  1/s      ,  2/s ,  ...  ,  (n-1)/2s, 
#                        -n/2s,  -(n-1)/2s,  ... ,   -2/s,  -1/s    ])
# this isn't very conveninet if we want to plot it, so for the array
# deltak_for_plot, I use the scipy function fftshift to rearange this numpy
# array into a more intuitive form:
# fftshift(k) = [ -n/2s,  -(n-1)/2s,  ...,  -2/s,  -1/s    ,
#                 0.   ,  1/s      ,  2/s,  ... ,  (n-1)/2s])

# Gravitational potential field Phi [Mpc^2 s^-2]
Phi = scipy.fft.ifftn( deltak / (kX**2+kY**2+kZ**2))
Phi = -4*np.pi * G * rho * np.real(Phi)
# In Lagrangian space, we drop the scale factor a(t)=1

# Restore kX,kY,kZ
kscale[np.argwhere(kscale==np.inf)] = 0.
kX,kY,kZ       = np.meshgrid(kscale,kscale,kscale,indexing='ij')

# Fourier transform of Phi
Phik           = scipy.fft.fftn(Phi)

# delta from Phi0^-1 * iFFT(-k^2 FFT(Phik))
delta_from_Phi = ( (4*np.pi*G*rho)**(-1) *
    np.real(scipy.fft.ifftn( -(kX**2+kY**2+kZ**2) * Phik )) )

# Gradient of Phi given by iFFT(k FFT(Phik)) where k is a vector
#dPhi = np.array([ -1.j * np.real( scipy.fft.ifftn(kX*Phik) ) ,
#                  -1.j * np.real( scipy.fft.ifftn(kY*Phik) ) ,
#                  -1.j * np.real( scipy.fft.ifftn(kZ*Phik) ) ])
#dPhi_mag = np.sqrt(np.abs( dPhi[0]**2 + dPhi[1]**2 + dPhi[2]**2 ))
dPhi = np.gradient( Phi )
dPhi_mag = np.sqrt( dPhi[0]**2 + dPhi[1]**2 + dPhi[2]**2 )


#############################################
### Finding local minima of the Phi field ###
#############################################

Phi_local_minima   = np.array([ [], [], [] ])
Phi_local_min_mags = np.array([])
Phi_local_maxima   = np.array([ [], [], [] ])
Phi_local_max_mags = np.array([])
dPhi_threshold   = np.min(dPhi_mag) + np.max(dPhi_mag)/5

for a in range(0,n):
    for b in range(0,n):
        for c in range(0,n):
            if dPhi_mag[a,b,c] < dPhi_threshold:
                if ( Phi[a,b,c] > Phi[ (a+1)%n, b,       c       ] and
                     Phi[a,b,c] > Phi[ a-1,     b,       c       ] and
                     Phi[a,b,c] > Phi[ a,       (b+1)%n, c       ] and
                     Phi[a,b,c] > Phi[ a,       b-1,     c       ] and
                     Phi[a,b,c] > Phi[ a,       b,       (c+1)%n ] and
                     Phi[a,b,c] > Phi[ a,       b,       c-1     ] ):
                    
                    Phi_local_maxima = np.array([
                        np.append( Phi_local_maxima[0],a ),
                        np.append( Phi_local_maxima[1],b ),
                        np.append( Phi_local_maxima[2],c )])

                    Phi_local_max_mags = np.append(
                        Phi_local_max_mags, Phi[a,b,c] )
                
                elif ( Phi[a,b,c] < Phi[ (a+1)%n, b,       c       ] and
                       Phi[a,b,c] < Phi[ a-1,     b,       c       ] and
                       Phi[a,b,c] < Phi[ a,       (b+1)%n, c       ] and
                       Phi[a,b,c] < Phi[ a,       b-1,     c       ] and
                       Phi[a,b,c] < Phi[ a,       b,       (c+1)%n ] and
                       Phi[a,b,c] < Phi[ a,       b,       c-1     ] ):
                    
                    Phi_local_minima = np.array([
                        np.append( Phi_local_minima[0],a ),
                        np.append( Phi_local_minima[1],b ),
                        np.append( Phi_local_minima[2],c )])

                    Phi_local_min_mags = np.append(
                        Phi_local_min_mags, Phi[a,b,c] )



##################################
### Plotting Peak Patch Fields ###
##################################
# For compactness, I've confined each figure to a function

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


def plot_delta_deltak_Phi_slices(delta,Phi,s,n,k_for_plot,deltak_for_plot,
                                 display=True):
    ''' INPUTS/OUTPUTS:
        All inputs except 'display' are share the class of like-named
        variables in the main code body. If display=True, there is no
        output, and the function simply displays the figure generated. If
        display=False, then 'fig' (of class matplotlib.figure.Figure) is
        returned.
        
        DESCRIPTION:
        Plots an {s/n}-thick slice of the initial density field 'delta',
        its Fourier conjugate 'deltak', and the corresponding gravitational 
        potential field 'Phi' given by the inverse FFT of deltak times k^-2
        (up to constants).
        
        USAGE:
        This function can be used to directly display a plot (as described
        above) to the screen, e.g. with
        >>> plot_delta_deltak_Phi_slices(delta,Phi,s,n,k_for_plot,
            deltak_for_plot)
        
        Or to return a matplotlib.figure.Figure class object, which allows
        us to display multiple figures at once, e.g. with
        >>> fig1 = plot_delta_deltak_Phi_slices(delta,Phi,s,n,k_for_plot,
            deltak_for_plot,display=False)
        >>> plt.show()
    '''
    # Creates figure object with three subplots in separate columns
    fig, axs = plt.subplots(1,3,figsize=(13,5))

    # Defines bin edges for real-space and phase-space plots
    edges    = np.linspace( -s/2 , s/2 , n+1 )
    X,Y,Z    = np.meshgrid(edges,edges,edges,indexing='ij')
    kedges   = np.append( k_for_plot    -(k_for_plot[1]-k_for_plot[0])/2 , 
                          k_for_plot[-1]+(k_for_plot[1]-k_for_plot[0])/2 )
    kx,ky,kz = np.meshgrid(kedges,kedges,kedges,indexing='ij')

    # Plots delta in first subplot
    ax1 = axs[0]
    p1  = ax1.pcolormesh( X[:,:,0], Y[:,:,0], delta[:,:,n//4],
                          cmap='viridis' )
    ax1.set_title(r'$\delta(\mathbf{x})$',fontsize=16)
    ax1.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax1.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax1.set_aspect('equal')
    ax1_pos = make_axes_locatable(ax1)
    ax1_colorbar_pos = ax1_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p1, cax=ax1_colorbar_pos)
    ax1_colorbar_pos.yaxis.set_label_position('right')
    ax1_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
                               labelpad=15, rotation=270, fontsize=14)

    # Plots deltak in second subplot
    ax2 = axs[1]
    p2  = ax2.pcolormesh( kx[:,:,0], ky[:,:,0], deltak_for_plot[:,:,n//4],
                          cmap='viridis' )
    ax2.set_title(r'$\delta(\mathbf{k}) = \mathcal{F}[\delta(\mathbf{x})]$'
                  +r'$(\mathbf{k})$',fontsize=16)
    ax2.set_xlabel(r'$k_x ~ [h \, \mathrm{Mpc}^{-1}$]',fontsize=14)
    ax2.set_ylabel(r'$k_y ~ [h \, \mathrm{Mpc}^{-1}$]',fontsize=14)
    ax2.set_aspect('equal')
    ax2_pos = make_axes_locatable(ax2)
    ax2_colorbar_pos = ax2_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p2, cax=ax2_colorbar_pos)
    ax2_colorbar_pos.yaxis.set_label_position('right')
    ax2_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
                               labelpad=15, rotation=270, fontsize=14)

    # Plots Phi in third subplot
    ax3 = axs[2]
    p3  = ax3.pcolormesh( X[:,:,0], Y[:,:,0], Phi[:,:,n//4],
                          cmap='viridis' )
    ax3.set_title(r'$ \Phi(\mathbf{x}) = -4\pi G \bar{\rho}_m $'
                  +r'$ \mathcal{F}^{-1} $' 
                  +r'$ \left[ \frac{ \delta(\mathbf{k}) }{ k^2} \right] $'
                  +r'$ (\mathbf{x}) $',fontsize=16)
    ax3.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax3.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax3.set_aspect('equal')
    ax3_pos = make_axes_locatable(ax3)
    ax3_colorbar_pos = ax3_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p3, cax=ax3_colorbar_pos)
    ax3_colorbar_pos.yaxis.set_label_position('right')
    ax3_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
                               labelpad=15, rotation=270, fontsize=14)

    # Displays plot or returns figure object
    plt.tight_layout(pad=1., w_pad=0., h_pad=0.)
    if display==True:
        plt.show()
    else:
        return fig


def check_Phi(delta,Phi,delta_from_Phi,s,n,display=True):
    ''' INPUTS/OUTPUTS:
        All inputs except 'display' are share the class of like-named
        variables in the main code body. If display=True, there is no
        output, and the function simply displays the figure generated. If
        display=False, then 'fig' (of class matplotlib.figure.Figure) is
        returned.
        
        DESCRIPTION:
        Plots an {s/n}-thick slice of the initial density field 'delta',
        the corresponding gravitational potential field 'Phi', and the
        density field up to a numerical uncertainty found by taking
        iFFT(k^2 FFT(Phi)) times a constant.
        
        USAGE:
        This function can be used to directly display a plot (as described
        above) to the screen, e.g. with
        >>> check_Phi(delta,Phi,delta_from_Phi,s,n)
        
        Or to return a matplotlib.figure.Figure class object, which allows
        us to display multiple figures at once, e.g. with
        >>> fig2 = check_Phi(delta,Phi,delta_from_Phi,s,n,display=False)
        >>> plt.show()
    '''
    # Creates figure object with three subplots in separate columns
    fig, axs = plt.subplots(1,3,figsize=(13,5))
   
    # Defines bin edges for plots
    edges    = np.linspace( -s/2 , s/2 , n+1 )
    X,Y,Z    = np.meshgrid(edges,edges,edges,indexing='ij')

    # Plots delta in first subplot
    ax1 = axs[0]
    p1  = ax1.pcolormesh( X[:,:,0], Y[:,:,0], delta[:,:,n//4],
                          cmap='viridis' )
    ax1.set_title(r'$\delta(\mathbf{x})$',fontsize=16)
    ax1.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax1.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax1.set_aspect('equal')
    ax1_pos = make_axes_locatable(ax1)
    ax1_colorbar_pos = ax1_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p1, cax=ax1_colorbar_pos)
    #ax1_colorbar_pos.yaxis.set_label_position('right')
    #ax1_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
    #                           labelpad=15, rotation=270, fontsize=14)

    # Plots deltak in second subplot
    ax2 = axs[1]
    p2  = ax2.pcolormesh( X[:,:,0], Y[:,:,0], Phi[:,:,n//4],
                          cmap='viridis' )
    ax2.set_title(r'$ \Phi(\mathbf{x}) = \Phi_0 $'
                  +r'$ \mathcal{F}^{-1} $'
                  +r'$ \left[ \frac{ \delta(\mathbf{k}) }{ k^2} \right] $'
                  +r'$ (\mathbf{x}) $',fontsize=16)
    ax2.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    #ax2.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax2.set_aspect('equal')
    ax2_pos = make_axes_locatable(ax2)
    ax2_colorbar_pos = ax2_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p2, cax=ax2_colorbar_pos)
    #ax2_colorbar_pos.yaxis.set_label_position('right')
    #ax2_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
    #                           labelpad=15, rotation=270, fontsize=14)

    # Plots Phi in third subplot
    ax3 = axs[2]
    p3  = ax3.pcolormesh( X[:,:,0], Y[:,:,0], delta_from_Phi[:,:,n//4],
                          cmap='viridis' )
    ax3.set_title(r'$ \delta(\mathbf{x}) = \Phi_0^{-1} \mathcal{F}^{-1} $'
                  +r'$ [ k^2 \Phi(\mathbf{k}) ] (\mathbf{x}) $',
                  fontsize=16)
    ax3.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    #ax3.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax3.set_aspect('equal')
    ax3_pos = make_axes_locatable(ax3)
    ax3_colorbar_pos = ax3_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p3, cax=ax3_colorbar_pos)
    ax3_colorbar_pos.yaxis.set_label_position('right')
    ax3_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
                               labelpad=15, rotation=270, fontsize=14)

    # Displays plot or returns figure object
    plt.tight_layout(pad=0., w_pad=0., h_pad=0.)
    if display==True:
        plt.show()
    else:
        return fig


def plot_delta_Phi_slices(delta,Phi,s,n,display=True):
    ''' INPUTS/OUTPUTS:
        All inputs except 'display' are share the class of like-named
        variables in the main code body. If display=True, there is no
        output, and the function simply displays the figure generated. If
        display=False, then 'fig' (of class matplotlib.figure.Figure) is
        returned.
        
        DESCRIPTION:
        Plots an {s/n}-thick slice in the x-y plane of the initial density
        field 'delta', and the corresponding gravitational potential field
        'Phi' at a series of points in the z axis. The 'delta' fields are
        subplots in the first row and the 'Phi' fields are subplots in the
        second row.
        
        USAGE:
        This function can be used to directly display a plot (as described
        above) to the screen, e.g. with
        >>> plot_delta_Phi_slices(delta,Phi,s,n)
        
        Or to return a matplotlib.figure.Figure class object, which allows
        us to display multiple figures at once, e.g. with
        >>> fig3 = plot_delta_Phi_slices(delta,Phi,s,n,display=False) 
        >>> plt.show()
    '''
    # Creates figure object with three subplots in separate columns
    fig, axs = plt.subplots(2,3,figsize=(13,10))
    
    # Defines bin edges for plots
    edges    = np.linspace( -s/2 , s/2 , n+1 )
    X,Y,Z    = np.meshgrid(edges,edges,edges,indexing='ij')
    
    # Plots delta(x,y,z=-s/4)
    ax1 = axs[0,0]
    p1  = ax1.pcolormesh( X[:,:,0], Y[:,:,0], delta[:,:,n//4],
                          cmap='viridis' )
    ax1.set_title(r'$\delta(x,y,z=$'+str(-s/4)+r'$\,\rm{Mpc})$',
                  fontsize=16)
    ax1.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax1.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax1.set_aspect('equal')
    ax1_pos = make_axes_locatable(ax1)
    ax1_colorbar_pos = ax1_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p1, cax=ax1_colorbar_pos)
    #ax1_colorbar_pos.yaxis.set_label_position('right')
    #ax1_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
    #                           labelpad=15, rotation=270, fontsize=14)

    # Plots delta(x,y,z=0)
    ax2 = axs[0,1]
    p2  = ax2.pcolormesh( X[:,:,0], Y[:,:,0], delta[:,:,n//2],
                          cmap='viridis' )
    ax2.set_title(r'$\delta(x,y,z=0 \, \rm{ Mpc})$',fontsize=16)
    ax2.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    #ax2.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax2.set_aspect('equal')
    ax2_pos = make_axes_locatable(ax2)
    ax2_colorbar_pos = ax2_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p2, cax=ax2_colorbar_pos)
    #ax2_colorbar_pos.yaxis.set_label_position('right')
    #ax2_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
    #                           labelpad=15, rotation=270, fontsize=14)

    # Plots delta(x,y,z=+s/4)
    ax3 = axs[0,2]
    p3  = ax3.pcolormesh( X[:,:,0], Y[:,:,0], delta[:,:,3*n//4],
                          cmap='viridis' )
    ax3.set_title(r'$\delta(x,y,z=$'+str(s/4)+r'$\,\rm{Mpc})$',fontsize=16)
    ax3.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    #ax3.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax3.set_aspect('equal')
    ax3_pos = make_axes_locatable(ax3)
    ax3_colorbar_pos = ax3_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p3, cax=ax3_colorbar_pos)
    ax3_colorbar_pos.yaxis.set_label_position('right')
    ax3_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
                               labelpad=15, rotation=270, fontsize=14)
    
    # Plots Phi(x,y,z=-s/4)
    ax4 = axs[1,0]
    p4  = ax4.pcolormesh( X[:,:,0], Y[:,:,0], Phi[:,:,n//4],
                          cmap='viridis' )
    ax4.set_title(r'$\Phi(x,y,z=$'+str(-s/4)+r'$\,\rm{Mpc})$',fontsize=14)
    ax4.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax4.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax4.set_aspect('equal')
    ax4_pos = make_axes_locatable(ax4)
    ax4_colorbar_pos = ax4_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p4, cax=ax4_colorbar_pos)
    #ax4_colorbar_pos.yaxis.set_label_position('right')
    #ax4_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
    #                           labelpad=15, rotation=270, fontsize=14)
    
    # Plots Phi(x,y,z=0)
    ax5 = axs[1,1]
    p5  = ax5.pcolormesh( X[:,:,0], Y[:,:,0], Phi[:,:,n//2],
                          cmap='viridis' )
    ax5.set_title(r'$\Phi(x,y,z=0 \, \rm{ Mpc})$',fontsize=16)
    ax5.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    #ax5.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax5.set_aspect('equal')
    ax5_pos = make_axes_locatable(ax5)
    ax5_colorbar_pos = ax5_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p5, cax=ax5_colorbar_pos)
    #ax5_colorbar_pos.yaxis.set_label_position('right')
    #ax5_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
    #                           labelpad=15, rotation=270, fontsize=14)

    # Plots Phi(x,y,z=+s/4)
    ax6 = axs[1,2]
    p6  = ax6.pcolormesh( X[:,:,0], Y[:,:,0], Phi[:,:,3*n//4],
                          cmap='viridis' )
    ax6.set_title(r'$\Phi(x,y,z=$'+str(s/4)+r'$\, \rm{ Mpc})$',fontsize=16)
    ax6.set_xlabel(r'$x ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    #ax6.set_ylabel(r'$y ~ [h^{-1} \mathrm{Mpc}]$',fontsize=14)
    ax6.set_aspect('equal')
    ax6_pos = make_axes_locatable(ax6)
    ax6_colorbar_pos = ax6_pos.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(p6, cax=ax6_colorbar_pos)
    ax6_colorbar_pos.yaxis.set_label_position('right')
    ax6_colorbar_pos.set_ylabel(r'$[\mathrm{Mpc}^2 \, \mathrm{s}^{-2}]$',
                               labelpad=15, rotation=270, fontsize=14)

    # Displays plot or returns figure object
    plt.tight_layout(pad=0., w_pad=0., h_pad=0.)
    if display==True:
        plt.show()
    else:
        return fig


def plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,PhiL,
                      xspan=[0,n],yspan=[0,n],zspan=[n//2,n//2+1],
                      display=True):
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

    def seek(xL,yL,zL,Rth,w):
        ''' NumPy arrays xL,yL,zL and Rth contain 3 lagrangian position
            coordinates and approximate radius of each peak patch respect-
            ively. The first element of each array define one peak, the
            second element of each define a second peak, etc. For this
            plot, we are only interested in a subset of these peaks that
            lie within the volume slice given by 
                -s/2 + s/n*xspan[0] < xL < -s/2 + s/n*xspan[1],
                -s/2 + s/n*yspan[0] < yL < -s/2 + s/n*yspan[1],
                -s/2 + s/n*zspan[0] < zL < -s/2 + s/n*zspan[1].
            To find the relevant peaks, we use a series of calls to the
            function np.argwhere(), i.e. np.argwhere(xL>0) returns an array
            of the elements of xL for which xL>0, and passing 
            xL[np.argwhere(xL>0)] returns those elements of xL as a column
            array. This function then quickly reshapes those column arrays
            into row arrays.
        '''
        return ( xL[w].reshape(len(w)) , yL[w].reshape(len(w))  ,
                 zL[w].reshape(len(w)) , Rth[w].reshape(len(w)) )

    def sek(Px,Py,Pz,w):
        ''' Here we do essentailly the same procedure as the former seek(),
            but now instead of determining the peak patches that lie within
            our volume slice, we are finding the local maxima/minima of Phi
            that lie within this volume slice.
        '''
        return ( Px[w].reshape(len(w)) , Py[w].reshape(len(w)) ,
                 Pz[w].reshape(len(w)) )

    ## Dealing with edge effects
    #
    # This was an added feature that I was pursuing but ended up giving up
    # on. The idea is that if you want to display the Phi field slice at 
    # the boundaries then you end up having half as many peak patches
    # because they get cut off. Because the volume has periodic boundary
    # conditions, you could pretty easily take the peaks at the other
    # extreme and copy them to have ie positions >s/2 or <-s/2.
    #
    #if xspan[0]<0:
    #    w=np.argwhere(xL>-s/2+s/n*(n+xspan[0])); e=seek(xL,yL,zL,Rth,w)
    #    xL,yL,zL,Rth=(np.append(xL,e[0]-s),
    #        np.append(yL,e[1]),np.append(zL,e[2]),np.append(Rth,e[3]))
    #    w=np.argwhere(Px>n+xspan[0]); f=sek(Px,Py,Pz,w)
    #    Px,Py,Pz=(np.append(Px,f[0]-n),
    #        np.append(Py,f[1]),np.append(Pz,f[2]))
    #if xspan[1]>=n:
    #    w=np.argwhere(xL<-s/2+s/n*(xspan[1]-n)); e=seek(xL,yL,zL,Rth,w)
    #    xL,yL,zL,Rth=(np.append(xL,e[1]+s),
    #        np.append(yL,e[0]),np.append(zL,e[2]),np.append(Rth,e[3]))
    #    w=np.argwhere(Px<xspan[0]-n); f=sek(Px,Py,Pz,w)
    #    Px,Py,Pz=(np.append(Px,f[0]-n),
    #        np.append(Py,f[1]),np.append(Pz,f[2]))
    #if yspan[0]<0:
    #    w=np.argwhere(yL>-s/2+s/n*(n+yspan[0])); e=seek(xL,yL,zL,Rth,w)
    #    xL,yL,zL,Rth=(np.append(yL,e[2]-s),
    #        np.append(xL,e[0]),np.append(zL,e[1]),np.append(Rth,edge[3]))
    #if yspan[1]>=n:
    #    w=np.argwhere(yL<-s/2+s/n*(yspan[1]-n)); e=seek(xL,yL,zL,Rth,w)
    #    xL,yL,zL,Rth=(np.append(yL,e[1]+s),
    #        np.append(xL,e[0]),np.append(zL,e[2]),np.append(Rth,e[3]))
    #if zspan[0]<0:
    #    w=np.argwhere(zL>-s/2+s/n*(n+zspan[0])); e=seek(xL,yL,zL,Rth,w)
    #    xL,yL,zL,Rth=(np.append(zL,e[0]-s),
    #        np.append(xL,e[0]),np.append(yL,e[1]),np.append(Rth,e[3]))
    #if zspan[1]>=n:
    #    w=np.argwhere(zL<-s/2+s/n*(zspan[1]-n)); e=seek(xL,yL,zL,Rth,w)
    #    xL,yL,zL,Rth=(np.append(zL,e[1]+s),
    #        np.append(xL,e[0]),np.append(yL,e[1]),np.append(Rth,e[3])) 
    
    # Trimming
    w=np.argwhere(xL>=-s/2+s/n*xspan[0]); xL,yL,zL,Rth=seek(xL,yL,zL,Rth,w)
    w=np.argwhere(xL< -s/2+s/n*xspan[1]); xL,yL,zL,Rth=seek(xL,yL,zL,Rth,w)
    w=np.argwhere(yL>=-s/2+s/n*yspan[0]); xL,yL,zL,Rth=seek(xL,yL,zL,Rth,w)
    w=np.argwhere(yL< -s/2+s/n*yspan[1]); xL,yL,zL,Rth=seek(xL,yL,zL,Rth,w)
    w=np.argwhere(zL>=-s/2+s/n*zspan[0]); xL,yL,zL,Rth=seek(xL,yL,zL,Rth,w)
    w=np.argwhere(zL< -s/2+s/n*zspan[1]); xL,yL,zL,Rth=seek(xL,yL,zL,Rth,w)
    Px,Py,Pz = PhiL[0],PhiL[1],PhiL[2]       ; del(PhiL)
    w        = np.argwhere(Px >= xspan[0])   ; Px,Py,Pz = sek(Px,Py,Pz,w)
    w        = np.argwhere(Px <  xspan[1])   ; Px,Py,Pz = sek(Px,Py,Pz,w)
    w        = np.argwhere(Py >= yspan[0])   ; Px,Py,Pz = sek(Px,Py,Pz,w)
    w        = np.argwhere(Py <  yspan[1])   ; Px,Py,Pz = sek(Px,Py,Pz,w)
    w        = np.argwhere(Pz >= zspan[0])   ; Px,Py,Pz = sek(Px,Py,Pz,w)
    w        = np.argwhere(Pz <  zspan[1])   ; Px,Py,Pz = sek(Px,Py,Pz,w)
    PhiL     = -s/2+s/n*np.array([Px,Py,Pz]) ; del(Px,Py,Pz,w)
    
    # Dependent variable used to plot circular peak patches
    theta = np.linspace(0, 2*np.pi, 100, endpoint=False)
    
    # Define bin edges for Phi plot
    edges    = np.linspace( -s/2 , s/2 , n+1 )
    X,Y,Z    = np.meshgrid(edges,edges,edges,indexing='ij')
    
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

    # Plots peak patches as circles centred at xL,yL,zL with radius Rth
    for a in range(len(xL)):
        mpp, = axs.plot(Rth[a]*np.cos(theta)+xL[a], Rth[a]*np.sin(theta)+yL[a], 
                        color='white')

    # Plots minima of the Phi field    
    axs.plot( PhiL[0], PhiL[1], ls='none', marker='+', mfc='white', 
              mec='white')#, label=r'Minima of $\Phi(\mathbf{x})$' )

    # Legend
    axs.legend(loc='center', bbox_to_anchor=(0., 1., 1., .1), ncol=2, #facecolor='light grey', 
        title=(r'For a slice '+legend_string+' projected onto the '+
            the_plane+r' plane'),
        handles=(
            Line2D([0],[0],ls='none',marker='o',mfc='none',mec='k',
                markersize=10,label=r'Peak patches'),
            Line2D([0],[0],ls='none',marker='+',mec='k',
                label=r'Minima of $\Phi(\mathbf{x})$')
            ),
        fancybox=True)

    # Displays plot or returns figure object
    if display==True:
         plt.show()
    else:
         return fig


def plot_peaks_hist(Phi_min,Phi_max,display=True):
    ''' INPUTS/OUTPUTS:
        'Phi_min' is a NumPy array containing the values of Phi at each
        minima located in a sequential nearest-neighbours search, 'Phi_max'
        is a similar array for the values of Phi at each maxima. If
        display=False, then 'fig' (of the class 'matplotlib.figure.Figure')
        is returned. 
        
        DESCRIPTION:
        This function plots a histogram showing the distribution of local
        maxima and minimia of Phi as a function of their value in Phi.

        USAGE:
        This function can be used to directly display a plot (as described
        above) to the screen, e.g. with
        >>> plot_peaks_hist(Phi_local_minima,Phi_local_maxima)
        
        Or to return a matplotlib.figure.Figure class object, which allows
        us to display multiple figures at once, e.g. with
        >>> figure4 = plot_peaks_hist(Phi_local_minima,Phi_local_maxima,
                display=False)
        >>> plt.show()
    '''
    # Convert Phi from Mpc^2 s^-2 to km^2 s^-2
    Phi_min = 9.5214061369e38*Phi_min
    Phi_max = 9.5214061369e38*Phi_max

    fig,ax = plt.subplots(1)#, figsize=(12,10))

    bin1 = max( np.abs( min(np.min(Phi_min),np.min(Phi_max)) ) ,
                np.abs( max(np.max(Phi_min),np.max(Phi_max)) ) )
    bins = np.linspace(-bin1,bin1,50)

    ax.hist(Phi_min, bins=bins, alpha=.5, label=r'$\Phi$ minima')
    ax.hist(Phi_max, bins=bins, alpha=.5, label=r'$\Phi$ maxima')
    ax.set_xlabel(r'$\Phi$ bin centres $[\mathrm{km}^2~\mathrm{s}^{-2}]$')
    ax.set_ylabel(r'Counts')
    ax.legend()
    
    # Displays plot or returns figure object
    if display==True:
         plt.show()
    else:
         return fig



###############################
### Plotting Function Calls ###
###############################

# # Plot slice of delta, deltak and Phi
# fig1 = plot_delta_deltak_Phi_slices(delta,Phi,s,n,k_for_plot,
#                                     deltak_for_plot,display=False)

# # Plot slice of delta, Phi and delta from iFFT(k^2 Phi)
# fig2 = check_Phi(delta,Phi,delta_from_Phi,s,n,display=False)

# Plot various slices of delta and Phi
fig3 = plot_delta_Phi_slices(delta,Phi,s,n,display=False)

# Plot various slices of Phi and grad Phi
fig4 = plot_delta_Phi_slices(Phi,dPhi_mag,s,n,display=False)

# Plot a slice of Phi with corresponding minima and peak patches
#fig5 = plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,PhiL=Phi_local_minima,
#                      xspan=[0,n],yspan=[0,n],zspan=[n//2-3,n//2+3],
#                      display=False)

fig6 = plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,PhiL=Phi_local_minima,
                      xspan=[0,n],yspan=[0,n],zspan=[n//4-3,n//4+3],
                      display=False)

#fig7 = plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,PhiL=Phi_local_minima,
#                      xspan=[0,n],yspan=[0,n],zspan=[3*n//4-3,3*n//4+3],
#                      display=False)

fig8 = plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,PhiL=Phi_local_minima,
                      xspan=[0,n],yspan=[0,n],zspan=[n//4-9,n//4-3],
                      display=False)

fig9 = plot_peaks_maxPhi(xL,yL,zL,Rth,Phi,s,n,PhiL=Phi_local_minima,
                      xspan=[0,n],yspan=[0,n],zspan=[n//4+3,n//4+9],
                      display=False)


fig10 = plot_peaks_hist(Phi_local_min_mags, Phi_local_max_mags,
                      display=False)

# Display matplotlib figures
plt.show()

# # Plot the peak pathces in the simmulation volume
# maya1 = plot_peaks(xL,yL,zL,Rth,display=False)

# # Display Mayavi figures
# mlab.show()

'''
mlab.figure(1, bgcolor=(1,1,1), fgcolor=(.5,.5,.5))
mlab.imshow(xLscale, xLscale, delta[n//4,:,:], colormap='viridis')
mlab.gcf().scene.parallel_projection = True
mlab.colorbar(orientation='vertical')
mlab.axes(nb_labels=4, xlabel=r'x', ylabel='[Mpc]', z_axis_visibility=False)
mlab.view(azimuth=0., elevation=0.)

mlab.figure(2, bgcolor=(1,1,1), fgcolor=(.5,.5,.5))
mlab.imshow(kscale, kscale, deltak[n//4,:,:], colormap='viridis')
mlab.gcf().scene.parallel_projection = True
mlab.colorbar(orientation='vertical')
mlab.axes(nb_labels=4, xlabel='x [Mpc^-1]', ylabel='y [Mpc^-1]', z_axis_visibility=False)
mlab.view(azimuth=0., elevation=0.)

mlab.figure(3, bgcolor=(1,1,1), fgcolor=(.5,.5,.5))
mlab.imshow(xLscale, xLscale, delta2[n//4,:,:], colormap='viridis')
mlab.gcf().scene.parallel_projection = True
mlab.colorbar(orientation='vertical')
mlab.axes(nb_labels=4, xlabel=r'x', ylabel='[Mpc]', z_axis_visibility=False)
mlab.view(azimuth=0., elevation=0.)

mlab.show()
'''

