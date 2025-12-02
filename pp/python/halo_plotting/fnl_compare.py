import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import signal
from mayavi import mlab # Mayavi free 3D plotting module


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



if len(sys.argv) != 2:
    # Error message displayed if you pass the wrong number of arguments.
    print('You\'ve passed me the wrong number of arguments at the command '
        +'line.\nCommand line should take the form\ncd <peak_patch_directo'
        +'ry>\npython3.8 python/halo_plotting/plot-great-attract.py runs/n'
        +'s128')
    sys.exit(2)

fnl_range  = [ -1e5 , -1e4 , -1e3 , -1e2 , -1e1 , 1 , 0 , 1 , 1e1 , 1e2 , 1e3 , 1e4 , 1e5 ]
fnl_label  = ['-1e5','-1e4','-1e3','-1e2','-1e1','1','0','1','1e1','1e2','1e3','1e4','1e5']
run_dir    = [''    ,''    ,''    ,''    ,''    ,'' ,'' ,'' ,''   ,''   ,''   ,''   ,''   ]
catalogue_files = run_dir

Non    = np.zeros( len(run_dir) )
RTHmax = np.zeros( len(run_dir) )
zin    = np.zeros( len(run_dir) )




runs_dir = str(sys.argv[1])
if runs_dir[-1] != '/':
    runs_dir += '/'

for i in range(len(fnl_label)):
    run_dir[i] = runs_dir+'fnl='+fnl_label[i]+'/'

with open(run_dir[0]+'param/param.params') as f:
    params = [i.strip() for i in f.readlines()]

for line in params:
    exec(line)

for i in range(len(Non)):
    catalogue_files[i] = '{0}output/{1}_merge.pksc.{2}'.format( run_dir[i], run_name, seed )

for i in range(len(Non)):
    in_catalgoue = open( catalogue_files[i], 'rb' )
    Non[i]    = np.fromfile(in_catalogue[i],dtype=np.int32,count=1)[0]
    RTHmax[i] = np.fromfile(in_catalogue[i],dtype=np.int32,count=1)[0]
    zin[i]    = np.fromfile(in_catalogue[i],dtype=np.int32,count=1)[0]
    catalogue = np.fromfile(in_catalogue[i],dtype=np.float32,count=33*Non[i])
    catalogue = np.reshape(catalogue,(Non[i],33))

    xL     = catalogue[:,7]
    yL     = catalogue[:,8]
    zL     = catalogue[:,9]
    x      = catalogue[:,0]
    y      = catalogue[:,1]
    z      = catalogue[:,2]
    Rth    = catalogue[:,6]
    Fpk    = catalogue[:,10]
    e_v    = catalogue[:,11]
    p_v    = catalogue[:,12]
    # Do some shit
M = 4*np.pi/3*Rth**3*rho
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

    del(catalogue)

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


