import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from mayavi import mlab
plt.rcParams['text.usetex'] = True

def dV(phi,chi):
    return (phi**2-1)**2*((chi**2-1)**2-1)/4

# MatPlotLib Plot
fig,ax = plt.subplots(subplot_kw={"projection": "3d"})

X,Y = np.arange(-1,1,.001),np.arange(-2,2,.002)
X,Y = np.meshgrid(X,Y)
Z   = dV(X,Y)
mpl_surf = ax.plot_surface(X,Y,Z,cmap=cm.coolwarm,linewidth=0,antialiased=False)
ax.set_xlabel(r'$(\phi-\phi_p)/\phi_w$')
ax.set_ylabel(r'$\chi/\nu$')
ax.set_zlabel(r'$\frac{\Delta V(\phi,\chi)}{\nu^4\lambda_\chi}$')
ax.zaxis.set_rotate_label(False)
fig.colorbar(mpl_surf, shrink=0.67, aspect=10)
plt.show()

# MayaVi plot
x,y = np.mgrid[ -1:1:0.01 , -2:2:0.02 ]
s = mlab.surf( x,y,dV )
mlab.axes()
mlab.show()
