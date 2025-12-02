# USAGE: at command line, run
# 


import numpy as np
import matplotlib.pyplot as plt
import os,sys
plt.rcParams['text.usetex']=True

# Read halo catalogue file from command line
cat_file = sys.argv[1]
in_cat   = open( cat_file, 'rb' )

# Read catalogue file header
N_halos = np.fromfile( in_cat, dtype=np.int32  , count=1 )[0]
RTHmax  = np.fromfile( in_cat, dtype=np.float32, count=1 )[0]
zin     = np.fromfile( in_cat, dtype=np.float32, count=1 )[0]

# Determine number of columns in halo catalogue table
size_bytes = os.path.getsize(cat_file)
cols       = int( (size_bytes-12)/4/N_halos )

# Read halo catalogue
catalogue = np.reshape( np.fromfile( in_cat, dtype=np.float32, count=N_halos*cols ), (N_halos,cols) )

# Halo positions
x = catalogue[:,0]
y = catalogue[:,1]
z = catalogue[:,2]

#zmin = np.min(z)
#zmax = np.max(z)
#z    -= zmin
#zmax -= zmin
#where = np.argwhere(z<=zmax/5 )
#x=x[where]
#y=y[where]

# Make a scatter plot of halos
fig,ax=plt.subplots()
ax.scatter( x, y, s=.1, alpha=.01, c='k', marker='.' )
ax.set_xlabel(r'comoving distance $x$ [Mpc]')
ax.set_ylabel(r'comoving distance $y$ [Mpc]')
fig.savefig(cat_file+'.png', dpi=200,bbox_inches='tight')
#fig.savefig(cat_file+'.pdf',bbox_inches='tight')
