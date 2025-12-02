# USAGE: at command line, pass the raw file name 
# python raw_halo_scatter.py <...>/<run_name>_raw.pksc.<seed>
# I will then loop over the ntile^3 individual raw files from the largerun and plot all as a scatter
# 
# <...>/<run_name>_raw.pksc.<seed>_000001
# ...
# <...>/<run_name>_raw.pksc.<seed>_<ntile^3}


import numpy as np
import matplotlib.pyplot as plt
import os,sys
plt.rcParams['text.usetex']=True

raw_file = sys.argv[1]

if not os.path.exists(raw_file+'_000001'):
    sys.exit('This script is meant for plotting scatters of halos for runs'
            +'with\n`largerun=1` in the parameter file. To plot largerun=0'
            +'catalogues, use the\nscript `halo_scatter.py`.')

# Count the number of raw files
num_raw_files = 1
while os.path.exists(raw_file+'_%06d'%(num_raw_files+1)):
    num_raw_files+=1

# Make plot object
fig,ax=plt.subplots()

# Loop over raw files for each parallelization tile
for j in range(1,num_raw_files+1):

    # Read halo catalogue file from command line
    cat_file = raw_file+'_%06d'%j
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
    x,y,z = catalogue[:,0:3].T
    del(cat_file,in_cat,catalogue)

    # Report
    print('Tile {0}/{1}, {2} halos read, RTHmax={3}, zin={4}'.format(j,num_raw_files,N_halos,RTHmax,zin))

    #zmin = np.min(z)
    #zmax = np.max(z)
    #z    -= zmin
    #zmax -= zmin
    #where = np.argwhere(z<=zmax/5 )
    #x=x[where]
    #y=y[where]

    # Make a scatter plot of halos
    ax.scatter( x, y, s=.1, alpha=.01, c='k', marker='.' )

ax.set_xlabel(r'comoving distance $x$ [Mpc]')
ax.set_ylabel(r'comoving distance $y$ [Mpc]')
fig.savefig(raw_file+'_largerun.png', dpi=200,bbox_inches='tight')
#fig.savefig(cat_file+'.pdf',bbox_inches='tight')

