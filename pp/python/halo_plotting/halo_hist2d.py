# USAGE: at command line, run
# 
#     python halo_hist2d.py <catalogue_file> <axis>=0 <bins>=1000
#         <mass-weighted>=0 <mass_cut>='none' <mass_cutoff>=1e11
#         <xmin>=-inf <xmax>=inf
#
# where
# - <catalogue_file> is the halo catalogue file
# - <axis> is the axis projected along to make a 2D histogram (using
#       python numbering so allowed values are 0, 1 or 2, default 0)
# - <bins> is the square root of the number of bins (default 1000)
# - <mass-weighted> set to 1 to weight halos by mass (default 0)
# - <mass_cut> is 'low' to take only halos below <mass_cutoff>, 'high' to
# -     take only halos above, or 'none' to plot all halos (default 'none')
# - <mass_cuttoff> is cutoff mass in solar masses
# - <xmin> in Mpc; if greater than -inf, halos with position on the <axis>
#       axis less than <xmin> are ignored
# - <xmax> in Mpc; if less than inf, halos with poistion on the <axis> axis
#       greater than <xmax> are ignored
# 
# The plot is saved as <catalogue_file>_hist.pdf (in the same directory as
# <catalogue_file>. 

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os,sys
plt.rcParams['text.usetex']=True

# Read halo catalogue file from command line
cat_file = sys.argv[1]
in_cat   = open( cat_file, 'rb' )

# Set number of bins in hist
axis,bins,mass_weighted,mass_cut,mass_cutoff,xmin,xmax = 0,1000,0,'none',1.e11,-np.inf,np.inf
if len(sys.argv)>2: axis          = int(        sys.argv[2] )
if len(sys.argv)>3: bins          = int(        sys.argv[3] )
if len(sys.argv)>4: mass_weighted = int(        sys.argv[4] )
if len(sys.argv)>5: mass_cut      =             sys.argv[5]
if len(sys.argv)>6: mass_cutoff   = float(      sys.argv[6] )
if len(sys.argv)>7: xmin          = np.float32( sys.argv[7] )
if len(sys.argv)>8: xmax          = np.float32( sys.argv[8] )

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
if   axis==0: x,y = catalogue[:,1:3  ].T
elif axis==1: x,y = catalogue[:,0:3:2].T
else        : x,y = catalogue[:,0:2  ].T

# Halo masses
h       = .6735
Omega_m = .3138
rhocrit = 2.775e11*h**2
rho_m   = Omega_m * rhocrit
Rth     = catalogue[:,6]
M       = 4*np.pi/3*Rth**3*rho_m

# sliceing
if xmin > -np.inf or xmax < np.inf: z = catalogue[:,axis]
if xmin > -np.inf:
    scut = np.where( z > xmin )
    x,y,z,M = x[scut],y[scut],z[scut],M[scut]
if xmax < np.inf:
    scut = np.where( z < 100. )
    x,y,z,M = x[scut],y[scut],z[scut],M[scut]

# Cut halos
if mass_cut=='low':
    Mcut = np.where( M<mass_cutoff )
    x, y, M = x[Mcut], y[Mcut], M[Mcut]
elif mass_cut=='high':
    Mcut = np.where( M>mass_cutoff )
    x, y, M = x[Mcut], y[Mcut], M[Mcut]

# Make 2D histogram
fig,ax = plt.subplots()

#if mass_weighted: hist = ax.hist2d(x,y,bins=bins,range=[[-1250.,1250.],[6350.,7850.]],weights=M,norm=mpl.colors.LogNorm(vmin=1.e10,vmax=1.e14))
#else:             hist = ax.hist2d(x,y,bins=bins,range=[[-1250.,1250.],[6350.,7850.]])
if mass_weighted: hist = ax.hist2d(x,y,bins=bins,weights=M,norm=mpl.colors.LogNorm(vmin=1.e10,vmax=1.e14))
else:             hist = ax.hist2d(x,y,bins=bins,cmap='binary')

# Make colourbar and labels
cbar = fig.colorbar( hist[3], ax=ax )
if mass_weighted: cbar.set_label(r'projected mass $[M_\odot]$')
else:             cbar.set_label(r'projected halo count')
ax.set_xlabel(r'comoving distance [Mpc]')
ax.set_ylabel(r'comoving distance [Mpc]')

# Convert to scientific notation x = a * 10 ** b
def scinot(x):
    b = int(np.floor(np.log10(np.abs(x))))
    a = x/10**b
    return (a,b)

# Add a title noting mass cut if one is made
mc = scinot(mass_cutoff)
if mc[0]==1            : title = r'10^{{{:0}}}'.format(mc[1])
elif mc[0]==int(mc[0]) : title = r'{:0}\times10^{{{:0}}}'.format(mc[0],mc[1])
else                   : title = r'{:0.4}\times10^{{{:0}}}'.format(mc[0],mc[1])
if mass_cut=='low'   : ax.set_title(r'$M<{0}M_\odot$'.format(title))
if mass_cut=='high'  : ax.set_title(r'$M>{0}M_\odot$'.format(title))

# Write histogram to file
fileout = cat_file
if   axis==1         : fileout += '_xz'
elif axis==2         : fileout += '_xy'
if mass_weighted     : fileout += '_massweighted'
if   mass_cut=='high': fileout += '_highmass'
elif mass_cut=='low' : fileout += '_lowmass'
if xmin > -np.inf    : fileout += '_xmin'+str(xmin)
if xmax <  np.inf    : fileout += '_xmax'+str(xmax)
fileout += '_hist2d.png'
fig.savefig(fileout, dpi=300, bbox_inches='tight')

# 1D hists
s_box = 3000. # Mpc
cenz  = 7500. # Mpc
try:
    start = cat_file.find('_nt')+3
    end = cat_file[start:].find('_')
    ntile = int( cat_file[start:start+end] )
except:
    ntile = 14

ybins = np.where( (x>0) & (x<s_box/ntile) )
xbins = np.where( (y>cenz-s_box/2+5*s_box/ntile) & (y < cenz-s_box/2+6*s_box/ntile) )
    
fig2,ax2 = plt.subplots()
ax2.hist( x[xbins], bins )
ax2.set_xlabel(r'comoving distance $x$ [Mpc]')
ax2.set_ylabel(r'counts')
fig2.savefig(fileout[:-6]+'x.pdf',bbox_inches='tight')

fig3,ax3 = plt.subplots()
ax3.hist( y[ybins],bins )#,range=[6350.,7850.] )
vline_x = np.arange( cenz-s_box/2+s_box/ntile , cenz+s_box/ntile , s_box/ntile )
for j in vline_x: ax3.axvline( j , 0, 20000 , ls='--')
ax3.set_xlabel(r'comoving distance $z$ [Mpc]')
ax3.set_ylabel(r'counts')
ax3.set_xlim([6350.,7850.])
fig3.savefig(fileout[:-6]+'z.pdf',bbox_inches='tight')

