# This script contains functions for reading and plotting Peak Patch halo catalogues.
# 
# USAGE:
# 
# $ python
# >>> import sys.path.insert( 0, '<...>/peakpatch/python/halo_plotting' )
# >>> import pkp_halo_tools as ht

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os,sys
plt.rcParams['text.usetex']=True

class HaloCatalogue:

    def __init__(self, catalogue_file, mass_cut=None, mass_cutoff=1e11, 
            lim=[ None , np.array( [[-np.inf,np.inf]]*3, dtype=np.float32 ) ]
            ):
        
        # Open Peak Patch halo catalogue file
        c_in    = open( catalogue_file , 'rb' )

        # Read halo catalogue header
        self.N_halos  = np.fromfile( c_in, dtype=np.int32  , count=1 )[0]
        self.R_th_max = np.fromfile( c_in, dtype=np.float32, count=1 )[0]
        self.z_obs    = np.fromfile( c_in, dtype=np.float32, count=1 )[0]

        # Read number of columns in catalogue
        self.cols  = int( ( os.path.getsize(catalogue_file)-12)/4/N_halos )

        # Read in Peak Patch halo catalogue (dimension cols x N_halos)
        c = np.reshape( np.fromfile( c_in, dtype=np.float32,
                        count=N_halos*cols ), (N_halos,cols) )

        # Throw out halos outside Lagrangian domain
        if lim[0]:
            if   lim[0][0].lower()[0]=='l': q = 0
            elif lim[0][0].lower()[0]=='e': q = 7
            for j in range(3):
                if lim[1][j,0] > -np.inf: c=c[ np.where( c[:,j+q] > lim[1][j,0] )[0] , : ]
                if lim[1][j,1] < +np.inf: c=c[ np.where( c[:,j+q] > lim[1][j,1] )[0] , : ]

        # NEED TO LOAD PARAMS TO GET A MASS FROM RTH
        # # Filter masses
        # if mass_cut.lower()[0] == 'h': c=c[ np.where( c[:,

        # Allocate halo arrays
        if cols >=  3: self.x ,self.y ,self.z  = np.ndarray.view( c[:,  : 3].T )
        if cols >=  6: self.dx,self.dy,self.dz = np.ndarray.view( c[:, 3: 6].T )
        if cols >=  7: self.R_th               = np.ndarray.view( c[:,  6  ]   )
        if cols >= 10: self.xL,self.yL,self.zL = np.ndarray.view( c[:, 7:10].T )
        if cols >= 11: self.f_pk               = np.ndarray.view( c[:,  10 ]   )
        if cols >= 13: self.e_v,self.p_v       = np.ndarray.view( c[:,11:13].T )
        if cols >= 19:(self.epsxx,self.epsyy,self.epszz,
                       self.epsyz,self.epsxz,self.epsxy) = np.ndarray.view( c[:,13:19].T )
        if cols >= 32:(self.F_d2    , self.zform ,
                       self.gradx   , self.grady , self.gradz,
                       self.gradfx  , self.gradfy, self.gradfz , self.Rfclv  ,
                       self.FcollvRf, self.F_d2Rf, self.gradrfx, self.gradrfy,
                       self.gradrfz ) = np.ndarray.vew( c[:,19:33]

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

